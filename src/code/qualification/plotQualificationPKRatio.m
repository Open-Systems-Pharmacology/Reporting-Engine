function [fig_handle, PKRatioTable, GMFE] = plotQualificationPKRatio(WSettings,figureHandle,PKParameter,PKRatioPlot,ObservedDataSets, SimulationMappings, AxesOptions, PlotSettings, REInputPath)
%PLOTQUALIFICATIONPKRATIO Plots PK ratio from qualification workflow
%
% plotQualificationPKRatio(WSettings,figureHandle,PKParameter,PKRatios,ObservedDataSets, SimulationMappings, AxesOptions, PlotSettings)
%
% Inputs:
%       WSettings (structure)    definition of properties used in all
%                   workflow functions see GETDEFAULTWORKFLOWSETTINGS
%   figureHandle ( handle) handle of figure
%   PKParameter (string) name of the PK parameter to be evaluated
%   PKRatios (struct) with direction to plot the results
%   ObservedDataSets (structure) of all loaded observed data
%   SimulationMappings (structure) to map simulations with observations
%   AxesOptions (structure) to set axes options
%   PlotSettings (structure) to set plot options
% Output
%   csv (cellarray) table with numeric information to the plot

% Open Systems Pharmacology Suite;  http://open-systems-pharmacology.org

%---------------------------------------------

% Loop on the Ratios to be plotted by PK Ratio plot
for i=1:length(PKRatioPlot.PKRatios)
    
    PKRatio=PKRatioPlot.PKRatios(i);
    
    % Get the observed data
    for j=1:length(ObservedDataSets)
        if strcmp(PKRatio.ObservedData, ObservedDataSets(j).Id)
            ObservedData = ObservedDataSets(j).y;
            break
        end
    end
    if isempty(ObservedData.ID(ObservedData.ID==PKRatio.ObservedDataRecordId))
        ME = MException('plotQualificationPKRatio:notFoundInPath', ...
            'In PK Ratio plot %d, Ratio %d, Study ID "%d" was not found in Observed Dataset', figureHandle, i, PKRatio.ObservedDataRecordId);
        throw(ME);
    end
    
    Study = table2array(ObservedData(ObservedData.ID==PKRatio.ObservedDataRecordId,'Study'));
    Result.Study{i,1} = Study{1};
    
    % Load the mapped Time Profile Simulation Results
    [csvSimFile, xmlfile] = getSimFile(PKRatio, SimulationMappings, REInputPath);
    if isempty(csvSimFile)
        ME = MException('plotQualificationPKRatio:notFoundInPath', ...
            'In PK Ratio plot %d, Ratio %d, Project "%s" or Simulation "%s" was not found in SimulationMappings', figureHandle, i, PKRatio.Project, PKRatio.Simulation);
        throw(ME);
    end
    SimResult = loadSimResultcsv(csvSimFile, PKRatio.Simulation);
    
    % All the output are kept so far, may be removed if not necessary
    [AGE, BW, MW, drugmass] = getInfofromSimulation(xmlfile, PKRatio.Output);
    
    Result.AGE(i, 1)=AGE;
    Result.BW(i, 1)=BW;
    Result.MW(i, 1)=MW;
    Result.drugmass(i, 1)=drugmass(1);
    
    % Get the right PK Output
    for j=1:length(SimResult.outputPathList)
        findPathOutput = strfind(SimResult.outputPathList{j}, PKRatio.Output);
        if ~isempty(findPathOutput)
            Pred=SimResult.y{j};
            PredUnit=SimResult.outputUnit{j};
            SimTime=SimResult.time;
            SimTimeUnit=SimResult.timeUnit;
            
            PredDimension = findDimensionfromUnit(PredUnit);
            break
        end
    end
    if isempty(findPathOutput)
        ME = MException('plotQualificationPKRatio:notFoundInPath', ...
            'In PK Ratio plot %d, Ratio %d, Output "%s" was not found in Project "%s" or Simulation "%s"', figureHandle, i, PKRatio.Output, PKRatio.Project, PKRatio.Simulation);
        throw(ME);
    end
    
    % Get the observation time and PK profile in the right unit for study
    ObsTime = ObservedData.Time(ObservedData.ID==PKRatio.ObservedDataRecordId);
    ObsTimeUnit = ObservedData.TimeUnit(ObservedData.ID==PKRatio.ObservedDataRecordId);
    ObsTimeUnit = ObsTimeUnit(1);
    Obs = ObservedData.Avg(ObservedData.ID==PKRatio.ObservedDataRecordId);
    if iscell(Obs)
        Obs = str2double(Obs);
    end
    ObsUnit = ObservedData.AvgUnit(ObservedData.ID==PKRatio.ObservedDataRecordId);
    ObsUnit = ObsUnit{1};
    
    % Get comparable units to pred
    Timefactor=getUnitFactor(ObsTimeUnit,SimTimeUnit,'time');
    % Because of encoding issues, some time ? is read as ? or ?
    try
        Obsfactor=getUnitFactor(ObsUnit,PredUnit,PredDimension, 'MW', MW);
    catch
        ObsUnit = convertPKSimUnit(ObsUnit);
        Obsfactor=getUnitFactor(ObsUnit,PredUnit,PredDimension, 'MW', MW);
    end
    ObsTime = ObsTime.*Timefactor;
    Obs = Obs.*Obsfactor;
    
    % Get PK parameters from the curves
    % For simulation
    allPKpred=getPKParametersForConcentration(SimTime, Pred, 'Dose', drugmass);
    
    % For Observations
    allPKobs=getPKParametersForConcentration(ObsTime, Obs, 'Dose', drugmass);
    
    for k=1:length(PKParameter)
        % Get the PK parameters requested in PKParameter
        if strcmpi(PKParameter{k}, 'AUC')
            PKpredField{k}= 'AUC_last';
        elseif strcmpi(PKParameter{k}, 'CMAX')
            PKpredField{k}= 'cMax';
        else
            PKpredField{k}=PKParameter{k};
        end
        if isfield(allPKpred, PKpredField{k}) && isfield(allPKobs, PKpredField{k})
            
            % Get the observation PK
            Result.obsPK(i, k) = getfield(allPKobs, PKpredField{k});
            % Get the predicted PK
            Result.predPK(i, k) = getfield(allPKpred, PKpredField{k});
            % Get the Ration
            Result.RatioPK(i, k) = Result.predPK(i, k)./Result.obsPK(i, k);
            
        else
            ME = MException('PKRatio:notFoundInField', ...
                'Requested PK Parameter "%s" not found in parameters extracted using getPKParametersForConcentration', PKParameter{k});
            throw(ME);
        end
    end
end

%--------------------------------------
% Plot Section
% Get X parameter:
for k=1:length(AxesOptions)
    if strcmp(AxesOptions(k).Type, 'X')
        Xparam = AxesOptions(k).Dimension;
        % Field is saved as AGE, similar process can be
        % implmented for different X parameters
        if strcmpi(Xparam, 'AGE')
            Xparam = 'AGE';
        end
        break
    end
end
Xvalues = getfield(Result, Xparam);

Xrange=[0.8*min(Xvalues) 1.2*max(Xvalues)]; Yrange=[1 1];

for k=1:length(PKParameter)
    
    % create figure for Obs vs Pred
    [ax, fig_handle(k).PKRatio] = getReportFigureQP(WSettings,1,1,k+figureHandle,PlotSettings);
    setFigureOptions(AxesOptions);
    % Ratio limits
    figure(fig_handle(k).PKRatio);
    plot(Xrange, Yrange, '-k', 'LineWidth', 1, 'HandleVisibility','off');
    plot(Xrange, Yrange/2, '--r', 'LineWidth', 1, 'HandleVisibility','off');
    plot(Xrange, Yrange*2, '--r', 'LineWidth', 1, 'HandleVisibility','off');
    plot(Xrange, Yrange/1.5, '--b', 'LineWidth', 1, 'HandleVisibility','off');
    plot(Xrange, Yrange*1.5, '--b', 'LineWidth', 1, 'HandleVisibility','off');
    
    pp=plot(Xvalues, Result.RatioPK(:, k), 'o', 'Linewidth',1);
    setCurveOptions(pp, PKRatioPlot);
    
    ylabel(sprintf('Predicted %s / Observed %s', PKParameter{k}, PKParameter{k}));
    axis([Xrange 0.8*min(min(Result.RatioPK(k,:)), 0.5) 1.2*max(max(Result.RatioPK(k,:)), 2)]);
    
    legend('off');
    
end

% --------------------------------------
% Table and Qualification Section

% Get the DDI Ratio Table
PKRatioHeader={};
PKRatioResults=[];

for k=1:length(PKParameter)
    PKRatioHeader = [PKRatioHeader {sprintf('Predicted %s', PKParameter{k}), ...
        sprintf('Observed %s', PKParameter{k}), sprintf('Pred/Obs %s Ratio', PKParameter{k})}];
    
    PKRatioResults = [PKRatioResults Result.predPK(:,k), Result.obsPK(:,k), Result.RatioPK(:,k)];
    
end

% Definition of the PKRatio table
PKRatioHeader = [{'Study ID', 'Age (y)', 'BodyWeight (kg)'}, PKRatioHeader];

PKRatioTable = [PKRatioHeader ; ...
    Result.Study num2cell([Result.AGE, Result.BW, PKRatioResults])];

disp(PKRatioTable);

% Calculation of GMFE
GMFE = 10.^(sum(abs(log(Result.RatioPK)))./length(Result.obsPK));
for k=1:length(PKParameter)
    fprintf('%s: GMFE = %f \n', PKParameter{k}, GMFE(k));
end

function [AGE, BW, MW, drugmass] = getInfofromSimulation(xmlfile, Output)

initSimulation(xmlfile,'none');

MW = getMolecularWeightForPath(Output);

drugmass = getParameter('*Application_*|ProtocolSchemaItem|DrugMass',1,'parametertype','readonly');

AGE = getParameter('*|Organism|Age',1,'parametertype','readonly');
BW = getParameter('*|Organism|Weight',1,'parametertype','readonly');

function UnitOut = convertPKSimUnit(UnitIn)

UnitOut = UnitIn;
UnitOut(UnitIn=='L')='l';
UnitOut(1)='µ';
