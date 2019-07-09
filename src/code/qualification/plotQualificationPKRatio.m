function [fig_handle, PKRatioTable, PKRatioQuali, GMFE] = plotQualificationPKRatio(WSettings,plotIndex,PKParameter,PKRatioPlot,ObservedDataSets, SimulationMappings, AxesOptions, PlotSettings, REInputPath)
%PLOTQUALIFICATIONPKRATIO Plots PK Ratios from Configuration Plan
%
% [fig_handle, DDIRatioTable, DDIRatioQuali] = plotQualificationDDIRatio(WSettings,plotIndex,
%   PKParameter,DDIRatioGroups,ObservedDataSets, SimulationMappings, AxesOptions, PlotSettings, REInputPath)
%
% Inputs:
%   WSettings (structure)    definition of properties used in all
%                   workflow functions see GETDEFAULTWORKFLOWSETTINGS
%   plotIndex (integer) index of plot
%   PKParameter (cells) name of the PK parameter to be evaluated
%   PKRatioPlot (structure) PK Ratio plot information
%   ObservedDataSets (structure) Observed data
%   SimulationMappings (structure) Map simulation results to project
%   AxesOptions (structure) to set axes options
%   PlotSettings (structure) to set plot options
%   REInputPath (string) path of RE input files and folders
% Output
%   fig_handle (handle) handle of output figures
%   PKRatioTable (cells) Table of DDI Ratio information
%   GMFE (cells)  Geometric Mean Fold Error
%

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
    
    % Get the matching Record ID in the Table
    ID = (ObservedData.ID==PKRatio.ObservedDataRecordId);
    
    if max(ID)==0
        ME = MException('plotQualificationPKRatio:notFoundInPath', ...
            'In PK Ratio Plot %d, Ratio %d, Study ID "%d" was not found in Observed Dataset', plotIndex, i, PKRatio.ObservedDataRecordId);
        throw(ME);
    end
    
    % Get the Study
    Result.Study(i,1) = ObservedData.Study(ID);
    
    % Get the requested PK Parameters, their unit and dimension
    for k=1:length(PKParameter)
        
        Result.obsPK(i, k) = table2array(ObservedData(ID, sprintf('%sAvg', PKParameter{k})));
        obsPKUnit = table2array(ObservedData(ID, sprintf('%sAvgUnit', PKParameter{k})));
        Result.obsPKUnit{i, k} = convertPKSimUnit(obsPKUnit);
        Result.obsPKDimension{i, k} = findDimensionfromUnit(Result.obsPKUnit{i, k});
        
        if isempty(Result.obsPKDimension{i, k})
            ME = MException('plotQualificationPKRatio:unknownUnit', ...
                'In PK Ratio Plot %d, Ratio %d, Dimension unknown for Unit "%s" in Observed Data Record ID "%d" \n', plotIndex, i, char(Result.obsPKUnit{i, k}), PKRatio.ObservedDataRecordId);
            throw(ME);
        end
        
    end
    
    % Load the mapped Time Profile Simulation Results
    [csvSimFile, xmlfile] = getSimFile(PKRatio, SimulationMappings, REInputPath);
    if isempty(csvSimFile)
        ME = MException('plotQualificationPKRatio:notFoundInPath', ...
            'In PK Ratio Plot %d, Ratio %d, Project "%s" or Simulation "%s" was not found in SimulationMappings', plotIndex, i, PKRatio.Project, PKRatio.Simulation);
        throw(ME);
    end
    SimResult = loadSimResultcsv(csvSimFile, PKRatio.Simulation);
    
    % All the output are kept so far, may be removed if not necessary
    [AGE, BW, MW, drugmass, drugmassUnit] = getInfofromSimulation(xmlfile, PKRatio.Output);
    
    Result.AGE(i, 1)=AGE;
    Result.BW(i, 1)=BW;
    Result.MW(i, 1)=MW;
    Result.drugmass(i, 1)=drugmass(1);
    Result.drugmassUnit{i, 1}=drugmassUnit(1);
    
    % Get the right PK Output
    if isempty(SimResult.outputPathList)
        ME = MException('plotQualificationPKRatio:emptyOutputPathInSimulation', ...
            'In PK Ratio Plot %d, Ratio %d, OutputPath is empty in Project "%s" or Simulation "%s"', plotIndex, i, PKRatio.Project, PKRatio.Simulation);
        throw(ME);
    end
    for j=1:length(SimResult.outputPathList)
        findPathOutput = strfind(SimResult.outputPathList{j}, PKRatio.Output);
        if ~isempty(findPathOutput)
            % Get Time and Concentration in PK Sim internal units
            % Concentration in µmol/l and
            % Time in min
            SimTime=SimResult.time;
            Pred=SimResult.y{j};
            
            break
        end
    end
    if isempty(findPathOutput)
        ME = MException('plotQualificationPKRatio:notFoundInPath', ...
            'In PK Ratio Plot %d, Ratio %d, Output "%s" was not found in Project "%s" or Simulation "%s"', plotIndex, i, PKRatio.Output, PKRatio.Project, PKRatio.Simulation);
        throw(ME);
    end
    
    % Get PK parameters from the curves
    % For simulation
    allPKpred=getPKParametersForConcentration(SimTime, Pred, 'Dose', drugmass);
    
    for k=1:length(PKParameter)
        % Get the PK parameters requested in PK Parameters
        % Internal Units are assumed for PK parameters
        % according to Obs Unit
        if strcmp(Result.obsPKDimension(i, k), 'AUC (mass)')
            % Internal Unit for AUC is µmol*min/l and MW in g/mol
            AUCpred = getfield(allPKpred, 'AUC_last')*MW;
            AUCpredUnitFactor = getUnitFactor('µg*min/l', Result.obsPKUnit{i, k}, 'AUC (mass)');
            Result.predPK(i, k) = AUCpred.*AUCpredUnitFactor;
            
        elseif strcmp(Result.obsPKDimension(i, k), 'AUC (molar)')
            % Internal Unit for AUC is µmol*min/l
            AUCpred = getfield(allPKpred, 'AUC_last');
            AUCpredUnitFactor = getUnitFactor('µmol*min/l', Result.obsPKUnit{i, k}, 'AUC (molar)');
            Result.predPK(i, k) = AUCpred.*AUCpredUnitFactor;
            
        elseif strcmp(Result.obsPKDimension(i, k), 'Concentration')
            % Internal Unit for Cmax is µmol/l
            CMAXpred = getfield(allPKpred, 'cMax');
            CMAXpredUnitFactor = getUnitFactor('µmol/l', Result.obsPKUnit{i, k}, 'Concentration', 'MW', MW);
            Result.predPK(i, k) = CMAXpred.*CMAXpredUnitFactor;
            
        elseif strcmp(Result.obsPKDimension(i, k), 'Flow')
            % Internal Unit for CL is l/min
            CLpred = getfield(allPKpred, 'CL');
            CLpredUnitFactor = getUnitFactor('l/min', Result.obsPKUnit{i, k}, 'Flow');
            Result.predPK(i, k) = CLpred.*CLpredUnitFactor;
            
        elseif strcmp(Result.obsPKDimension(i, k), 'Flow per weight')
            % Internal Units for CL is l/min and for BW is kg
            CLpred = getfield(allPKpred, 'CL')/BW;
            CLpredUnitFactor = getUnitFactor('l/min/kg', Result.obsPKUnit{i, k}, 'Flow per weight');
            Result.predPK(i, k) = CLpred.*CLpredUnitFactor;
            
        else
            ME = MException('plotQualificationPKRatio:unknownDimension', ...
                'In PK Ratio Plot %d, Ratio %d, Observed Study ID "%d", PK Parameter "%s" \n Unknown dimension for observed unit "%s" ', plotIndex, i, PKRatio.ObservedDataRecordId, PKParameter{k}, char(Result.obsPKUnit{i, k}));
            throw(ME);
        end
        
        % Get the Ratio
        Result.RatioPK(i, k) = Result.predPK(i, k)./Result.obsPK(i, k);
        
    end
end

%--------------------------------------
% Plot Section
% Check that X parameter is AGE
for k=1:length(AxesOptions)
    if strcmp(AxesOptions(k).Type, 'X')
        Xparam = AxesOptions(k).Dimension;
        % Field is saved as AGE, similar process can be
        % implmented for different X parameters
        if strcmpi(Xparam, 'AGE')
            Xparam = 'AGE';
        else
            ME = MException('plotQualificationPKRatio:XaxisDimension', ...
                'In PK Ratio Plot %d,, X-axis does not correspond to Age', plotIndex);
            throw(ME);
        end
        break
    end
end
Xvalues = getfield(Result, Xparam);

Xrange=[0.8*min(Xvalues) 1.2*max(Xvalues)]; Yrange=[1 1];

for k=1:length(PKParameter)
    
    % create figure for Obs vs Pred
    [ax, fig_handle(k).PKRatio] = getReportFigureQP(WSettings,1,1,[],PlotSettings);
    setFigureOptions(AxesOptions);
    % Ratio limits
    plot(Xrange, Yrange, '-k', 'LineWidth', 1, 'HandleVisibility','off');
    plot(Xrange, Yrange/2, '--r', 'LineWidth', 1, 'HandleVisibility','off');
    plot(Xrange, Yrange*2, '--r', 'LineWidth', 1, 'HandleVisibility','off');
    plot(Xrange, Yrange/1.5, '--b', 'LineWidth', 1, 'HandleVisibility','off');
    plot(Xrange, Yrange*1.5, '--b', 'LineWidth', 1, 'HandleVisibility','off');
    
    pp=plot(Xvalues, Result.RatioPK(:, k), 'o', 'Linewidth',1);
    setCurveOptions(pp, PKRatioPlot);
    
    ylabel(sprintf('Predicted %s / Observed %s', PKParameter{k}, PKParameter{k}));
    axis([Xrange 0.8*min(min(Result.RatioPK(:,k)), 0.5) 1.2*max(max(Result.RatioPK(:,k)), 2)]);
    
    legend('off');
    
end

% --------------------------------------
% Table and Qualification Section

% Calculation of GMFE
GMFE = 10.^(sum(abs(log10(Result.RatioPK)))./length(Result.obsPK));
for k=1:length(PKParameter)
    fprintf('%s: GMFE = %f \n', PKParameter{k}, GMFE(k));
end

% Get the PK Ratio Table
PKRatioHeader={};
PKRatioResults=[];

for k=1:length(PKParameter)
    PKRatioHeader = [PKRatioHeader {sprintf('Predicted %s [%s]', PKParameter{k}, Result.obsPKUnit{1, k}), ...
        sprintf('Observed %s [%s]', PKParameter{k}, Result.obsPKUnit{1, k}), sprintf('Pred/Obs %s Ratio', PKParameter{k})}];
    
    PKRatioResults = [PKRatioResults Result.predPK(:,k), Result.obsPK(:,k), Result.RatioPK(:,k)];
    
end

% Definition of the PKRatio table
PKRatioHeader = [{'Study ID', 'Age [y]', 'BodyWeight [kg]'}, PKRatioHeader];

PKRatioTable = [PKRatioHeader ; ...
    Result.Study num2cell([Result.AGE, Result.BW, PKRatioResults])];

disp(PKRatioTable);

% Get the PK Ratio Qualification
for k=1:length(PKParameter)
    % Update the number of points for each group
    conditionPoints = ~isnan(Result.RatioPK(:, k));
    QualiMeasure(k).PointsTotal = length(Result.RatioPK(conditionPoints,k));
    condition15fold = (Result.RatioPK(:,k) >= 1/1.5 & Result.RatioPK(:,k) <= 1.5);
    QualiMeasure(k).Points15fold= length(Result.RatioPK(condition15fold,k));
    condition2fold = (Result.RatioPK(:,k) >= 0.5 & Result.RatioPK(:,k) <= 2);
    QualiMeasure(k).Points2fold= length(Result.RatioPK(condition2fold,k));
    
    PKRatioQualiHeader = {PKParameter{k}, 'Number', 'Ratio [%]'};
    PKRatioQuali_1st_Column = {'Points total'; 'Points within 1.5 fold'; 'Points within 2-fold'};
    PKRatioQuali_Other_Columns = num2cell([QualiMeasure(k).PointsTotal NaN; ...
        QualiMeasure(k).Points15fold 100.*QualiMeasure(k).Points15fold./QualiMeasure(k).PointsTotal; ...
        QualiMeasure(k).Points2fold 100.*QualiMeasure(k).Points2fold./QualiMeasure(k).PointsTotal]);
    PKRatioQuali_Other_Columns{1,2}='-';
    
    
    PKRatioQuali(k).Output = [PKRatioQualiHeader  ; ...
        PKRatioQuali_1st_Column PKRatioQuali_Other_Columns];
    
    fprintf('Qualification Measures for PK parameter %s \n', PKParameter{k});
    disp(PKRatioQuali(k).Output);
    
end

function [AGE, BW, MW, drugmass, drugmassUnit] = getInfofromSimulation(xmlfile, Output)

initSimulation(xmlfile,'none');

MW = getMolecularWeightForPath(Output);

drugmass = getParameter('*Application_*|ProtocolSchemaItem|DrugMass',1,'parametertype','readonly');
drugmassUnit = getParameter('*Application_*|ProtocolSchemaItem|DrugMass',1,'parametertype','readonly', 'property', 'Unit');

AGE = getParameter('*|Organism|Age',1,'parametertype','readonly');
BW = getParameter('*|Organism|Weight',1,'parametertype','readonly');

function UnitOut = convertPKSimUnit(UnitIn)

UnitOut = string(UnitIn);
UnitOut = replace(UnitOut,'L','l');

%UnitOut(1)='µ';
