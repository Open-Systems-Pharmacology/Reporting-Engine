function [PKRatioTable, GMFE] = plotQualificationPKRatio(WSettings,figureHandle,PKParameter,PKRatios,ObservedDataSets, SimulationMappings, AxesOptions, PlotSettings, CurveOptions)
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
% Create figure with first setting from WSettings using getReportFigure
% To be updated using the Configuration plan Settings as optional arguments
[ax, fig_handle1] = getReportFigureQP(WSettings,1,1,figureHandle,PlotSettings);

[xAxesOptions, yAxesOptions] = setFigureOptions(AxesOptions);

% Initialize error for computing GMFE
Error=[];

% Get dimensions for scaling plots
dimensionList=getDimensions;

% Loop on the Ratios to be plotted by PK Ratio plot
for i=1:length(PKRatios)
    
    PKRatio=PKRatios(i);
    
    % Get the observed data
    for j=1:length(ObservedDataSets)
        if strcmp(PKRatio.ObservedData, ObservedDataSets(j).Id)
            ObservedData = ObservedDataSets(j).y;
            break
        end
    end
    Study(i)=table2array(ObservedData(ObservedData.ID==PKRatio.ObservedDataRecordId,'Study'));
    
    % Load the mapped Time Profile Simulation Results
    [csvSimFile, xmlfile] = getSimFile(PKRatio, SimulationMappings);
    SimResult = loadSimResultcsv(csvSimFile, PKRatio.Simulation);
    
    initSimulation(xmlfile,'none');
    drugmass = getParameter('*Application_*|ProtocolSchemaItem|DrugMass',1,'parametertype','readonly');
    Compound{i}=PKRatio.Project;
    AGE(i) = getParameter('*|Organism|Age',1,'parametertype','readonly');
    BW(i) = getParameter('*|Organism|Weight',1,'parametertype','readonly');
    MW(i) = getParameter(sprintf('*|%s|Molecular weight',Compound{i}),1,'parametertype','readonly');
    
    X(i) = getParameter(sprintf('*|Organism|%s', xAxesOptions.Dimension),1,'parametertype','readonly');
    
    % Use the right Output
    for j=1:length(SimResult.outputPathList)
        findPathOutput = strfind(SimResult.outputPathList{j}, PKRatio.Output);
        if ~isempty(findPathOutput)
            predicted=SimResult.y{j};
            predictedTime=SimResult.time;
            predictedTimeUnit=SimResult.timeUnit;
            break
        end
    end
    
    % Get the PK parameter out of the simulation
    allPKpred=getPKParametersForConcentration(predictedTime,predicted/1000,'Dose', drugmass);
    
    XDimension=dimensionList{strContains(dimensionList, xAxesOptions.Dimension)};
    if strContains({xAxesOptions.Dimension}, 'Age')
        XUnit = getUnitsForDimension('Age in years');
        XUnit = XUnit{1};
    end
    if strContains({xAxesOptions.Dimension}, 'BodyWeight')
        XUnit = getUnitsForDimension('Mass');
        XUnit = XUnit{1};
    end
    
    if isfield(allPKpred, PKParameter)
        
        % ISSUE TO BE FIXED
        % PK obs unit is mL/min/kg whereas PK seems to be total clearance in L/min ?
        PKobs(i) = BW(i)*table2array(ObservedData(ObservedData.ID==PKRatio.ObservedDataRecordId,...
            (strcmp(ObservedData.Properties.VariableNames, [PKParameter 'Avg']))));
        
        PKobsUnit(i) = table2array(ObservedData(ObservedData.ID==PKRatio.ObservedDataRecordId,...
            (strcmp(ObservedData.Properties.VariableNames, [PKParameter 'AvgUnit']))));
        
        % PKDimension=dimensionList{strContains(PKParameter, dimensionList)};
        % PKfactor =getUnitFactor(PKobsUnit(i),stdUnit,PKDimension, 'BW',BW(i), 'MW', MW(i));
        % PKobs(i) = PKobs(i)*PKfactor;
        
        Xfactor=getUnitFactor(XUnit,xAxesOptions.Unit,XDimension);
        
        PKpred(i)=getfield(allPKpred, PKParameter);
        pp=plot(X(i).*Xfactor, PKpred(i)/PKobs(i), 'Linewidth',1);
        setCurveOptions(pp, CurveOptions);
        
    else
        ME = MException('PKRatio:notFoundInPath', ...
            'PKParameter %s not found using getPKParametersForConcentration', PKParameter);
        throw(ME);
    end
    
end

% Definition of the PKRatio table
PKRatioHeader = {'StudyID', 'Age (y)', 'BodyWeight (kg)', ['Report' PKParameter], ['PKSim' PKParameter], 'Ratio'};
PKRatioTable = [PKRatioHeader ; ...
    Study' num2cell([AGE', BW', PKobs', PKpred', PKpred'./PKobs'])];

disp(PKRatioTable);

% Calculation of GMFE
GMFE = 10.^(sum(abs(log(PKpred./PKobs)))/length(PKobs));
fprintf('GMFE = %f \n', GMFE);

% Perform the plot based on Curves indications
Xrange=[0.8*min(AGE) 1.2*max(AGE)]; Yrange=[1 1];

plot(Xrange, Yrange, '-k', 'LineWidth', 1);
plot(Xrange, Yrange/2, '--r', 'LineWidth', 1); plot(Xrange, Yrange*2, '--r', 'LineWidth', 1);
plot(Xrange, Yrange/1.5, '--b', 'LineWidth', 1); plot(Xrange, Yrange*1.5, '--b', 'LineWidth', 1);

ylim([min([0.1 PKpred./PKobs]) max([10 PKpred./PKobs])]);
ylabel(sprintf('PK ratio %s Pred / %s Obs', PKParameter, PKParameter));
% Placement of legend without masking the plot has not been figure out
% Update to be perform
legend('off')