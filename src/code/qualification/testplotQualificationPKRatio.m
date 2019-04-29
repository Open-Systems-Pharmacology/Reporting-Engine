function testplotQualificationPKRatio(WSettings,figureHandle,PKParameter,PKRatios,ObservedDataSets, SimulationMappings, AxesOptions, PlotSettings)
%PLOTQUALIFICATIONTIMEPROFILE Plots the time profile of a population in comparison to a reference population
%
% plotQualificationTimeProfile(WSettings,figureHandle,SimTL,DataTP, Curves, AxesOptions,PlotSettings)
%
% Inputs:
%       WSettings (structure)    definition of properties used in all
%                   workflow functions see GETDEFAULTWORKFLOWSETTINGS
%   figureHandle ( handle) handle of figure
%   SimTL (structure) with simulated results
%       .time (double vector): time vector (in units of time_unit)
%       .y  (double matrix (nIndividual x n timepoints)): time profile of population
%       .timeRef (double vector): time vector of reference population (in units of time_unit)
%           if empty no reference population is drawn
%       .yRef  (double matrix (nIndividual x n timepoints)): time profile of of reference population
%           if empty no reference population is drawn
%       .meanModel  (double matrix (1 x n timepoints)): time profile of population
%   DataTP (structure) dataprofile
%
% Output
%   csv (cellarray) table with numeric information to the plot

% Open Systems Pharmacology Suite;  http://open-systems-pharmacology.org

%---------------------------------------------
% Create figure with first setting from WSettings using getReportFigure
% To be updated using the Configuration plan Settings as optional arguments
ax = getReportFigureQP(WSettings,1,1,figureHandle,PlotSettings);

setFigureOptions(AxesOptions);

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
    
    % TO BE FIXED
    % PK obs unit is L/h/kg whereas PK seems to be total clearance in L/min ?
    BW(i)=10;
    AGE(i)=1;
    
    PKobs = BW(i)*table2array(ObservedData(ObservedData.ID==PKRatio.ObservedDataRecordId,...
        (strcmp(ObservedData.Properties.VariableNames, [PKParameter 'Avg']))));
    
    % Load the mapped Time Profile Simulation Results
    [csvSimFile, xmlfile] = getSimFile(PKRatio, SimulationMappings);
    SimResult = loadSimResultcsv(csvSimFile, PKRatio.Simulation);
    
    initSimulation(xmlfile,'none');
    drugmass = getParameter('*Application_*|ProtocolSchemaItem|DrugMass',1,'parametertype','readonly');
    
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
    
    allPKpred=getPKParametersForConcentration(predictedTime/60,predicted,'Dose', drugmass);
    
    if isfield(allPKpred, PKParameter)
        PKpred=getfield(allPKpred, PKParameter);
        pp=plot(AGE(i), PKpred/PKobs, 'ok', 'Linewidth',1);
    else
        disp('Error')
    end
    
end

% Perform the plot based on Curves indications
Xrange=[0 5]; Yrange=[1 1];
plot(Xrange, Yrange, '-k', 'LineWidth', 1);
plot(Xrange, Yrange/2, '--r', 'LineWidth', 1); plot(Xrange, Yrange*2, '--r', 'LineWidth', 1);
plot(Xrange, Yrange/1.5, '--b', 'LineWidth', 1); plot(Xrange, Yrange*1.5, '--b', 'LineWidth', 1);

legend('off')