function plotQualificationTimeProfile(WSettings,figureHandle,TimeProfile,ObservedDataSets,SimulationMappings, Curves, AxesOptions, PlotSettings)
%PLOTQUALIFICATIONTIMEPROFILE Plots the time profile of a population in comparison to a reference population
%
% plotQualificationTimeProfile(WSettings,figureHandle,SimTL,DataTP, Curves, AxesOptions,PlotSettings)
%
% Inputs:
%       WSettings (structure)    definition of properties used in all
%                   workflow functions see GETDEFAULTWORKFLOWSETTINGS
%   figureHandle ( handle) handle of figure
%   SimTL (structure) with simulated results
%   DataTP (structure) of all loaded observed data
%   Curves (structure) to map simulations with observations
%   AxesOptions (structure) to set plot options
%   PlotSettings (structure) to set plot options
%
% Open Systems Pharmacology Suite;  http://open-systems-pharmacology.org

%---------------------------------------------

% Create figure with first setting from WSettings using getReportFigure
% To be updated using the Configuration plan Settings as optional arguments

ax = getReportFigureQP(WSettings,1,1,figureHandle,PlotSettings);

[xAxesOptions, yAxesOptions, yyAxesOptions] = setFigureOptions(AxesOptions);

% Load the mapped Time Profile Simulation Results
[csvSimFile, xmlfile] = getSimFile(TimeProfile, SimulationMappings);
SimResult = loadSimResultcsv(csvSimFile, TimeProfile);

initSimulation(xmlfile,'none');
Compound=TimeProfile.Project;
MW = getParameter(sprintf('*|%s|Molecular weight',Compound),1,'parametertype','readonly');


% Perform the plot based on Curves indications
legendLabels={};
hold on;

% Get dimensions for scaling plots
dimensionList=getDimensions;

% Plot the Curves indicated by the array of structures Curves
for i=1:length(Curves)
    
    % For simulations
    for j=1:length(SimResult.outputPathList)
        % Get the right simulation curve
        findPathOutput = strfind(Curves(i).Y,SimResult.outputPathList{j});
        if ~isempty(findPathOutput)
            legendLabels{length(legendLabels)+1}=Curves(i).Name;
            
            if isfield(Curves(i).CurveOptions, 'yAxisType')
                if strcmp(Curves(i).CurveOptions.yAxisType, 'Y2')
                    yyaxis right
                    YDimension=dimensionList{strContains(yyAxesOptions.Dimension, dimensionList)};
                    Yfactor=getUnitFactor(SimResult.outputUnit{j},yyAxesOptions.Unit,YDimension);
                else
                    yyaxis left
                end
            else
                YDimension=dimensionList{strContains(yAxesOptions.Dimension, dimensionList)};
                Yfactor=getUnitFactor(SimResult.outputUnit{j},yAxesOptions.Unit,YDimension, 'MW',MW*1e10);
            end
            
            % Convert units to reference unit
            XDimension=dimensionList{strContains(xAxesOptions.Dimension, dimensionList)};
            Xfactor=getUnitFactor(SimResult.timeUnit,xAxesOptions.Unit,XDimension);
            
            pp = plot(SimResult.time.*Xfactor, SimResult.y{j}.*Yfactor);
            break
        end
    end
    % For observations
    for j=1:length(ObservedDataSets)
        % Get the right observation profile
        findPathOutput = strfind(Curves(i).Y,ObservedDataSets(j).Id);
        
        if ~isempty(findPathOutput)
            legendLabels{length(legendLabels)+1}=Curves(i).Name;
            
            if isfield(Curves(i).CurveOptions, 'yAxisType')
                if strcmp(Curves(i).CurveOptions.yAxisType, 'Y2')
                    yyaxis right
                    YDimension=dimensionList{strContains(yyAxesOptions.Dimension, dimensionList)};
                    Yfactor=getUnitFactor(ObservedDataSets(j).outputUnit,yyAxesOptions.Unit,YDimension);
                else
                    yyaxis left
                end
            else
                YDimension=dimensionList{strContains(yAxesOptions.Dimension, dimensionList)};
                % Caution: Code to be updated
                % Observed concentration was massic whereas output was
                % molar
                Yfactor=getUnitFactor(ObservedDataSets(j).outputUnit{1},yAxesOptions.Unit,YDimension, 'MW',MW*1e10);
                
            end
            XDimension=dimensionList{strContains(xAxesOptions.Dimension, dimensionList)};
            Xfactor=getUnitFactor(ObservedDataSets(j).timeUnit,xAxesOptions.Unit,XDimension);
            pp = plot(ObservedDataSets(j).time.*Xfactor, ObservedDataSets(j).y{1}.*Yfactor);
            break
        end
    end
    
    % If the output was not found
    if ~exist('pp')
        ME = MException('TimeProfile:notFoundInPath', ...
            'Curve %s not found within simulation or observation files', Curves(i).Y);
        throw(ME);
    end
    setCurveOptions(pp, Curves(i).CurveOptions);
end
%legend(legendLabels);
legend('off')