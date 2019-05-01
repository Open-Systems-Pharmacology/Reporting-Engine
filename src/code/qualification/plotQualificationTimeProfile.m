function plotQualificationTimeProfile(WSettings,figureHandle,TimeProfile,ObservedDataSets,SimulationMappings, Curves, AxesOptions, PlotSettings)
%PLOTQUALIFICATIONTIMEPROFILE Plots the time profile of a population in comparison to a reference population
%
% plotQualificationTimeProfile(WSettings,figureHandle,SimTL,DataTP, Curves, AxesOptions,PlotSettings)
%
% Inputs:
%       WSettings (structure)    definition of properties used in all
%                   workflow functions see GETDEFAULTWORKFLOWSETTINGS
%   figureHandle ( handle) handle of figure
%   TimeProfile (structure) with Project and Simulation to be mapped
%   ObservedDataSets (structure) of all loaded observed data
%   SimulationMappings (structure) to mapping linking simulations to there path
%   Curves (structure) to internal path of simulations or observations
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

% Initialize simulation, and get Molecular Weight in g/mol for correct use of getUnitFactor
initSimulation(xmlfile,'none');
Compound=TimeProfile.Project;
MW = getParameter(sprintf('*|%s|Molecular weight',Compound),1,'parametertype','readonly');
MWUnit = getParameter(sprintf('*|%s|Molecular weight',Compound),1,'parametertype','readonly', 'property', 'Unit');
MW = MW.*getUnitFactor(MWUnit, 'g/mol', 'Molecular weight');

% Perform the plot based on Curves indications
legendLabels={};
hold on;

% Plot the Curves indicated by the array of structures Curves
for i=1:length(Curves)
    
    % For simulations: Get the right simulation curve
    [p_handle, legLabel] = testandplotSimResults(Curves(i), SimResult, MW, xAxesOptions, yAxesOptions, yyAxesOptions);
    
    if isempty(p_handle)
        % For observation: Get the right observation curve with right unit
        [p_handle, legLabel] = testandplotObservations(Curves(i), ObservedDataSets, MW, xAxesOptions, yAxesOptions, yyAxesOptions);
        
        % If the output was not found
        if isempty(p_handle)
            ME = MException('plotQualificationTimeProfile:notFoundInPath', ...
                'Curves %d : %s not found', i, Curves(i).Y);
            throw(ME);
        end
    end
    legendLabels{length(legendLabels)+1}=legLabel;
    
end
legend(legendLabels, 'Location', 'northoutside');
%legend('off')

% ------------------------- Auxiliary functions -------------------------
% For simulations: Get the right simulation curve with right unit
function [p_handle, legendLabels] = testandplotSimResults(Curves, SimResult, MW, xAxesOptions, yAxesOptions, yyAxesOptions)

dimensionList=getDimensions;

for j = 1:length(SimResult.outputPathList)
    
    findPathOutput = contains(Curves.Y,SimResult.outputPathList{j});
    if findPathOutput
        legendLabels=Curves.Name;
        
        if isfield(Curves.CurveOptions, 'yAxisType')
            if strcmp(Curves.CurveOptions.yAxisType, 'Y2')
                yyaxis right
                YDimension=dimensionList{strContains(yyAxesOptions.Dimension, dimensionList)};
                Yfactor=getUnitFactor(SimResult.outputUnit{j},yyAxesOptions.Unit,YDimension, 'MW',MW);
                
                % Convert units to reference unit
                XDimension=dimensionList{strContains(xAxesOptions.Dimension, dimensionList)};
                Xfactor=getUnitFactor(SimResult.timeUnit,xAxesOptions.Unit,XDimension);
                
                p_handle = plot(SimResult.time.*Xfactor, SimResult.y{j}.*Yfactor);
                yyaxis left
                break
            end
        else
            
            YDimension=dimensionList{strContains(yAxesOptions.Dimension, dimensionList)};
            Yfactor=getUnitFactor(SimResult.outputUnit{j},yAxesOptions.Unit,YDimension, 'MW',MW);
            
            % Convert units to reference unit
            XDimension=dimensionList{strContains(xAxesOptions.Dimension, dimensionList)};
            Xfactor=getUnitFactor(SimResult.timeUnit,xAxesOptions.Unit,XDimension);
            
            p_handle = plot(SimResult.time.*Xfactor, SimResult.y{j}.*Yfactor);
            break
        end
    end
end
if exist('p_handle')
    setCurveOptions(p_handle, Curves.CurveOptions);
else
    p_handle=[];
    legendLabels=[];
end

% For observation: Get the right observation curve with right unit
function [p_handle, legendLabels] = testandplotObservations(Curves, ObservedDataSets, MW, xAxesOptions, yAxesOptions, yyAxesOptions)

dimensionList=getDimensions;

for j = 1:length(ObservedDataSets)
    findPathOutput = contains(Curves.Y,ObservedDataSets(j).Id);
    if findPathOutput
        for k=1:length(ObservedDataSets(j).outputPathList)
            findPathOutputj = contains(Curves.Y,ObservedDataSets(j).outputPathList{k});
            if findPathOutputj
                
                legendLabels=Curves.Name;
                
                if isfield(Curves.CurveOptions, 'yAxisType')
                    if strcmp(Curves.CurveOptions.yAxisType, 'Y2')
                        yyaxis right
                        YDimension=dimensionList{strContains(yyAxesOptions.Dimension, dimensionList)};
                        Yfactor=getUnitFactor(ObservedDataSets(j).outputUnit{j},yyAxesOptions.Unit,YDimension, 'MW',MW);
                        
                        % Convert units to reference unit
                        XDimension=dimensionList{strContains(xAxesOptions.Dimension, dimensionList)};
                        Xfactor=getUnitFactor(ObservedDataSets(j).timeUnit,xAxesOptions.Unit,XDimension);
                        
                        p_handle = plot(ObservedDataSets(j).time.*Xfactor, ObservedDataSets(j).y{k}.*Yfactor);
                        
                        % Check if error bars to be plotted
                        if length(ObservedDataSets(j).outputPathList)>1
                            for kk=1:length(ObservedDataSets(j).outputPathList)
                                errorObs = contains(ObservedDataSets(j).outputPathList{kk}, 'SD');
                                if errorObs
                                    errorfactor=getUnitFactor(ObservedDataSets(j).outputUnit{kk},yyAxesOptions.Unit,YDimension, 'MW',MW);
                                    p_handle2=errorbar(ObservedDataSets(j).time.*Xfactor, ObservedDataSets(j).y{k}.*Yfactor, ObservedDataSets(j).y{kk}.*errorfactor);
                                    setCurveOptions(p_handle2, Curves.CurveOptions);
                                end
                            end
                        end
                        
                        yyaxis left
                        break
                    end
                else
                    
                    YDimension=dimensionList{strContains(yAxesOptions.Dimension, dimensionList)};
                    Yfactor=getUnitFactor(ObservedDataSets(j).outputUnit{k},yAxesOptions.Unit,YDimension, 'MW',MW);
                    
                    %fprintf('%d %d %s %s %s %f \n', i, j, YDimension, SimResult.outputUnit{j}, yAxesOptions.Unit, Yfactor);
                    %fprintf('%s \n %s \n\n', Curves(i).Y, SimResult.outputPathList{j});
                    
                    % Convert units to reference unit
                    XDimension=dimensionList{strContains(xAxesOptions.Dimension, dimensionList)};
                    Xfactor=getUnitFactor(ObservedDataSets(j).timeUnit,xAxesOptions.Unit,XDimension);
                    
                    p_handle = plot(ObservedDataSets(j).time.*Xfactor, ObservedDataSets(j).y{k}.*Yfactor);
                    
                    % Check if error bars to be plotted
                    if length(ObservedDataSets(j).outputPathList)>1
                        for kk=1:length(ObservedDataSets(j).outputPathList)
                            errorObs = contains(ObservedDataSets(j).outputPathList{kk}, 'SD');
                            if errorObs
                                errorfactor=getUnitFactor(ObservedDataSets(j).outputUnit{kk},yAxesOptions.Unit,YDimension, 'MW',MW);
                                p_handle2=errorbar(ObservedDataSets(j).time.*Xfactor, ObservedDataSets(j).y{k}.*Yfactor, ObservedDataSets(j).y{kk}.*errorfactor);
                                setCurveOptions(p_handle2, Curves.CurveOptions);
                            end
                        end
                    end
                    break
                end
                break
            end
        end
    end
end
if exist('p_handle')
    setCurveOptions(p_handle, Curves.CurveOptions);
else
    p_handle=[];
    legendLabels=[];
end