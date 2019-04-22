function plotQualificationTimeProfile(WSettings,figureHandle,SimTL,DataTP, Curves, AxesOptions, PlotSettings)
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


% Perform the plot based on Curves indications
legendLabels={};
hold on;

% Get dimensions for scaling plots
dimensionList=getDimensions;

% Plot the Curves indicated by the array of structures Curves
for i=1:length(Curves)
    % Check if data is from simulation or observation
    pathOrigin = strfind(Curves(i).Y,SimTL.Id);
    
    % For simulations
    if ~isempty(pathOrigin)
        for j=1:length(SimTL.outputPathList)
            % Get the right simulation curve
            findPathOutput = strfind(Curves(i).Y,SimTL.outputPathList{j});
            if ~isempty(findPathOutput)
                legendLabels{length(legendLabels)+1}=Curves(i).Name;
                
                if isfield(Curves(i).CurveOptions, 'yAxisType')
                    if strcmp(Curves(i).CurveOptions.yAxisType, 'Y2')
                        yyaxis right
                        YDimension=dimensionList{strContains(yyAxesOptions.Dimension, dimensionList)};
                        Yfactor=getUnitFactor(SimTL.outputUnit{j},yyAxesOptions.Unit,YDimension);
                    else
                        yyaxis left
                    end
                else
                    YDimension=dimensionList{strContains(yAxesOptions.Dimension, dimensionList)};
                    Yfactor=getUnitFactor(SimTL.outputUnit{j},yAxesOptions.Unit,YDimension);
                end
                
                % Convert units to reference unit
                XDimension=dimensionList{strContains(xAxesOptions.Dimension, dimensionList)};
                Xfactor=getUnitFactor(SimTL.timeUnit,xAxesOptions.Unit,XDimension);
                
                pp = plot(SimTL.time.*Xfactor, SimTL.y{j}.*Yfactor);
                break
            end
        end
        % For observations
    else
        for j=1:length(DataTP)
            % Get the right observation profile
            findPathOutput = strfind(Curves(i).Y,DataTP(j).Id);
            
            if ~isempty(findPathOutput)
                legendLabels{length(legendLabels)+1}=Curves(i).Name;
                
                if isfield(Curves(i).CurveOptions, 'yAxisType')
                    if strcmp(Curves(i).CurveOptions.yAxisType, 'Y2')
                        yyaxis right
                        YDimension=dimensionList{strContains(yyAxesOptions.Dimension, dimensionList)};
                        Yfactor=getUnitFactor(DataTP(j).outputUnit,yyAxesOptions.Unit,YDimension);
                    else
                        yyaxis left
                    end
                else
                    YDimension=dimensionList{strContains(yAxesOptions.Dimension, dimensionList)};
                    % Caution: Code to be updated
                    % Observed concentration was massic whereas output was
                    % molar
                    Yfactor=getUnitFactor(DataTP(j).outputUnit{1},yAxesOptions.Unit,YDimension, 'MW',325.78);
                    
                end
                XDimension=dimensionList{strContains(xAxesOptions.Dimension, dimensionList)};
                Xfactor=getUnitFactor(DataTP(j).timeUnit,xAxesOptions.Unit,XDimension);
                pp = plot(DataTP(j).time.*Xfactor, DataTP(j).y{1}.*Yfactor);
                break
            end
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
legend(legendLabels);
