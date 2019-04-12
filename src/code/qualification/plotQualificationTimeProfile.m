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

xAxisUnit=setFigureOptions(AxesOptions);

% Perform the plot based on Curves indications
legendLabels={};
hold on;

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
                    else
                        yyaxis left
                    end
                end
                % Convert units to reference unit
                SimTL.time = ConvertTimeUnit(SimTL.time, SimTL.timeUnit, xAxisUnit);
                pp = plot(SimTL.time, SimTL.y{j});
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
                    else
                        yyaxis left
                    end
                end
                DataTP(j).time = ConvertTimeUnit(DataTP(j).time, DataTP(j).timeUnit, xAxisUnit);
                pp = plot(DataTP(j).time, DataTP(j).y{1});
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
