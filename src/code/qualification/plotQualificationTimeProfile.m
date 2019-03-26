function plotQualificationTimeProfile(WSettings,figureHandle,SimTL,DataTP, Curves, AxesOptions)
%PLOTQUALIFICATIONTIMEPROFILE Plots the time profile of a population in comparison to a reference population
%
% csv =  =  plotQualificationTimeProfile(WSettings,figureHandle,SimTL,DataTP, Curves, AxesOptions)
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
ax = getReportFigure(WSettings,1,1,figureHandle,'figureformat','landscape');

setFigureOptions(AxesOptions);

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
                pp = plot(SimTL.time, SimTL.y{j});
            end
        end
        % For observations
    else
        for j=1:length(DataTP)
            % Get the right observation profile
            findPathOutput = strfind(Curves(i).Y,DataTP(j).Id);
            if ~isempty(findPathOutput)
                legendLabels{length(legendLabels)+1}=Curves(i).Name;
                pp = plot(DataTP(j).time, DataTP(j).y{1});
            end
        end
    end
    setCurveOptions(pp, Curves(i).CurveOptions);
end
legend(legendLabels);

%------------------------------------------------------------------------
% Auxiliary functions

function setFigureOptions(AxesOptions)

% Axes definitions and options
for iAxes=1:length(AxesOptions)
    
    if iscell(AxesOptions(iAxes))
        AxesOptionsStr=AxesOptions{iAxes};
    else
        AxesOptionsStr=AxesOptions(iAxes);
    end
    switch AxesOptionsStr.Type
        case 'X'
            % x axis options
            xAxesOptions=AxesOptionsStr;
        case 'Y'
            % y axis options
            yAxesOptions=AxesOptionsStr;
        case 'Y2'
            % y2 axis options (right axis)
            yyAxesOptions=AxesOptionsStr;
    end
end

% Set all axis options based on TimeProfile.Plot.Axes

% Set x boundaries
if isfield(xAxesOptions, {'Min', 'Max'})
    xlim([xAxesOptions.Min xAxesOptions.Max]);
elseif isfield(xAxesOptions, 'Min')
    xlim([xAxesOptions.Min Inf]);
elseif isfield(xAxesOptions, 'Max')
    xlim([-Inf xAxesOptions.Max]);
end

% Set y boundaries
if isfield(yAxesOptions, {'Min', 'Max'})
    ylim([yAxesOptions.Min yAxesOptions.Max]);
elseif isfield(yAxesOptions, 'Min')
    ylim([yAxesOptions.Min Inf]);
elseif isfield(yAxesOptions, 'Max')
    ylim([-Inf yAxesOptions.Max]);
end

% Set x label
if isfield(xAxesOptions, {'Dimension', 'Unit'})
    xlabel([xAxesOptions.Dimension ' [' xAxesOptions.Unit ']']);
elseif isfield(xAxesOptions, 'Dimension')
    xlabel(xAxesOptions.Dimension);
end

% Set y label
if isfield(yAxesOptions, {'Dimension', 'Unit'})
    ylabel([yAxesOptions.Dimension ' [' yAxesOptions.Unit ']']);
elseif isfield(yAxesOptions, 'Dimension')
    ylabel(yAxesOptions.Dimension);
end

% Set Grid options
if isfield(xAxesOptions, 'GridLines')
    if xAxesOptions.GridLines==1
        set(gca, 'XGrid', 'on');
    end
end
if isfield(yAxesOptions, 'GridLines')
    if yAxesOptions.GridLines==1
        set(gca, 'YGrid', 'on');
    end
end

% Set scaling
if isfield(xAxesOptions, 'Scaling')
    if strcmp(xAxesOptions.Scaling, 'Log') || strcmp(xAxesOptions.Scaling, 'log')
        set(gca, 'XScale', 'log')
    end
end

if isfield(yAxesOptions, 'Scaling')
    if strcmp(yAxesOptions.Scaling, 'Log') || strcmp(yAxesOptions.Scaling, 'log')
        set(gca, 'YScale', 'log')
    end
end

function setCurveOptions(pp, CurveOptions)

if isfield(CurveOptions, 'Color')
    rgb= hex2rgb(CurveOptions.Color);
    set(pp, 'Color', rgb);
end

if isfield(CurveOptions, 'LineStyle')
    if strcmp(CurveOptions.LineStyle, 'Dash')
        CurveOptions.LineStyle='--';
    end
    set(pp, 'LineStyle', CurveOptions.LineStyle);
end

if isfield(CurveOptions, 'Symbol')
    % Some symbols are not accounted by matlab
    if strcmp(CurveOptions.Symbol, 'Circle')
        CurveOptions.Symbol='o';
    elseif strcmp(CurveOptions.Symbol, 'Square')
        CurveOptions.Symbol='s';
    elseif strcmp(CurveOptions.Symbol, 'Dot')
        CurveOptions.Symbol='.';
    end
    
    set(pp, 'Marker', CurveOptions.Symbol);
end
