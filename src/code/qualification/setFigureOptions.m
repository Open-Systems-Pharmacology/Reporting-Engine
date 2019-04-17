function [xUnit, yUnit, y2Unit] = setFigureOptions(AxesOptions)

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

xUnit=[];
yUnit=[];
y2Unit=[];

% Set all axis options based on TimeProfile.Plot.Axes

% Set x boundaries
if isfield(xAxesOptions, {'Min', 'Max'})
    xlim([xAxesOptions.Min xAxesOptions.Max]);
elseif isfield(xAxesOptions, 'Min')
    xlim([xAxesOptions.Min Inf]);
elseif isfield(xAxesOptions, 'Max')
    xlim([-Inf xAxesOptions.Max]);
end

% Check if 2 y-axes
if exist('yyAxesOptions')
    yyaxis right
    % Set yy boundaries
    if isfield(yyAxesOptions, {'Min', 'Max'})
        ylim([yyAxesOptions.Min yyAxesOptions.Max]);
    elseif isfield(yyAxesOptions, 'Min')
        ylim([yyAxesOptions.Min Inf]);
    elseif isfield(yyAxesOptions, 'Max')
        ylim([-Inf yyAxesOptions.Max]);
    end
    yyaxis left
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
    xUnit = xAxesOptions.Unit;
elseif isfield(xAxesOptions, 'Dimension')
    xlabel(xAxesOptions.Dimension);
end

% Set y label
if exist('yyAxesOptions')
    yyaxis right
    if isfield(yyAxesOptions, {'Dimension', 'Unit'})
        ylabel([yyAxesOptions.Dimension ' [' yyAxesOptions.Unit ']']);
    elseif isfield(yyAxesOptions, 'Dimension')
        ylabel(yyAxesOptions.Dimension);
    end
    yyaxis left
end
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
if exist('yyAxesOptions')
    yyaxis right
    if isfield(yyAxesOptions, 'GridLines')
        if yyAxesOptions.GridLines==1
            set(gca, 'YGrid', 'on');
        end
    end
    yyaxis left
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

if exist('yyAxesOptions')
    yyaxis right
    if isfield(yyAxesOptions, 'Scaling')
        if strcmp(yyAxesOptions.Scaling, 'Log') || strcmp(yyAxesOptions.Scaling, 'log')
            set(gca, 'YScale', 'log')
        end
    end
    yyaxis left
end
if isfield(yAxesOptions, 'Scaling')
    if strcmp(yAxesOptions.Scaling, 'Log') || strcmp(yAxesOptions.Scaling, 'log')
        set(gca, 'YScale', 'log')
    end
end
