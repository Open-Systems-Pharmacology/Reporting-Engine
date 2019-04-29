function [xAxesOptions, yAxesOptions, yyAxesOptions] = setFigureOptions(AxesOptions)

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

% Associate empty yyAxesOptions if does not exist
if ~exist('yyAxesOptions')
    yyAxesOptions=[];
end

% Set x boundaries
if isfield(xAxesOptions, {'Min', 'Max'})
    xlim([xAxesOptions.Min xAxesOptions.Max]);
elseif isfield(xAxesOptions, 'Min')
    xlim([xAxesOptions.Min Inf]);
elseif isfield(xAxesOptions, 'Max')
    xlim([-Inf xAxesOptions.Max]);
end

% Check if 2 y-axes
if ~isempty(yyAxesOptions)
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
if ~isfield(xAxesOptions, 'Unit')
    xAxesOptions.Unit=[];
end
if isfield(xAxesOptions, 'Dimension')
    xLabelFinal = getLabelWithUnit(xAxesOptions.Dimension,xAxesOptions.Unit);
    xlabel(xLabelFinal);
end

% Set y label
if ~isempty(yyAxesOptions)
    yyaxis right
    if ~isfield(yyAxesOptions, 'Unit')
        yyAxesOptions.Unit=[];
    end
    if isfield(yyAxesOptions, 'Dimension')
        yyLabelFinal = getLabelWithUnit(yyAxesOptions.Dimension,yyAxesOptions.Unit);
        ylabel(yyLabelFinal);
    end
    yyaxis left
end
if ~isfield(yAxesOptions, 'Unit')
    yAxesOptions.Unit=[];
end
if isfield(yAxesOptions, 'Dimension')
    yLabelFinal = getLabelWithUnit(yAxesOptions.Dimension,yAxesOptions.Unit);
    ylabel(yLabelFinal);
end


% Set Grid options
if isfield(xAxesOptions, 'GridLines')
    if xAxesOptions.GridLines==1
        set(gca, 'XGrid', 'on');
    end
end
if ~isempty(yyAxesOptions)
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
        set(gca, 'XScale', 'log');
    end
end

if ~isempty(yyAxesOptions)
    yyaxis right
    if isfield(yyAxesOptions, 'Scaling')
        if strcmp(yyAxesOptions.Scaling, 'Log') || strcmp(yyAxesOptions.Scaling, 'log')
            set(gca, 'YScale', 'log');
        end
    end
    yyaxis left
end
if isfield(yAxesOptions, 'Scaling')
    if strcmp(yAxesOptions.Scaling, 'Log') || strcmp(yAxesOptions.Scaling, 'log')
        set(gca, 'YScale', 'log');
    end
end
