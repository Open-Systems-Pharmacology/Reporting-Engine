function setCurveOptions(pp, CurveOptions)
% setCurveOptions Auxiliary to Set propeerties of plotted curves

% If CurveOptions is directly a plot handle, 
% Get the plot handle properties as a structure

if ~isstruct(CurveOptions)
    CurveOptions = struct(CurveOptions);
end

% Set plot color
if isfield(CurveOptions, 'Color')
    % If color is in hex translate to RGB, else use RGB
    if ischar(CurveOptions.Color)
        rgb = hex2rgb(CurveOptions.Color);
    else
        rgb = CurveOptions.Color;
    end
    set(pp, 'Color', rgb);
end

if isfield(CurveOptions, 'LineStyle')
    if strcmp(CurveOptions.LineStyle, 'Dash')
        CurveOptions.LineStyle='--';
    elseif strcmp(CurveOptions.LineStyle, 'Dot')
        CurveOptions.LineStyle=':';
    elseif strcmp(CurveOptions.LineStyle, 'DashDot')
        CurveOptions.LineStyle='-.';
    elseif strcmp(CurveOptions.LineStyle, 'Solid')
        CurveOptions.LineStyle='-';
    end
    set(pp, 'LineStyle', CurveOptions.LineStyle);
end

if isfield(CurveOptions, 'Symbol')
    % Some symbols are not accounted by matlab
    if strcmp(CurveOptions.Symbol, 'Circle')
        CurveOptions.Symbol='o';
    elseif strcmp(CurveOptions.Symbol, 'Asterisk')
        CurveOptions.Symbol='*';
    elseif strcmp(CurveOptions.Symbol, 'Point')
        CurveOptions.Symbol='.';
    elseif strcmp(CurveOptions.Symbol, 'Dot')
        CurveOptions.Symbol='.';
    elseif strcmp(CurveOptions.Symbol, 'Cross')
        CurveOptions.Symbol='x';
    elseif strcmp(CurveOptions.Symbol, 'Triangle')
        CurveOptions.Symbol='^';
    end
    
    set(pp, 'Marker', CurveOptions.Symbol);
end

if isfield(CurveOptions, 'Marker')
    set(pp, 'Marker', CurveOptions.Marker);
end
