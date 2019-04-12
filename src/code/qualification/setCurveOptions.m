function setCurveOptions(pp, CurveOptions)

if isfield(CurveOptions, 'Color')
    rgb= hex2rgb(CurveOptions.Color);
    set(pp, 'Color', rgb);
end

if isfield(CurveOptions, 'LineStyle')
    if strcmp(CurveOptions.LineStyle, 'Dash')
        CurveOptions.LineStyle='--';
    elseif strcmp(CurveOptions.LineStyle, 'Dot')
        CurveOptions.LineStyle=':';
	elseif strcmp(CurveOptions.LineStyle, 'DashDot')
        CurveOptions.LineStyle='.-';
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
