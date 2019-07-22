function lgd = reshapeQualificationLegend(lgd, figureHandle)
% Algorithm that break lines if is too long

% Get figure position as [x0 y0 width height]
figPosition = get(figureHandle, 'Position');

% If caption is too wide split line in half 
LegendLengths = cellfun(@(x)size(x, 2), lgd.String);
% Get the number of line breaks:
WidthRatio = lgd.Position(3)/figPosition(3);

% Get the legend length ratio compared to Figure Width
% Add an offset of 10% to ensure that all the legend is visible
LegendLengthsRatios = ceil(WidthRatio.*LegendLengths./max(LegendLengths)+0.10);

for LegendIndex = 1:length(lgd.String)
    switch LegendLengthsRatios(LegendIndex)
        case 2
            % Since character '-' was often use to split different
            % information, it was used as splitting character
            LegendSplit = strsplit(lgd.String{LegendIndex}, '-');
            if length(LegendSplit)<2
                lgd.String{LegendIndex} = [lgd.String{LegendIndex}(1:ceil(LegendLengths(LegendIndex)/2)) ...
                newline LegendSplit(ceil(LegendLengths(LegendIndex)/2)+1:end)];
            else
                lgd.String{LegendIndex} = strjoin([LegendSplit(1:ceil(length(LegendSplit)/2)) ...
                strjoin([newline LegendSplit(ceil(length(LegendSplit)/2)+1:end)],'')], '-');
            end
        case 3
            LegendSplit = strsplit(lgd.String{LegendIndex}, '-');
            if length(LegendSplit)<3
                lgd.String{LegendIndex} = [LegendSplit(1:ceil(LegendLengths(LegendIndex)/3)) ...
                newline LegendSplit(ceil(LegendLengths(LegendIndex)/3)+1:ceil(2*LegendLengths(LegendIndex)/3)) ...
                newline LegendSplit(ceil(2*LegendLengths(LegendIndex)/3)+1:end)];
            else
                lgd.String{LegendIndex} = strjoin([LegendSplit(1:ceil(length(LegendSplit)/3)) ...
                strjoin([newline LegendSplit(ceil(length(LegendSplit)/3)+1:ceil(2*length(LegendSplit)/3))],'') ...
                strjoin([newline LegendSplit(ceil(2*length(LegendSplit)/3)+1:end)]),''], '-');
            end
        case 4
            LegendSplit = strsplit(lgd.String{LegendIndex}, '-');
            if length(LegendSplit)<4
                lgd.String{LegendIndex} = [LegendSplit(1:ceil(LegendLengths(LegendIndex)/4)) ...
                newline LegendSplit(ceil(LegendLengths(LegendIndex)/4)+1:ceil(2*LegendLengths(LegendIndex)/4)) ...
                newline LegendSplit(ceil(2*LegendLengths(LegendIndex)/4)+1:ceil(3*LegendLengths(LegendIndex)/4)) ...
                newline LegendSplit(ceil(3*LegendLengths(LegendIndex)/4)+1:end)];
            else
                lgd.String{LegendIndex} = strjoin([LegendSplit(1:ceil(length(LegendSplit)/4)) ...
                strjoin([newline LegendSplit(ceil(length(LegendSplit)/4)+1:ceil(2*length(LegendSplit)/4))],'') ...
                strjoin([newline LegendSplit(ceil(2*length(LegendSplit)/4)+1:ceil(3*length(LegendSplit)/4))],'') ...
                strjoin([newline LegendSplit(ceil(3*length(LegendSplit)/4)+1:end)],'')], '-');
            end
    end
end

lgd.Location = 'South';
set(figureHandle, 'Position', get(lgd, 'Position'));

CurvesProperties = get(figureHandle,'children'); % Get the curves properties
AnnotationTags = findall(gcf,'Type', 'Textbox'); % Get the Watermarks and other annotations

for AnnotationIndex=1:length(AnnotationTags)
    set(AnnotationTags(AnnotationIndex), 'visible', 'off');
end

for CurveIndex=1:length(CurvesProperties)
    CurvesProperties(CurveIndex).XData = NaN*CurvesProperties(CurveIndex).XData;
end
axis off;