function ImageLegend = Legend2Image(lgd)
% Algorithm that break lines if is too long

% Get figure position as [x0 y0 width height]
figHandle = gcf;
figPosition = figHandle.Position;

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
                newline lgd.String{LegendIndex}(ceil(LegendLengths(LegendIndex)/2)+1:end)];
            else
                lgd.String{LegendIndex} = strjoin([LegendSplit(1:ceil(length(LegendSplit)/2)) ...
                strjoin([newline LegendSplit(ceil(length(LegendSplit)/2)+1:end)],'')], '-');
            end
        case 3
            LegendSplit = strsplit(lgd.String{LegendIndex}, '-');
            if length(LegendSplit)<3
                lgd.String{LegendIndex} = [lgd.String{LegendIndex}(1:ceil(LegendLengths(LegendIndex)/3)) ...
                newline lgd.String{LegendIndex}(ceil(LegendLengths(LegendIndex)/3)+1:ceil(2*LegendLengths(LegendIndex)/3)) ...
                newline lgd.String{LegendIndex}(ceil(2*LegendLengths(LegendIndex)/3)+1:end)];
            else
                lgd.String{LegendIndex} = strjoin([LegendSplit(1:ceil(length(LegendSplit)/3)) ...
                strjoin([newline LegendSplit(ceil(length(LegendSplit)/3)+1:ceil(2*length(LegendSplit)/3))],'') ...
                strjoin([newline LegendSplit(ceil(2*length(LegendSplit)/3)+1:end)],'')], '-');
            end
        case 4
            LegendSplit = strsplit(lgd.String{LegendIndex}, '-');
            if length(LegendSplit)<4
                lgd.String{LegendIndex} = [lgd.String{LegendIndex}(1:ceil(LegendLengths(LegendIndex)/4)) ...
                newline lgd.String{LegendIndex}(ceil(LegendLengths(LegendIndex)/4)+1:ceil(2*LegendLengths(LegendIndex)/4)) ...
                newline lgd.String{LegendIndex}(ceil(2*LegendLengths(LegendIndex)/4)+1:ceil(3*LegendLengths(LegendIndex)/4)) ...
                newline lgd.String{LegendIndex}(ceil(3*LegendLengths(LegendIndex)/4)+1:end)];
            else
                lgd.String{LegendIndex} = strjoin([LegendSplit(1:ceil(length(LegendSplit)/4)) ...
                strjoin([newline LegendSplit(ceil(length(LegendSplit)/4)+1:ceil(2*length(LegendSplit)/4))],'') ...
                strjoin([newline LegendSplit(ceil(2*length(LegendSplit)/4)+1:ceil(3*length(LegendSplit)/4))],'') ...
                strjoin([newline LegendSplit(ceil(3*length(LegendSplit)/4)+1:end)],'')], '-');
            end
    end
end

% getframe uses units in pixels
set(lgd, 'Units', 'pixels');

% Get a snapshot of legend only
% Caution, getframe works only if legend is within plot limits
FrameLegendPosition = lgd.Position;
FrameLegend = getframe(gcf, FrameLegendPosition);
ImageLegend = FrameLegend.cdata;

