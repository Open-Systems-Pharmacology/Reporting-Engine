function lgd = reshapeQualificationLegend(lgd)
% Algorithm that check legend size of current figure and reshape if necessary

% Minimum FontSize allowed
MinFontSize = 6;

% Maximum number of characters allowed with minimum FontSize
LegendLengthMax = 100;

% Get figure position as [x0 y0 width height]
figPosition = get(gca, 'Position');

% If Legend caption is bigger than 10% of figure width
% Legend Fontisize is decreased until minimum allowed FontSize
while lgd.Position(3) >= 1.10*figPosition(3) && lgd.FontSize >= MinFontSize
    lgd.FontSize = lgd.FontSize-0.5;
end

% Update figure position
figPosition = get(gca, 'Position');
    
% If Legend caption is bigger than 60% of plot height
% Legend Fontisize is decreased until minimum allowed FontSize
while lgd.Position(4) >= 0.6*figPosition(4) && lgd.FontSize >= MinFontSize
    lgd.FontSize = lgd.FontSize-0.5;
end

% If Minimum of FontSize is reached and the caption is still too wide
% Split line in 2 
LegendLengths = cellfun(@(x)size(x, 2), lgd.String);
if lgd.Position(3) >= 1.10*figPosition(3)
    for LegendIndex = 1:length(lgd.String)
        if LegendLengths(LegendIndex) > LegendLengthMax
            % Split legend entry
            splitLegend = strsplit(lgd.String{LegendIndex});
            % Break line of the legend
            strjoin([splitLegend(1:ceil(length(splitLegend)/2)), ...
                newline , splitLegend(ceil(length(splitLegend)/2+1:end))]);
        end
    end
end

% If Minimum of FontSize is reached and the caption still shrinks the plot
% Set legend within the plot
% Update figure position
figPosition = get(gca, 'Position');
if lgd.Position(4) >= 0.6*figPosition(4)
    lgd.Location = 'northwest';
end