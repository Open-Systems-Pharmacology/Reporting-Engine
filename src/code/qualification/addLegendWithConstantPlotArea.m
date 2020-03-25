function addLegendWithConstantPlotArea(currentFig,currentAxis,leg_labels)
% ADDS LEGEND BY KEEPING THE PLOTTING AREA CONSTANT

%set(currentFig,'visible','on');

legend(currentAxis,'off');
set(currentFig, 'unit', 'inches');
% get the original size of figure before the legends are added
figure_size =  get(currentFig, 'position');
% add legends and get its handle
h_legend = legend(leg_labels);
set(h_legend, 'location', 'northeastoutside');
% set unit for legend size to inches
set(h_legend, 'unit', 'inches');
% get legend size and font size
legend_size = get(h_legend, 'position');
legend_FontSize = get(h_legend,'FontSize');

% decrease font size of legend if vertical length of legend is too large
while 1.1*legend_size(4) > figure_size(4)
    set(h_legend,'FontSize',legend_FontSize-0.1)
    legend_size = get(h_legend,'position');
    legend_FontSize = get(h_legend,'FontSize');
end

% new figure width
figure_size(3) = figure_size(3) + legend_size(3);

% set new figure size
set(currentFig, 'position', figure_size)

end