function setLogarithmicYticks(ax)
% set yticks to an axes, where the conten where plotted as log10
% 
%  setLogarithmicYticks(ax)
% 
% Inputs ax = axes handle

% Open Systems Pharmacology Suite;  http://open-systems-pharmacology.org


yt = get(ax,'ytick');
yl = get(ax,'ylim');

% get new integer ticks
ytNew = intersect(yt,ceil(yt(1)):floor(yt(end)));

% if none or only one add outer limits
if  length(ytNew) <=1
    ytNew = unique([floor(yt(1)),ytNew,ceil(yt(end))]);
    yl = [ytNew(1) ytNew(end)];
    set(ax,'ylim',yl)
end

% switch to string if numbers to large
if any(abs(ytNew)>=5)
    ytl = cellfun(@(x) sprintf('%.0e',x),num2cell(10.^ytNew),'uniformoutput',false);
    set(ax,'ytick',ytNew,'yticklabel',ytl);
else
    set(ax,'ytick',ytNew,'yticklabel',10.^ytNew);
end

% get Minor ticks of yt
yMinorTicks = 10.^(ytNew(1)).*[0.1:0.1:1];
for iT = 1:length(ytNew)
    yMinorTicks = [yMinorTicks 10.^(ytNew(iT)).*[2:1:10]]; %#ok<AGROW>
end
jj = log10(yMinorTicks) >= yl(1) &  log10(yMinorTicks) <= yl(end);
yMinorTicks = log10(yMinorTicks(jj));
yMinorTicks = setdiff(yMinorTicks,ytNew);

xl = get(ax,'xlim');
for iT = 1:length(yMinorTicks)
    plot(xl(1)+[0 0.008*range(xl)],yMinorTicks(iT)*[1 1],'k-');
    plot(xl(2)-[0 0.008*range(xl)],yMinorTicks(iT)*[1 1],'k-');
end

