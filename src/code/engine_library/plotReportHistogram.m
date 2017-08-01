function  [csvArray] = plotReportHistogram(WSettings,figureHandle,y,sourceIndex,yData,xLabel,xUnit,xCategory,legendEntries)
% PLOTREPORTHISTOGRAM creates histogramm 
%
% [csvArray] = plotReportHistogram(WSettings,figureHandle,y,yRef,xLabel,legendEntries)
% 
% Inputs: 
%       WSettings (structure)    definition of properties used in all
%                   workflow functions see GETDEFAULTWORKFLOWSETTINGS
%       y  ( double vector nInd x 1) property to plot
%      sourceIndex (double vector) (nInd x 1) index of source population
%       yData ( double vector nData x 1) property of the data population
%       xLabel (string) name of property 
%       xUnit (string) name of unit
%       legendEntries (cell array of strings) {name of population for
%                       legend, name of reference population for legend}
% Output
%   csv (cellarray) table with numeric information to the plot

% Open Systems Pharmacology Suite;  http://open-systems-pharmacology.org

% Initialize outputs
csvArray = {};

% get min and max values
minY = min([y;yData]);
maxY = max([y;yData]);

% get bins
if isempty(xCategory)
    nBin = min(20,floor(length(y)/3));
    dBin = (maxY-minY)/(nBin);
    if dBin>0
        bins = [(minY-dBin/2):dBin:(maxY+dBin/2)];
        mBins = [(minY):dBin:(maxY)];
    else
        bins = [maxY-0.5:1:maxY+0.5];
        mBins = [maxY maxY+2];
    end
else
    mBins = unique(y);
    bins = [mBins(1)-0.5; mBins+0.5];
end

% get percentage of individual per bin for population
l = zeros(length(bins),max(sourceIndex));
lnorm = zeros(length(bins),max(sourceIndex));
lgtxt = cell(1,max(sourceIndex));
for iSource = 1:max(sourceIndex)
    jj = sourceIndex== iSource;
    l(:,iSource) = histc(y(jj),bins);
    lnorm(:,iSource) = l(:,iSource)./sum(jj).*100;
    lgtxt{iSource} = legendEntries{iSource};
end
[colMap] = getcolmarkForMap(WSettings.colormapSimulation,max(sourceIndex));


% get percentage of individual per bin for dataPopulation
if ~isempty(yData)
    l(:,end+1) = histc(yData,bins);
    lnorm(:,end+1) = l(:,end)./length(yData).*100;
    lgtxt{end+1} = legendEntries{3};
    [colMap] = [colMap;getcolmarkForMap(WSettings.colormapData,1)];

end


% do the Plot
ax = getReportFigure(WSettings,1,1,figureHandle);

if isempty(xCategory)
    bar(mBins,lnorm(1:length(mBins),:));
    set(ax,'xlim',bins([1 end]));
else
    bar(1:length(mBins),lnorm(1:length(mBins),:));
    set(ax,'xlim',[0.5 length(mBins)+0.5]);
end

colormap(ax,colMap);

% set labels 
if isempty(xCategory)
    xlabel(getLabelWithUnit(xLabel,xUnit));
    legend(lgtxt,'location','northoutside');
else
    xt = mBins(1:length(bins)-1);
    set(ax,'xtick',1:length(mBins),'xticklabel',xCategory(xt));
end
ylabel('Number of individuals [%]');

legend(ax,lgtxt,'location','northoutside');

% construct csvArray
csvArray(1,:) = {getLabelWithUnit(xLabel,xUnit),legendEntries{1}};
for iBin = 1:length(bins)-1
    if isempty(xCategory)
        csvArray(iBin+1,:) = {sprintf('%.3g <= %s < %.3g',max(minY,bins(iBin)),xLabel,min(maxY,bins(iBin+1))),sprintf('%d',l(iBin,1))}; %#ok<AGROW>
    else
        csvArray(iBin+1,:) = {sprintf('%s',xCategory{mBins(iBin)}),sprintf('%d',l(iBin,1))}; %#ok<AGROW>
    end
end

if ~isempty(yData)
    csvArray{1,end+1} = legendEntries{3};
    for iBin = 1:length(bins)-1
        csvArray{iBin+1,end} = sprintf('%d',l(iBin,end)); %#ok<AGROW>
    end
end

return
