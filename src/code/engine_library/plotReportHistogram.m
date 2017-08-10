function  csvArray = plotReportHistogram(WSettings,figureHandle,y,sourceIndex,yData,xLabel,xUnit,xCategory,legendEntries,flag)
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
%       flag (string) default 'none';
%                  'isResiduals' : histogram is shifted 0 is in the middle,
%                               a normal distribution is plotted
%
% Output
%   csv (cellarray) table with numeric information to the plot

% Open Systems Pharmacology Suite;  http://open-systems-pharmacology.org

% Initialize outputs
csvArray = {};

if ~exist('flag','var')
    flag = 'none';
end
switch flag
    case 'isResiduals'
        barType = 'stacked';
        ylabeltxt = 'Number of residuals [%]';
    otherwise
        barType = 'grouped';
        ylabeltxt = 'Number of individuals [%]';
end

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
    lgh = bar(mBins,lnorm(1:length(mBins),:),barType);
    set(ax,'xlim',bins([1 end]));
else
    lgh = bar(1:length(mBins),lnorm(1:length(mBins),:),barType);
    set(ax,'xlim',[0.5 length(mBins)+0.5]);
end
for iL = 1:length(lgh);
    set(lgh(iL),'displayname',lgtxt{iL});
end

colormap(ax,colMap);

% set labels 
if isempty(xCategory)
    xlabel(getLabelWithUnit(xLabel,xUnit));
else
    xt = mBins(1:length(bins)-1);
    set(ax,'xtick',1:length(mBins),'xticklabel',xCategory(xt));
end
ylabel(ylabeltxt);

% flag specific add ons
if strcmp(flag,'isResiduals')
        xl = get(ax,'xlim');
        xl = max(abs(xl)).*[-1 1];
        set(ax,'xlim',xl);
        yl = get(ax,'ylim');
        
        xNormal = xl(1):range(xl)/100:xl(2);
        yNormal = pdf('normal',xNormal,0,1);
        
        lgh(end+1) = plot(xNormal,yNormal.*max(sum(lnorm,2))/max(yNormal),'-k','linewidth',2,'displayname','normal distribution');
        plot([0 0],yl,'-k','linewidth',2)
        
end
legend(lgh,get(lgh,'displayname'),'location','northoutside');





%% construct csvArray
if strcmp(flag,'isResiduals')
%         no csv
        return
end

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
