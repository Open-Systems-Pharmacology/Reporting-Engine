function  [csvArray] = plotReportHistogram(Settings,figure_handle,y,yRef,xLabel,xUnit,popReportName,refPopReportName)
% PLOTREPORTHISTOGRAM creates histogramm 
%
% Inputs: 
%       Settings (structure)    definition of properties used in all
%                   workflow functions see GETDEFAULTWORKFLOWSETTINGS
%       y  ( double vector) property to plot
%       yRef ( double vector) property of the reference population
%       xLabel (string) name of property 
%       unit (name of unit)

% Open Systems Pharmacology Suite;  http://forum.open-systems-pharmacology.org
% Date: 26-July-2017

% Initialize outputs
csvArray = {};

% get min and max values
minY = min([y,yRef]);
maxY = max([y,yRef]);

% get bins
nBin = min(20,floor(length(y)/3));
dBin = (maxY-minY)/(nBin);
bins = [(minY-dBin/2):dBin:(maxY+dBin/2)];
mBins = [(minY):dBin:(maxY)];

% get percentage of individual per bin for population
l(:,1) = histc(y,bins);
nIndPerPop(1) = length(y);
lnorm(:,1) = l(:,1)./nIndPerPop(1).*100;
lgtxt = {sprintf('%s: n=%d',popReportName,nIndPerPop(1))};

% get percentage of individual per bin for referencePopulation
if ~isempty(yRef)
    l(:,2) = histc(yRef,bins);
    nIndPerPop(2) = length(yRef);
    lnorm(:,2) = l(:,2)./nIndPerPop(2).*100;
    lgtxt{2} = sprintf('%s: n=%d',refPopReportName,nIndPerPop(2));
end

% do the Plot
ax = getReportFigure(Settings,1,1,figure_handle);

bar(mBins,lnorm(1:length(mBins),:));

colormap(ax,Settings.colorVector(1:size(l,2),:));
set(ax,'xlim',bins([1 end]));

% set labels and color
if isempty(xUnit)
    xlabel_txt = xLabel;
else
    xlabel_txt = sprintf('%s [%s]',xLabel,xUnit);
end
xlabel(xlabel_txt);
ylabel('Number of individuals [%]');
legend(lgtxt,'location','northoutside');


% construct csvArray
csvArray(1,:) = {xlabel_txt,sprintf('Number of individuals of %s',popReportName)};
for iBin = 1:length(bins)-1
    csvArray(iBin+1,:) = {sprintf('%.3g <= %s < %.3g',bins(iBin),xLabel,bins(iBin+1)),sprintf('%d',l(iBin,1))}; %#ok<AGROW>
end

if ~isempty(yRef)
    csvArray{1,end+1} = sprintf('Number of individuals of %s',refPopReportName);
    for iBin = 1:length(bins)-1
        csvArray{iBin+1,end} = sprintf('%d',l(iBin,2)); %#ok<AGROW>
    end
end
    
return
