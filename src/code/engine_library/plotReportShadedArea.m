function csv = plotReportShadedArea(WSettings,figureHandle,x,y,yRef,xLabel,xUnit,yLabel,yUnit,yscale,legendEntries,yRefPopPK)
%PLOTREPORTSHADEDAREA Plots one property of a population vs another as shaed arae in comparison to a reference population
%
%   plotLayout_TimeProfile(time,Y,time_Ref,Y_Ref,...
%    time_label,time_unit,time_xlimit,...
%    Y_label,Y_unit,ylimit,yscale,...
%    shadedArea_percentiles,shadedArea_lineType,...
%    figure_handle,colorVector_patch,colorVector_line)
%
% Inputs:
%   WSettings (structure)    definition of properties used in all
%                   workflow functions see GETDEFAULTWORKFLOWSETTINGS
%   figure_handle handle of figure, 
%       empty a new figure is created, 
%       valid handle of a figure, the figure is cleared
%   x (double vector (nIndividual x 1)): X values of population
%   y  (double vector (nIndividual x 1)): Y values of population
%   yRef  (double matrix (nIndividual x 1)): Y values of population
%           if empty no reference population is drawn
%   yRefPopPK  (double vector (3 x 1)): Y values of popPKasComparison
%   xLabel (string)  text beside the X-axis
%   xUnit (string)   text beside the X-axis
%   yLabel (string)  text beside the Y-axis
%   yUnit (string)   text beside the Y-axis
%   yscale (string): may be 'lin' or 'log'
%   legendEntries (cell array of strings) {name of population for  legend,
%                       name of reference population for legend}

% Open Systems Pharmacology Suite;  http://open-systems-pharmacology.org

if ~exist('yRefPopPK','var')
    yRefPopPK = nan;
end


% Binning
% nBins  bins with equally distributed entries
nBins = 20;
% at least three entries per bin
nBins=floor(min(nBins,length(x)/3));

% calulate binBorders
binBorders=prctile(x,100/nBins.*[0:nBins]); %#ok<NBRAK>
binBorders(end) = ceil(binBorders(end));
binBorders(1) = floor(binBorders(1));

% Initialize arrays
binMedian = nan(nBins,1);
yMean = nan(nBins,1);
yMin = nan(nBins,1);
yMax = nan(nBins,1);
nEntry = nan(nBins,1);
    
for iBin=1:length(binBorders)-1
    jj = x >= binBorders(iBin) &  x < binBorders(iBin+1);
    
    [yMin(iBin),yMean(iBin),yMax(iBin),legendTextMean,rangeTxt,csvHeader] = getRangePlotPercentiles(WSettings,y(jj));
    
    binMedian(iBin) = median(x(jj));
    nEntry(iBin) = sum(jj); 
end


% Reference population
if ~isempty(yRef) 
    [yMinRef,yMeanRef,yMaxRef] = getRangePlotPercentiles(WSettings,yRef);
else
    yMinRef = nan;
    yMaxRef = nan;
    yMeanRef = nan;
end

% get Csv Array beforr y values are tarnsformed for y scale
csv = getCsvArray(WSettings,yMin,yMean,yMax,yLabel,yUnit,xLabel,xUnit,binBorders,nEntry,csvHeader,binMedian,yMinRef,yMaxRef,yMeanRef);


% for yscal log transfer values to logscale
if strcmp(yscale,'log')
    yMin = log10(yMin);
    yMax = log10(yMax);
    yMean = log10(yMean);
            
    if ~isnan(yMeanRef)
        yMinRef = log10(yMinRef);
        yMaxRef = log10(yMaxRef);
        yMeanRef = log10(yMeanRef);
        

    end
    
    if ~isnan(yRefPopPK)
        yRefPopPK = log10(yRefPopPK);
    end
    
end



% create figure
ax = getReportFigure(WSettings,1,1,figureHandle,'figureformat','landscape');
        
% exclude log(0)
ij = find(~isinf(yMin) & ~isinf(yMax));


% Patch of reference population
% Plot before main population, so it is in background
if ~isnan(yMinRef)
    
    lgh(4) = patch(binMedian([1 end end 1]),[yMinRef yMinRef yMaxRef yMaxRef],WSettings.colorVector(2,:),...
        'linestyle','-','LineWidth',0.5,'FaceAlpha',WSettings.FaceAlpha,...
        'EdgeColor',WSettings.colorVector(2,:),'displayname',strtrim(strrep(legendEntries{2},'<xxx>',rangeTxt)));

    % Mean of reference population
    lgh(3)=plot([binMedian(1) binMedian(end)],yMeanRef*[1 1],'-','color',WSettings.colorVector(2,:),'linewidth',2,...
        'displayname',strtrim(strrep(legendEntries{2},'<xxx>',legendTextMean)));

end


% plot PopPK Reference
if any(~isnan(yRefPopPK))
    lgh(5)=plot([binMedian(1) binMedian(end)],yRefPopPK(2)*[1 1],'-','color',WSettings.colorVector(3,:),'linewidth',2,...
        'displayname',strtrim(strrep(legendEntries{3},'<xxx>',legendTextMean)));
    plot([binMedian(1) binMedian(end)],yRefPopPK(1)*[1 1],'-','color',WSettings.colorVector(3,:),'linewidth',1);
    plot([binMedian(1) binMedian(end)],yRefPopPK(3)*[1 1],'-','color',WSettings.colorVector(3,:),'linewidth',1);
end

lgh(2) = patch([binMedian(ij); binMedian(ij(end:-1:1))],[yMin(ij); yMax(ij(end:-1:1))],WSettings.colorVector(1,:),...
    'linestyle','-','LineWidth',0.5,'FaceAlpha',WSettings.FaceAlpha,...
    'EdgeColor',WSettings.colorVector(1,:),'displayname',strtrim(strrep(legendEntries{1},'<xxx>',rangeTxt)));

% Plot mean
lgh(1)=plot(binMedian,yMean,'-','color',WSettings.colorVector(1,:),'linewidth',2,...
    'displayname',strtrim(strrep(legendEntries{1},'<xxx>',legendTextMean)));



% Set axes labels
xlabel(getLabelWithUnit(xLabel,xUnit));
ylabel(getLabelWithUnit(yLabel,yUnit));

% set(ax,'xlim',[floor(binMedian(1)) ceil(binMedian(end))]);
set(ax,'xlim',[binMedian(1) binMedian(end)]);

if strcmp(yscale,'log')
    setLogarithmicYticks(ax);    
end

% set legend
jj = isgraphics(lgh);
legend(lgh(jj),get(lgh(jj),'displayname'),'location','northoutside');
    

return

function csvArray = getCsvArray(WSettings,yMin,yMean,yMax,yLabel,yUnit,xLabel,xUnit,binBorders,nEntry,csvHeader,binMedian,yMinRef,yMaxRef,yMeanRef)

csvArray={};

csvArray{1,2} = getLabelWithUnit(yLabel,yUnit);
csvArray{1,1} = getLabelWithUnit(xLabel,xUnit);

csvArray(2,1:3+length(csvHeader)) = [{sprintf('Range of %s',xLabel)},{sprintf('Median of %s',xLabel)},csvHeader,{'N individuals'}];

binLower = binBorders(1:end-1);
binUpper = binBorders(2:end);

for iBin = 1:length(binLower)
    csvArray(2+iBin,:) = ...
        {sprintf([WSettings.csvFormat '-' WSettings.csvFormat],binLower(iBin),binUpper(iBin)),...
        sprintf(WSettings.csvFormat,binMedian(iBin)),sprintf(WSettings.csvFormat,yMin(iBin)),...
        sprintf(WSettings.csvFormat,yMean(iBin)),sprintf(WSettings.csvFormat,yMax(iBin)),sprintf('%d',nEntry(iBin))}; 
end

if ~isnan(yMinRef)
    iRow=size(csvArray,1)+2;
    csvArray{iRow,1}='reference values:';
    csvArray(iRow,2:4) = {sprintf(WSettings.csvFormat,yMinRef),sprintf(WSettings.csvFormat,yMeanRef),...
        sprintf(WSettings.csvFormat,yMaxRef)}; 
end

return
