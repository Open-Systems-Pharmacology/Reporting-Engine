function csv = plotReportShadedArea(WSettings,figureHandle,x,y,yRef,xLabel,xUnit,yLabel,yUnit,yscale,legendEntries)
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
%   xLabel (string)  text beside the X-axis
%   xUnit (string)   text beside the X-axis
%   yLabel (string)  text beside the Y-axis
%   yUnit (string)   text beside the Y-axis
%   yscale (string): may be 'lin' or 'log'
%   legendEntries (cell array of strings) {name of population for  legend,
%                       name of reference population for legend}

% Open Systems Pharmacology Suite;  http://open-systems-pharmacology.org




% Binning
% nBins  bins with equally distributed entries
nBins = 20;
% at least three entries per bin
nBins=floor(min(nBins,length(x)/3));

% calulate binBorders
binBorders=prctile(x,100/nBins.*[0:nBins]); %#ok<NBRAK>
binBorders(end) = binBorders(end)+0.01*(binBorders(end)-binBorders(end-1));
binBorders(1) = binBorders(1)-0.01*(binBorders(2)-binBorders(1));

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

% add edges
binMedian=[binBorders(1); binMedian; binBorders(end)];
yMean = [yMean(1); yMean; yMean(end)];
yMin = [yMin(1); yMin; yMin(end)];
yMax = [yMax(1); yMax; yMax(end)];


% Reference population
if ~isempty(yRef) 
    [yMinRef,yMeanRef,yMaxRef] = getRangePlotPercentiles(WSettings,yRef);

end

% get Csv Array beforr y values are tarnsformed for y scale
csv = getCsvArray(yMin,yMean,yMax,yLabel,yUnit,xLabel,xUnit,binBorders,nEntry,csvHeader,yMinRef,yMaxRef,yMeanRef);

% for yscal log transfer valiues to logscale
if strcmp(yscale,'log')
    yMin = log10(yMin);
    yMax = log10(yMax);
    yMean = log10(yMean);
            
    if ~isempty(yRef)
        yMinRef = log10(yMinRef);
        yMaxRef = log10(yMaxRef);
        yMeanRef = log10(yMeanRef);
        

    end
    
end



% create figure
ax = getReportFigure(WSettings,1,1,figureHandle);


% legend handle array for plots
lgh=[];

        
% exclude log(0)
ij = find(~isinf(yMin) & ~isinf(yMax));

lgh(2) = patch([binMedian(ij); binMedian(ij(end:-1:1))],[yMin(ij); yMax(ij(end:-1:1))],WSettings.colorVector(1,:),'linestyle','none','FaceAlpha',0.5,...
    'EdgeColor',WSettings.colorVector(1,:),'displayname',strtrim(sprintf('%s %s',legendEntries{1},rangeTxt)));

% Patch of referencep poulation
if ~isnan(yMinRef)
    
    lgh(4) = patch(binMedian([1 end end 1]),[yMinRef yMinRef yMaxRef yMaxRef],WSettings.colorVector(2,:),...
        'linestyle','none','FaceAlpha',0.5,...
        'EdgeColor',WSettings.colorVector(2,:),'displayname',strtrim(sprintf('%s %s',legendEntries{2},rangeTxt)));
end

% Plot mean
lgh(1)=plot(binMedian,yMean,'-','color',WSettings.colorVector(1,:),'linewidth',2,...
    'displayname',strtrim(sprintf('%s %s',legendEntries{1},legendTextMean)));

% Mean of referencep population
if ~isnan(yMinRef)
    lgh(3)=plot([binMedian(1) binMedian(end)],yMeanRef*[1 1],'-','color',WSettings.colorVector(2,:),'linewidth',2,...
        'displayname',strtrim(sprintf('%s %s',legendEntries{2},legendTextMean)));
end

% plot data


% set legend
jj = lgh>0;
legend(lgh(jj),get(lgh(jj),'displayname'),'location','northoutside');
    
% Set axes labels
xlabel(getLabelWithUnit(xLabel,xUnit));
ylabel(getLabelWithUnit(yLabel,yUnit));

set(ax,'xlim',[binMedian(1) binMedian(end)]);

if strcmp(yscale,'log')
    setLogarithmicYticks(ax);    
end


return

function csvArray = getCsvArray(yMin,yMean,yMax,yLabel,yUnit,xLabel,xUnit,binBorders,nEntry,csvHeader,yMinRef,yMaxRef,yMeanRef)

csvArray={};

csvArray{1,4} = getLabelWithUnit(yLabel,yUnit);
csvArray{1,1} = getLabelWithUnit(xLabel,xUnit);

csvArray(2,1:3+length(csvHeader)) = [{'min','max','N individuals'},csvHeader];

binLower = binBorders(1:end-1);
binUpper = binBorders(2:end);

csvArray(3:2+length(nEntry),:) = ...
    num2cell([binLower',binUpper',nEntry,yMin(2:end-1),yMean(2:end-1),yMax(2:end-1)]);

if ~isnan(yMinRef)
    iRow=size(csvArray,1)+2;
    csvArray{iRow,1}='comparsion values:';
    csvArray(iRow,4:6) = {yMinRef,yMeanRef,yMaxRef}; %#ok<NBRAK>
end
