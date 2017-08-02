function csv =  plotReportTimeProfile(WSettings,figureHandle,time,y,timeRef,yRef,timeUnit,yLabel,yUnit,legendEntries,yscale,lloq)
%PLOTREPORTTIMEPROFILE Plots the time profile of a population in comparison to a reference population
%
%  csv =  plotReportTimeProfile(WSettings,time,y,timeRef,yRef,timeUnit,yLabel,yUnit,yscale,LLOQ,figureHandle)
%
% Inputs: 
%       WSettings (structure)    definition of properties used in all
%                   workflow functions see GETDEFAULTWORKFLOWSETTINGS
%   figureHandle ( handle) handle of figure
%   time (double vector): time vector (in units of time_unit)
%   y  (double matrix (nIndividual x n timepoints)): time profile of population
%   timeRef (double vector): time vector of reference population (in units of time_unit)
%           if empty no reference population is drawn
%   yRef  (double matrix (nIndividual x n timepoints)): time profile of of reference population 
%           if empty no reference population is drawn
%   timeUnit (string)   text beside the X-axis
%   yLabel (string)  text beside the Y-axis
%   yUnit (string)   text beside the Y-axis
%   legendEntries (cell array of strings) {name of population for  legend,
%                       name of reference population for legend}
%   yscale (string): may be 'lin' or 'log' defines
%       number and Y-axes scale of the axes
%   lloq, Lower limit for quantification, if it is not nan it is plotted as line
%
% Output
%   csv (cellarray) table with numeric information to the plot

% Open Systems Pharmacology Suite;  http://forum.open-systems-pharmacology.org
% Date: 28-July-2017

% initialize XLS Array
csv = {}; 

% calculate percentiles of population
[yMin,yMean,yMax,legendTextMean,rangeTxt,csvHeader] = getPercentiles(WSettings,y);

% calculate percentiles of reference pop
if ~isempty(yRef) 
    [yMinRef,yMeanRef,yMaxRef] = getPercentiles(WSettings,yRef);
else
    yMinRef = nan;
end


% set csvArray
if ~isnan(yMin)   
    
    % description
    csv{1,1} = getLabelWithUnit(yLabel,yUnit);
    
    % time
    csv{2,1} = sprintf('Time [%s]',timeUnit);
    csv(3:2+length(time),1) = num2cell(time);
    
    % values
    csv(2,2:4) = csvHeader;
    csv(3:2+length(time),2:4) = num2cell([yMin,yMean,yMax]);
else
    
    % time
    csv{1,1} = sprintf('Time [%s]',timeUnit);
    csv(2:1+length(time),1) = num2cell(time);

    % description
    csv{1,2} = getLabelWithUnit(ylabel,yUnit);    
    csv(2:1+length(time),2) = num2cell(yMean);
    
end


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
    if ~isnan(lloq)
        lloq = log10(lloq);
    end
end

% create figure
ax = getReportFigure(WSettings,1,1,figureHandle);


% legend handle array for plots
lgh=[];
    
% if not single plot patch
if ~isnan(yMin)   
    
    % exclude log(0)
    ij = find(~isinf(yMin) & ~isinf(yMax)); 
    
    
    lgh(2) = patch([time(ij); time(ij(end:-1:1))],[yMin(ij); yMax(ij(end:-1:1))],WSettings.colorVector(1,:),'linestyle','none','FaceAlpha',0.5,...
        'EdgeColor',WSettings.colorVector(1,:),'displayname',strtrim(sprintf('%s %s',legendEntries{1},rangeTxt)));
end

% Patch of referencep poulation
if ~isnan(yMinRef)
    
    % exclude log(0)
    jj = ~isinf(yMin); 

    lgh(4) = patch([timeRef(jj); timeRef(end:-1:1)],[yMinRef(jj); yMaxRef(end:-1:1)],WSettings.colorVector(2,:),'linestyle','none','FaceAlpha',0.5,...
        'EdgeColor',WSettings.colorVector(2,:),'displayname',strtrim(sprintf('%s %s',legendEntries{2},rangeTxt)));
end

% Plot mean
lgh(1)=plot(time,yMean,'-','color',WSettings.colorVector(1,:),'linewidth',2,...
    'displayname',strtrim(sprintf('%s %s',legendEntries{1},legendTextMean)));

% Patch of referencep poulation
if ~isnan(yMinRef)
    lgh(4)=plot(timeRef,yMeanRef,'-','color',WSettings.colorVector(1,:),'linewidth',2,...
        'displayname',strtrim(sprintf('%s %s',legendEntries{2},legendTextMean)));
end
    
% plot data

% Plot LLOQ
if ~isnan(lloq)
    plot([time(1) time(end)],[1 1]*lloq,'k:','linewidth',2,'displayname','lower limit of quantification');
end

% set legend
legend(lgh,get(lgh,'displayname'),'location','northoutside');
    
% Set axes labels
xlabel(sprintf('Time [%s]',timeUnit));
ylabel(getLabelWithUnit(yLabel,yUnit));

setAxesScaling(ax,'timeUnit',timeUnit,'xlim',[time(1) time(end)]);
    
% set ytick to logscale if necessary
xl = get(ax,'xlim');
yl = get(ax,'ylim');

if strcmp(yscale,'log')
    yt = get(ax,'ytick');
    
    % get new integer ticks
    ytNew = intersect(yt,ceil(yt(1)):floor(yt(end)));
    
    % if none or only one add outer limits
    if  length(ytNew) <=1
        ytNew = [floor(yt(1)),ytNew,ceil(yt(end))];
        yl = [ytNew(1) ytNew(end)];
        set(ax,'ylim',yl)
    end
 
    set(ax,'ytick',ytNew,'yticklabel',10.^ytNew);
 
    % get ticks of yt
    ytTicks = 10.^(ytNew(1)).*[0.1:0.1:1];
    for iT = 1:length(ytNew)
        ytTicks = [ytTicks 10.^(ytNew(iT)).*[2:1:10]];
    end
    jj = log10(ytTicks) >= yt(1) &  log10(ytTicks) <= yt(end);
    ytTicks = log10(ytTicks(jj));
    
    for iT = 1:length(ytTicks)
        plot(xl(1)+[0 0.008*range(xl)],ytTicks(iT)*[1 1],'k-');
        plot(xl(2)-[0 0.008*range(xl)],ytTicks(iT)*[1 1],'k-');
    end
    
    
end
% plot box
plot(xl([1 1]),yl([1 2]),'k-');
plot(xl([2 2]),yl([1 2]),'k-');
plot(xl([1 2]),yl([1 1]),'k-');
plot(xl([1 2]),yl([2 2]),'k-');
    

return



function [yMin,yMean,yMax,legendTextMean,rangeTxt,csvHeader] = getPercentiles(WSettings,y)

if size(y,1) == 2 % only one individual
    yMin = nan;
    yMax = nan;
    yMean = y;
    legendTextMean = '';
    
    csvHeader = {};
else
    % exclude non processable individuals
    jjNonan = all(~isnan(y));
    yMin=prctile(y(:,jjNonan), WSettings.displayPercentiles(1),2);
    yMax=prctile(y(:,jjNonan), WSettings.displayPercentiles(end),2);
    switch WSettings.shadedAreaMeanType
        case 'median'
            yMean=median(y,2);
            legendTextMean='median';
        case 'geomean'
            yMean=geomean(y,2);
            legendTextMean='geo. mean';
        case 'mean'
            yMean=mean(y,2);
            legendTextMean='arith. mean';
    end

    % get text for range legend
    rangeTxt = sprintf('%s-%s percentile',...
        getPercentilePotenzText(WSettings.displayPercentiles(1)),...
        getPercentilePotenzText(WSettings.displayPercentiles(end)));

    csvHeader = {sprintf('%s precentile',getPercentilePotenzText(WSettings.displayPercentiles(1))),...
        legendTextMean,sprintf('%s precentile',getPercentilePotenzText(WSettings.displayPercentiles(end)))};
        
end
