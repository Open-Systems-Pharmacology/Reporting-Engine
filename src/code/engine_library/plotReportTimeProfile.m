function csv =  plotReportTimeProfile(WSettings,figureHandle,time,y,timeRef,yRef,DataTP,timeLabel,timeUnit,yLabel,yUnit,legendEntries,yscale,lloq)
%PLOTREPORTTIMEPROFILE Plots the time profile of a population in comparison to a reference population
%
% csv =  plotReportTimeProfile(WSettings,figureHandle,time,y,timeRef,yRef,DataTP,timeUnit,yLabel,yUnit,legendEntries,yscale,lloq)
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
%   DataTP (structure) dataprofile
%   timeLabel (string)   text beside the X-axis
%   timeUnit (string)   text beside the X-axis
%   yLimit (double vector 1 x 2)   limits of x axes
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

% Open Systems Pharmacology Suite;  http://open-systems-pharmacology.org

% calculate percentiles of population
[yMin,yMean,yMax,legendTextMean,rangeTxt,csvHeader] = getRangePlotPercentiles(WSettings,y');

% calculate percentiles of reference pop
if ~isempty(yRef) 
    [yMinRef,yMeanRef,yMaxRef] = getRangePlotPercentiles(WSettings,yRef');
else
    yMinRef = nan;
end

% get Csv Array beforr y values are tarnsformed for y scale
csv = getCsvArray(yMin,yMean,yMax,yLabel,yUnit,timeLabel,timeUnit,time,csvHeader);


% for yscal log transfer valiues to logscale
if strcmp(yscale,'log')
    yMin = log10(yMin);
    yMax = log10(yMax);
    yMean = log10(yMean);
        
    ylim = min([yMin(11:end), yMean(11:end)]);
    
    if ~isempty(yRef)
        yMinRef = log10(yMinRef);
        yMaxRef = log10(yMaxRef);
        yMeanRef = log10(yMeanRef);
        
        ylim = min([yMin(11:end), yMean(11:end),yMinRef(11:end),yMeanRef(11:end)]);

    end
    if ~isnan(lloq)
        lloq = log10(lloq);
    end
    
    if ~isempty(DataTP)
        for iInd = 1:length(DataTP)
            DataTP(iInd).y = log10(DataTP(iInd).y);
            jj = ~isinf(DataTP(iInd).y);
            ylim = min([ylim;DataTP(iInd).y(jj)]);
        end
    end
    
    ylim = ylim*0.8;
else
    ylim = 0;
    
end




% create figure
ax = getReportFigure(WSettings,1,1,figureHandle);


% legend handle array for plots
lgh=[];
    
% if not single plot patch
if ~isnan(yMin)   
    
    % exclude log(0)
    ij = find(~isinf(yMin) & ~isinf(yMax)); 
    
    
    lgh(2) = patch([time(ij), time(ij(end:-1:1))],[yMin(ij), yMax(ij(end:-1:1))],WSettings.colorVector(1,:),'linestyle','none','FaceAlpha',0.5,...
        'EdgeColor',WSettings.colorVector(1,:),'displayname',strtrim(sprintf('%s %s',legendEntries{1},rangeTxt)));
end

% Patch of referencep poulation
if ~isnan(yMinRef)
    
    % exclude log(0)
     ij = find(~isinf(yMinRef) & ~isinf(yMaxRef)); 

    lgh(4) = patch([timeRef(ij), timeRef(ij(end:-1:1))],[yMinRef(ij), yMaxRef(ij(end:-1:1))],WSettings.colorVector(2,:),'linestyle','none','FaceAlpha',0.5,...
        'EdgeColor',WSettings.colorVector(2,:),'displayname',strtrim(sprintf('%s %s',legendEntries{2},rangeTxt)));
end

% Plot mean
lgh(1)=plot(time,yMean,'-','color',WSettings.colorVector(1,:),'linewidth',2,...
    'displayname',strtrim(sprintf('%s %s',legendEntries{1},legendTextMean)));

% Mean of referencep population
if ~isnan(yMinRef)
    lgh(3)=plot(timeRef,yMeanRef,'-','color',WSettings.colorVector(2,:),'linewidth',2,...
        'displayname',strtrim(sprintf('%s %s',legendEntries{2},legendTextMean)));
end
    
% plot data
if ~isempty(DataTP)
    [col,mk] = getcolmarkForMap(WSettings.colormapData,length(DataTP));
    
    for iInd = 1:length(DataTP)
        lgh(5) = plot(DataTP(iInd).time,DataTP(iInd).y,mk(iInd),'color',col(iInd,:),'markerfacecolor',col(iInd,:),...
            'displayname',legendEntries{3});
        plot(DataTP(iInd).time,DataTP(iInd).lloq,mk(iInd),'color',col(iInd,:),'linewidth',2);
    end
end    

% Plot LLOQ
if ~isnan(lloq)
    plot([time(1) time(end)],[1 1]*lloq,'k:','linewidth',1,'displayname',legendEntries{end});
end

% set legend
jj = lgh>0;
legend(lgh(jj),get(lgh(jj),'displayname'),'location','northoutside');
    
% Set axes labels
xlabel(getLabelWithUnit(timeLabel,timeUnit));
ylabel(getLabelWithUnit(yLabel,yUnit));

setAxesScaling(ax,'timeUnit',timeUnit,'xlim',[time(1) time(end)]);
    
% set ytick to logscale if necessary
xl = get(ax,'xlim');
yl = get(ax,'ylim');
yl(1) = ylim;
set(ax,'ylim',yl);

if strcmp(yscale,'log')    
    setLogarithmicYticks(ax);
end
% plot box
plot(xl([1 1]),yl([1 2]),'k-');
plot(xl([2 2]),yl([1 2]),'k-');
plot(xl([1 2]),yl([1 1]),'k-');
plot(xl([1 2]),yl([2 2]),'k-');
    

return


function csv = getCsvArray(yMin,yMean,yMax,yLabel,yUnit,timeLabel,timeUnit,time,csvHeader)
csv = {}; 


% set csvArray
if ~isnan(yMin)   
    
    % description
    csv{1,1} = getLabelWithUnit(yLabel,yUnit);
    
    % time
    csv{2,1} = getLabelWithUnit(timeLabel,timeUnit);
    csv(3:2+length(time),1) = num2cell(time);
    
    % values
    csv(2,2:4) = csvHeader;
    csv(3:2+length(time),2:4) = num2cell([yMin',yMean',yMax']);
else
    
    % time
    csv{1,1} = getLabelWithUnit(timeLabel,timeUnit);
    csv(2:1+length(time),1) = num2cell(time);

    % description
    csv{1,2} = getLabelWithUnit(yLabel,yUnit);    
    csv(2:1+length(time),2) = num2cell(yMean);
    
end
