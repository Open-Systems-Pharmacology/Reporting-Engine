function [csv,goodness] =  plotReportTimeProfile(WSettings,figureHandle,SimTL,DataTP,timeLabel,timeUnit,yLabel,yUnit,legendEntries,yscale,lloq,...
    VPCresult,iVPCFlag)
%PLOTREPORTTIMEPROFILE Plots the time profile of a population in comparison to a reference population
%
% csv =  plotReportTimeProfile(WSettings,figureHandle,time,y,timeRef,yRef,DataTP,timeUnit,yLabel,yUnit,legendEntries,yscale,lloq)
%
% Inputs: 
%       WSettings (structure)    definition of properties used in all
%                   workflow functions see GETDEFAULTWORKFLOWSETTINGS
%   figureHandle ( handle) handle of figure
%   SimTL (structure) with simulated results 
%       .time (double vector): time vector (in units of time_unit)
%       .y  (double matrix (nIndividual x n timepoints)): time profile of population
%       .timeRef (double vector): time vector of reference population (in units of time_unit)
%           if empty no reference population is drawn
%       .yRef  (double matrix (nIndividual x n timepoints)): time profile of of reference population 
%           if empty no reference population is drawn
%       .meanModel  (double matrix (1 x n timepoints)): time profile of population
%   DataTP (structure) dataprofile
%   timeLabel (string)   text beside the X-axis
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

% Open Systems Pharmacology Suite;  http://open-systems-pharmacology.org

% initialize goodness
goodness = {'',nan}; %#ok<NASGU>


% get Range
if ~exist('VPCresult','var') || isempty(VPCresult)
    % calculate percentiles of population
    [yMin,yMean,yMax,legendTextMean,rangeTxt,csvHeader] = getRangePlotPercentiles(WSettings,SimTL.y');
    
    % calculate percentiles of reference pop
    if ~isempty(SimTL.yRef)
        [yMinRef,yMeanRef,yMaxRef] = getRangePlotPercentiles(WSettings,SimTL.yRef');
    else
        yMinRef = nan;
    end
    
else
    yMean = SimTL.y';
    yMax = VPCresult.yMax(:,iVPCFlag)';
    yMin = VPCresult.yMin(:,iVPCFlag)';
    legendTextMean = VPCresult.legendTextMean;
    rangeTxt = VPCresult.rangeTxt;
    csvHeader = VPCresult.csvHeader;
    
     yMinRef = nan;
     
end

% get Csv Array beforr y values are transformed for y scale
csv = getCsvArray(WSettings,yMin,yMean,yMax,yLabel,yUnit,timeLabel,timeUnit,SimTL.time,csvHeader);

% reshape time vector
SimTL.time = reshape(SimTL.time,size(yMean));
if ~isempty(SimTL.yRef)
    SimTL.timeRef = reshape(SimTL.timeRef,size(yMeanRef));
end
 
% for yscal log transfer values to logscale ( otherwise problems with patch
% are possible
if strcmp(yscale,'log')
    yMin = log10(yMin);
    yMax = log10(yMax);
    yMean = log10(yMean);
        
    if ~isnan(yMin)
        ylim = prctile(yMin,5);
    else
        ylim = prctile(yMean,5);
    end
    
    if ~isempty(SimTL.yRef)
        yMinRef = log10(yMinRef);
        yMaxRef = log10(yMaxRef);
        yMeanRef = log10(yMeanRef);
        
        ylim = min([yMin(11:end), yMean(11:end),yMinRef(11:end),yMeanRef(11:end)]);

    end
    
    if ~isempty(SimTL.yMeanModel)
        SimTL.yMeanModel = log10(SimTL.yMeanModel);
    end
    
    if ~isnan(lloq)
        lloq = log10(lloq);
    end
    
    if ~isempty(DataTP)
        ylimD = inf;
        for iInd = 1:length(DataTP)
            DataTP(iInd).y = log10(DataTP(iInd).y);
            DataTP(iInd).lloq = log10(DataTP(iInd).lloq);
            jj = ~isinf(DataTP(iInd).y);
            ylimD = min([ylimD;DataTP(iInd).y(jj)]);
        end
        
        ylim = min(ylim,ylimD);
    end
    
    ylim = ylim/1.2^sign(ylim);
    
else    
    ylim = 0;    
    
end

% create figure
ax = getReportFigure(WSettings,1,1,figureHandle,'figureformat','landscape');


% Patch of reference poulation
% Plot before main population, so it is in background
if ~isnan(yMinRef)
    
    % exclude log(0)
     ij = find(~isinf(yMinRef) & ~isinf(yMaxRef)); 

    lgh(4) = patch([SimTL.timeRef(ij), SimTL.timeRef(ij(end:-1:1))],[yMinRef(ij), yMaxRef(ij(end:-1:1))],WSettings.colorVector(2,:),...
        'linestyle','-','LineWidth',0.5,'FaceAlpha',WSettings.FaceAlpha,...
        'EdgeColor',WSettings.colorVector(2,:),'displayname',strtrim(strrep(legendEntries{2},'<xxx>',rangeTxt)));
    
    % Mean of referencep population 
    lgh(3)=plot(SimTL.timeRef,yMeanRef,'-','color',WSettings.colorVector(2,:),'linewidth',2,...
        'displayname',strtrim(strrep(legendEntries{2},'<xxx>',legendTextMean)));

end



% if not single plot patch
if ~isnan(yMin)   
    
    % exclude values below ylim
    jj = yMin < ylim;
    yMin(jj) = ylim;
    ij = find(yMax>=ylim); 
    
    
    lgh(2) = patch([SimTL.time(ij), SimTL.time(ij(end:-1:1))],[yMin(ij), yMax(ij(end:-1:1))],WSettings.colorVector(1,:),...
        'linestyle','-','LineWidth',0.5,'FaceAlpha',WSettings.FaceAlpha,...
        'EdgeColor',WSettings.colorVector(1,:),'displayname',strtrim(strrep(legendEntries{1},'<xxx>',rangeTxt)));
end


% Plot mean
lgh(1)=plot(SimTL.time,yMean,'-','color',WSettings.colorVector(1,:),'linewidth',2,...
    'displayname',strtrim(strrep(legendEntries{1},'<xxx>',legendTextMean)));

% plotMeanModel
if ~isempty(SimTL.yMeanModel)
    lgh(end+1) = plot(SimTL.time,SimTL.yMeanModel,'-.','color',WSettings.colorVector(1,:)*0.7,'linewidth',2,...
        'displayname',strtrim(legendEntries{4}));
    
    % resort legend
    lgh = [lgh(1:2) lgh(end) lgh(3:end-1)];
end

    
% plot data
if ~isempty(DataTP)
    [col,mk] = getcolmarkForMap(WSettings.colormapData,length(DataTP));    
    iLghData = length(lgh)+1;
    
    for iInd = 1:length(DataTP)
        if ~isempty(DataTP(iInd).time)
            lgh(iLghData) = plot(DataTP(iInd).time,DataTP(iInd).y,mk(iInd),'color',col(iInd,:),'markerfacecolor',col(iInd,:),...
                'displayname',legendEntries{3});
            plot(DataTP(iInd).time,DataTP(iInd).lloq,mk(iInd),'color',col(iInd,:),'linewidth',2,'displayname',legendEntries{3});
        end
    end
end    

% Plot LLOQ
if ~isnan(lloq)
    lgh(end+1) = plot([SimTL.time(1) SimTL.time(end)],[1 1]*lloq,'k-','linewidth',2,'displayname',legendEntries{end});
end

    
% Set axes labels
xlabel(getLabelWithUnit(timeLabel,timeUnit));
ylabel(getLabelWithUnit(yLabel,yUnit));

% set ytick to logscale if necessary
yl = get(ax,'ylim');
yl(1) = ylim;

setAxesScaling(ax,'timeUnit',timeUnit,'xlim',[SimTL.time(1) SimTL.time(end)]);
set(ax,'ylim',yl);

if strcmp(yscale,'log')    
    setLogarithmicYticks(ax);
end

% set legend
jj = isgraphics(lgh);
legend(lgh(jj),get(lgh(jj),'displayname'),'location','northoutside');


% getGoodness
if ~isempty(DataTP) && any(~isnan(yMin))
    nIn = 0;
    nAll = 0;
    for iInd = 1:length(DataTP)
        jjNAN = ~isnan(DataTP(iInd).y);
        nAll = sum(jjNAN) + nAll;
        % get predicted min and max
        DataMin = interp1(SimTL.time,yMin,DataTP(iInd).time);
        DataMax = interp1(SimTL.time,yMax,DataTP(iInd).time);
        jj = jjNAN & DataTP(iInd).y >= DataMin & DataTP(iInd).y <= DataMax;
        nIn = sum(jj) + nIn;
    end
    
    p =  (WSettings.displayPercentiles(5)-WSettings.displayPercentiles(1))/100;
    
    nMax = binoinv(0.95,nAll,p);
    nMin = binoinv(0.05,nAll,p);
    
    success = nMin <= nIn & nIn <= nMax;
    
    goodness = {sprintf('%d of %d datapoints within %g%% range of simulation, expected are values between %d %d',nIn,nAll,p*100,nMin,nMax),success};
    
    
else
    goodness = {'',nan};
end
return


function csv = getCsvArray(WSettings,yMin,yMean,yMax,yLabel,yUnit,timeLabel,timeUnit,time,csvHeader)
csv = {}; 


% set csvArray
if ~isnan(yMin)   
    
    % description
    csv{1,1} = getLabelWithUnit(yLabel,yUnit);
    
    % time
    csv{2,1} = getLabelWithUnit(timeLabel,timeUnit);
    csv(3:2+length(time),1) = cellfun(@(x) sprintf(WSettings.csvFormat,x),num2cell(time),'uniformoutput',false);
    
    % values
    csv(2,2:4) = csvHeader;
    csv(3:2+length(time),2:4) = cellfun(@(x) sprintf(WSettings.csvFormat,x),num2cell([yMin',yMean',yMax']),'uniformoutput',false);
else
    
    % time
    csv{1,1} = getLabelWithUnit(timeLabel,timeUnit);
    csv(2:1+length(time),1) = cellfun(@(x) sprintf(WSettings.csvFormat,x),num2cell(time),'uniformoutput',false);

    % description
    csv{1,2} = getLabelWithUnit(yLabel,yUnit);    
    csv(2:1+length(time),2) = cellfun(@(x) sprintf(WSettings.csvFormat,x),num2cell(yMean),'uniformoutput',false);
    
end
