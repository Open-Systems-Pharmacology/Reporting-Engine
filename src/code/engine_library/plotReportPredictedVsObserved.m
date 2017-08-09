function  plotReportPredictedVsObserved(WSettings,figureHandle,DataTP,yLabel,yUnit,legendEntries,yscale,lloq)
% PLOTREPORTPREDICTEDVSOBSERVED plot predcited vs observed
%
% csv =  plotReportPredictedVsObserved(WSettings,figureHandle,DataTP,yLabel,yUnit,legendEntries,yscale,lloq)
%
% Inputs: 
%       WSettings (structure)    definition of properties used in all
%                   workflow functions see GETDEFAULTWORKFLOWSETTINGS
%   figureHandle ( handle) handle of figure
%   DataTP (structure) dataprofile
%   yLimit (double vector 1 x 2)   limits of x axes
%   yLabel (string)  text beside the Y-axis
%   yUnit (string)   text beside the Y-axis
%   legendEntries (cell array of strings) {name of population for  legend,
%                       name of reference population for legend}
%   yscale (string): may be 'lin' or 'log' defines
%       number and Y-axes scale of the axes
%   lloq, Lower limit for quantification, if it is not nan it is plotted as line
%

% Open Systems Pharmacology Suite;  http://open-systems-pharmacology.org


% get maxima and minima
switch yscale
    case 'lin'
        yMin = min(min(cellfun(@min,{DataTP.y})),min(cellfun(@min,{DataTP.predicted})))*0.8;
        yMax = max(max(cellfun(@max,{DataTP.y})),max(cellfun(@max,{DataTP.predicted})))*1.2;
    case 'log'
        yMin = inf;
        for iData = 1:length(DataTP)
            jj = DataTP(iData).y >0 & DataTP(iData).predicted>0;
            yMin =  min([DataTP(iData).y(jj);DataTP(iData).predicted(jj);yMin]);
        end
        yMin = yMin/2;
        yMax = max(max(cellfun(@max,{DataTP.y})),max(cellfun(@max,{DataTP.predicted})))*2;
end

% create figure
ax = getReportFigure(WSettings,1,1,figureHandle);

% line of identity
lgh(2) = plot([yMin yMax],[yMin yMax],'k-','displayname',legendEntries{2});

[col,mk] = getcolmarkForMap(WSettings.colormapData,length(DataTP));
    
for iInd = 1:length(DataTP)
    lgh(1) = plot(DataTP(iInd).y,DataTP(iInd).predicted,mk(iInd),'color',col(iInd,:),'markerfacecolor',col(iInd,:),'displayname',legendEntries{1});
    plot(DataTP(iInd).lloq,DataTP(iInd).predicted,mk(iInd),'color',col(iInd,:),'linewidth',2);
end

set(ax,'xscale',yscale,'yscale',yscale,'DataAspectRatio',[1 1 1]);

set(ax,'ylim',[yMin yMax],'xlim',[yMin yMax]);

% Plot LLOQ
if ~isnan(lloq)
    lgh(3) = plot([1 1]*lloq,[yMin yMax],'k:','linewidth',1,'displayname',legendEntries{3});
end

% Set axes labels
xlabel({'observed',getLabelWithUnit(yLabel,yUnit)});
ylabel({'predicted',getLabelWithUnit(yLabel,yUnit)});


legend(lgh,get(lgh,'displayname'),'location','eastoutside');

return
