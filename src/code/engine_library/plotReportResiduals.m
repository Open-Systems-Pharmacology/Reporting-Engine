function   [goodness] = plotReportResiduals(WSettings,figureHandle,DataTP,timeLabel,timeUnit,yLabel,yUnit,yscale,xAxesFlag,legendEntries)
% PLOTREPORTRESIDUALS plots residuals vs time or predicted y
%
% plotReportResiduals(WSettings,figureHandle,DataTP,xAxesFlag,yLabel,yUnit,timeLabel,Def.timeDisplayUnit,legendEntries,yscale,lloq)
%
% Inputs: 
%       WSettings (structure)    definition of properties used in all
%                   workflow functions see GETDEFAULTWORKFLOWSETTINGS
%   figureHandle ( handle) handle of figure
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
%   xAxesFlag (string) 'vsTime' or 'vsY' defines what is plotted on x axes
%   lloq, Lower limit for quantification, if it is not nan it is plotted as line


% create figure
ax = getReportFigure(WSettings,1,1,figureHandle,'figureformat','landscape');

% get properties of xAxes
switch xAxesFlag
    case 'vsTime'
        xField = 'time';
        xLabelTxt = getLabelWithUnit(timeLabel,timeUnit);
        xscale = 'lin';
    case    'vsY'
        xField = 'predicted';
        xLabelTxt = {'predicted',getLabelWithUnit(yLabel,yUnit)};
        xscale = yscale;
    otherwise
        error('unknown flag');

end



[col,mk] = getcolmarkForMap(WSettings.colormapData,length(DataTP));

% initialize vectors for ggodnessfit
gX = [];
gY = [];

for iInd = 1:length(DataTP)
    if ~isempty(DataTP(iInd).time)
    
        % calculate residulas according yscale % ToDo Normalisation
        switch yscale
            case 'lin'
                res = DataTP(iInd).y-DataTP(iInd).predicted;
            case 'log'
                res = log(DataTP(iInd).y)-log(DataTP(iInd).predicted);
            otherwise
                error('unknown scale');
                
        end
        
        lgh(1) = plot(DataTP(iInd).(xField),res,mk(iInd),'color',col(iInd,:),'markerfacecolor',col(iInd,:),'displayname',legendEntries{1});
        %     plot(DataTP(iInd).(xField),resLloq,mk(iInd),'color',col(iInd,:),'linewidth',2);
        
        % construct vectors for goodness fit
        gX = [gX; DataTP(iInd).(xField)]; %#ok<AGROW>
        gY = [gY; res]; %#ok<AGROW>
    end

end

set(ax,'xscale',xscale);

xl = get(ax,'xlim');
yl = get(ax,'ylim');
set(ax,'ylim',max(abs(yl))*[-1 1]);

plot(xl,[0 0],'k-','linewidth',2);

xlabel(xLabelTxt);
 switch yscale
     case 'lin'
         ylabel({'residual','(measured - predicted)'});
     case 'log'
         ylabel({'residual','(log(measured) - log(predicted))'});
     otherwise
         error('unknown scale');

 end
 
 
 legend(lgh,get(lgh,'displayname'),'location','northoutside');

 
 % getGoodness
 if strcmp(xscale,'log')
     gX = log(gX);
 end
 
 [b,bint] = regress(gY,[ones(size(gX)), gX]);
 % slope and intercept may be 0:
 success = all(sign(bint(:,1))<0) & all(sign(bint(:,2))>0);
    
 goodness = {sprintf('residuals fit a line with intercept %.3g [%.3g-%.3g] and slope %.3g [%.3g-%.3g]',b(1),bint(1,1),bint(2,2),b(1),bint(2,1),bint(2,2)),success};
    
 
 
 return
