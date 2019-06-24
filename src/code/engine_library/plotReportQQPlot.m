function  plotReportQQPlot(WSettings,figureHandle,res)
% PLOTREPORTQQPLOT creates qqPlot for residuals
%
% [csvArray] = plotReportQQPlot(WSettings,figureHandle,y,yRef,xLabel,legendEntries)
% 
% Inputs: 
%       WSettings (structure)    definition of properties used in all
%                   workflow functions see GETDEFAULTWORKFLOWSETTINGS
%       res  ( double vector nInd x 1) property to plot
%      sourceIndex (double vector) (nInd x 1) index of source population
%
% Output
%   csv (cellarray) table with numeric information to the plot

% Open Systems Pharmacology Suite;  http://open-systems-pharmacology.org

% do the Plot
ax = getReportFigure(WSettings,1,1,figureHandle);

lgh = qqplot(res);

set(lgh(1),'LineWidth',2);
set(lgh(2),'LineWidth',2);
set(lgh(3),'LineWidth',2);

ylabel('Quantiles of residuals');
title('');
grid on;

xl = get(ax,'xlim');
xl = max(abs(xl)).*[-1 1];
set(ax,'xlim',xl);

return