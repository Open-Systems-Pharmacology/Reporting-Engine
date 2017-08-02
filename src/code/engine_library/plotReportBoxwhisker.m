function [csv] = plotReportBoxwhisker(WSettings,figureHandle,y,yLabel,yUnit,yscale,xLabel)
%PLOTREPORTBOXWHISKER creates box whisker plot
%
%   plotReportBoxwhisker(X,Y,Y_ref,...
%    Y_label,Y_unit,ylimit,yscale,withExtrema,BW_Layout,...
%    figure_handle,colorVector_patch,colorVector_line)
%
%   y  (double matrix) (nIndividual x n Population)): parameters to plot
%   yLabel (string)  text beside the Y-axis
%   yUnit (string)   text beside the Y-axis
%   yscale (string): may be 'lin' or 'log' scale of the Y axes
%   xLabel: (cellarray of strings) description of groups on x axes
% Output
%   csv (cellarray) table with numeric information to the plot


% Open Systems Pharmacology Suite;  http://forum.open-systems-pharmacology.org
% Date: 1-Aug-2017


% calculate percentiles
outlier=[];
pct = prctile(y,WSettings.displayPercentiles);

if size(pct,1)==1
    pct = pct';
end

for iSet=1:size(y,2)
    
    if WSettings.withExtrema
        jjOut = y(:,iSet) < pct(1,iSet) | y(:,iSet) > pct(end,iSet);
        outlier=[outlier;[repmat(iSet,sum(jjOut),1), y(jjOut,iSet)]]; %#ok<AGROW>
    end
end

% X coordinates
x = 1:size(y,2);
widthBox = 0.2;
xL = x - widthBox;
xU = x + widthBox;
xlim=[xL(1)-widthBox xU(end) + widthBox];

% find columns whith data 
ji_plot=find(~all(isnan(pct)));


% do the Plot
ax = getReportFigure(WSettings,1,1,figureHandle);

% horizontal lines
for iY = 1:5
    for iX = ji_plot
        plot([xL(iX),xU(iX)],[pct(iY,iX) pct(iY,iX)],'-','color',WSettings.colorVector(1,:),'linewidth',2);
        hold on;
    end
end
    
% side of box
for iX = ji_plot
    plot([xL(iX),xL(iX)],[pct(2,iX) pct(4,iX)],'-','color',WSettings.colorVector(1,:),'linewidth',2);
    plot([xU(iX),xU(iX)],[pct(2,iX) pct(4,iX)],'-','color',WSettings.colorVector(1,:),'linewidth',2);
    plot([x(iX),x(iX)],[pct(1,iX) pct(2,iX)],'-','color',WSettings.colorVector(1,:),'linewidth',2);
    plot([x(iX),x(iX)],[pct(4,iX) pct(5,iX)],'-','color',WSettings.colorVector(1,:),'linewidth',2);
end

% outlier
if WSettings.withExtrema
    plot(outlier(:,1),outlier(:,2),'o','color',WSettings.colorVector(1,:),'linewidth',2,'markersize',6);
end

% set layout
setAxesScaling(ax,'xlim',xlim,'yscale',yscale);
    
% grid
set(ax,'Ygrid','on')
set(ax,'YMinorGrid','off')
    
% Ylabel
ylabel(getLabelWithUnit(yLabel,yUnit));
    
% Xlabels
set(ax,'xtick',x,'xticklabel',xLabel)

% csv    

csv{1,1} = getLabelWithUnit(yLabel,yUnit);
% percentiles
for iPrct=1:length(WSettings.displayPercentiles)
    csv{2,1+iPrct} = sprintf('%s precentile',getPercentilePotenzText(WSettings.displayPercentiles(iPrct)));
end

% get current size
iRow = size(csv,1);

    
% heder for differnt sets
csv((iRow+1:iRow+size(y,2)),1) = xLabel;

% percentiles
csv((iRow+1):(iRow+size(y,2)),2:(1+size(pct,1))) = num2cell(pct);

return


