function [csv] = plotReportBoxwhisker(WSettings,figureHandle,y,yLabel,yUnit,yscale,xLabel,popPKValues,popPKReference)
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
%   popPKValues: (cellarray of strings) description of groups on x axes
%   xLabel: (cellarray of strings) description of groups on x axes
% Output
%   csv (cellarray) table with numeric information to the plot

% Open Systems Pharmacology Suite;  http://open-systems-pharmacology.org


% calculate percentiles
outlier=[];
pct = prctile(y,WSettings.displayPercentiles);

if size(pct,1)==1
    pct = pct';
end

for iSet=1:size(y,2)
    
    if WSettings.boxwhiskerWithExtrema
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
ax = getReportFigure(WSettings,1,1,figureHandle,'figureformat','landscape');

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
if WSettings.boxwhiskerWithExtrema
    plot(outlier(:,1),outlier(:,2),'o','color',WSettings.colorVector(1,:),'linewidth',2,'markersize',6);
end

if ~isnan(popPKValues)
    posX = x(end)+0.75;
    errorbar(posX,popPKValues(2),popPKValues(2)-popPKValues(1),popPKValues(3)-popPKValues(2),'o','color','k',...
        'linewidth',2,'markersize',6);
    xlim(end) = posX+widthBox;
end

% set layout
setAxesScaling(ax,'xlim',xlim,'yscale',yscale);
    
% grid
set(ax,'Ygrid','on')
set(ax,'YMinorGrid','off')
    
% Ylabel
ylabel(getLabelWithUnit(yLabel,yUnit));
    
% Xlabels
if ~isfield(WSettings,'boxwhiskerXlabelRotation') || WSettings.boxwhiskerXlabelRotation ==0
    set(ax,'xtick',x,'xticklabel',xLabel)
else
    set(ax,'xtick',1:length(xLabel),'xticklabel',[]);
    yl = get(ax,'ylim');
    switch yscale
        case 'lin'
                ytext = yl(1)-0.05*(yl(2)-yl(1));
        case 'log'
                ytext = exp(log(yl(1))-0.05*(log(yl(2))-log(yl(1))));
        otherwise
            error('unknown scale')
    end
    if ~isnan(popPKValues)
        text(posX,ytext,popPKReference,'rotation',WSettings.boxwhiskerXlabelRotation,'horizontalAlignment','right');
        set(ax,'xtick',[1:length(xLabel) posX],'xticklabel',[]);
    end
    for iX = 1:length(xLabel) 
        tmp = xLabel{iX};
        if length(tmp)>1
            ix =  strfind(tmp,'with');
            if isempty(ix)
                ix =  strfind(tmp,' ');
            end
            [~,ixPos] = min(abs(ix-length(tmp)/2));
            tmp = {strtrim(tmp(1:ix(ixPos)-1)),strtrim(tmp(ix(ixPos):end))};
        end
        text(x(iX),ytext,tmp,'rotation',WSettings.boxwhiskerXlabelRotation,'horizontalAlignment','right');
    end
    p = get(ax,'position');
    set(ax,'position',p+[0 0.2 0 -0.2]);
end

% no legend
legend(ax,'off')

% csv    

csv{1,1} = getLabelWithUnit(yLabel,yUnit);
% percentiles
for iPrct=1:length(WSettings.displayPercentiles)
    csv{2,1+iPrct} = sprintf('%s prctl.',getPercentilePotenzText(WSettings.displayPercentiles(iPrct)));
end

% get current size
iRow = size(csv,1);

    
% header for differnt sets
winY = (iRow+1:iRow+size(y,2));
csv(winY,1) = xLabel;

% percentiles
csv(winY,2:(1+size(pct,1))) = num2cell(pct)';

% add additional numbers
csv{2,end+1} = 'arith. mean';
csv(winY,end) = num2cell(nanmean(y));
csv{2,end+1} = 'std';
csv(winY,end) = num2cell(nanstd(y));
csv{2,end+1} = 'geo. mean';
csv(winY,end) = num2cell(exp(nanmean(log(y))));
csv{2,end+1} = 'geo. CV';
s = nanstd(log(y));
csv(winY,end) = num2cell(sqrt(exp(s.^2)-1));



return


