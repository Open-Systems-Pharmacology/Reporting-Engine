function FP = plotSensListMostSensitive(WSettings,FP,sens,sortSumSensIx,iCut,parameterList,marker,colorVector,lgtxtPrctl,popLabels)
% PLOTSENSLISTMOSTSENSITIVE plots all sensitivity above cutoff value
%
% FP = plotSensListMostSensitive(WSettings,FP,sensPop,sortSumSensIx,iCut,parameterList,marker,colorVector,lgtxtPrctl,popLabels)
%
% Inputs:
%       WSettings (structure)    definition of properties used in all
%                   workflow functions see GETDEFAULTWORKFLOWSETTINGS
%       FP (ReportFigure) object to manage figure print
%       sens (cellarray of structures) with sensitivity informtaion
%               cell array correpsonds to OutputList, structur entries to
%               PK Parameter list
%   sortSumSensIx (doubleVector) index vector of sorted sensitivities
%   iCut (double) number of parameters above cut
%   parameterList (cellarray) list of senstivity parameter
%   marker (char) marker for different percentile individuals
%   colorVector (matrix) colarmatrix for different population 
%   lgtxtPrctl (cellarray of strings) legend text for differnt percentile individuals
%   popLabels (cellarray of strings) legend text for differnt populations
%

% Open Systems Pharmacology Suite;  http://open-systems-pharmacology.org




% create figure
ax = getReportFigure(WSettings,1,1,FP.figureHandle,'axes_position',[0.3965 0.11 0.5085 0.815]);

% get number of sensitivity vectors
nPop=size(sens,2);
nPrc=size(sens,1);

% plot sensitivity
lgh = [];
for iPop=1:nPop
    for iPrc=1:nPrc
        
        if nPop*nPrc >1
            tmp = plot([sens{iPop,iPrc}(sortSumSensIx(1:iCut)).slope],[1:iCut],...
                marker(iPrc),'color',colorVector(iPop,:),'markerfacecolor',colorVector(iPop,:)); %#ok<NBRAK>
        else
            for iPar = 1:iCut
                patch([0 1 1 0].*sens{iPop,iPrc}(sortSumSensIx(iPar)).slope,iPar + [-0.4 -0.4 0.4 0.4],colorVector(iPop,:),'facealpha',0.5)
            end
        end
        
        % set legend
        if ~isempty(popLabels)
            lgh = addToLegendPopulationSensitivity(lgh,tmp,popLabels,lgtxtPrctl,iPop,iPrc);
        end

        plot([[sens{iPop,iPrc}(sortSumSensIx(1:iCut)).slopeCILower];[sens{iPop,iPrc}(sortSumSensIx(1:iCut)).slopeCIUpper]],...
            repmat([1:iCut],2,1),...
            'color',colorVector(iPop,:),'linewidth',1); %#ok<NBRAK>


        
    end
end

% set axes scaling
xl = get(ax,'Xlim');
xl = [-1 1]*max(abs(xl));
yl = [0,iCut+1];
setAxesScaling(ax,'xlim',xl,'xscale','lin',...
    'yscale','lin','ylim',yl);
set(ax,'YDir','reverse')
grid on
xlabel('sensitivity')
% set ytxt
set(ax,'ytick',1:iCut)
set(ax,'yticklabel',[])
set(ax,'Xlim',xl);

for iPar=1:iCut
    text(xl(1)-0.05*(xl(2)-xl(1)),iPar,parameterList{iPar,1},'horizontalalignment','right','interpreter','none','fontsize',8)
end

if ~isempty(popLabels)
    legend(lgh,get(lgh,'displayname'),'location','northoutside','fontsize',8)
end

return
