function  lgh = addToLegendPopulationSensitivity(lgh,tmp,popLabels,lgtxtPrctl,iPop,iPrc)
% ADDTOLEGENDPOULATIONSENSITIVITY add legend entires to population sensitivity plots
%
%  lgh = addToLegendPoulationSensitivity(lgh,tmp,popLabels,lgtxtPrctl,iPop,iPrc)
%
%  lgh (vector) handles to lineseries objects, will be incremented
%  tmp current lineseries handle
%  popLabels ( cellarray of strings) description of populations
%  lgtxtPrctl ( cellarray of strings) description of percentile individual
%  iPop (double) number of population
% iPrc (double) number of percentile individual


% Open Systems Pharmacology Suite;  http://open-systems-pharmacology.org

% set legend
if iPop ==1
    lgh(end+1) = tmp;
    set(tmp,'displayname',sprintf('%s %s',popLabels{iPop},lgtxtPrctl{iPrc}));
elseif iPrc ==2
    lgh(end+1) = tmp;
    set(tmp,'displayname',sprintf('%s',popLabels{iPop}));
end

return