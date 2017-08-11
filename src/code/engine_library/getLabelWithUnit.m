function yLabelFinal = getLabelWithUnit(yLabel,yUnit)
% GETLABELWITHUNIT add unit to label if unit is not empty
% 
%  yLabel = getLabelWithUnit(yLabel,yUnit)
%
%   yLabel (string) label without unit
%   yLabel (string) unit
%   yLabelFinal (string) label with unit


% Open Systems Pharmacology Suite;  http://open-systems-pharmacology.org


if ~isempty(yUnit)
    yLabelFinal = sprintf('%s [%s]',yLabel,yUnit);
else
    yLabelFinal = yLabel;
end
