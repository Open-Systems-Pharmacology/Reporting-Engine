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

    % replace (s) to s at the end of the unit
    ij = strfind(yUnit,'(s)');
    if ~isempty(ij) && ij == length(yUnit)-2
        yUnit = [yUnit(1:end-3) 's'];
    end
    
    yLabelFinal = sprintf('%s [%s]',yLabel,yUnit);
else
    yLabelFinal = yLabel;
end

% set first letter to capital
yLabelFinal = [upper(yLabelFinal(1)),yLabelFinal(2:end)]; 
