function BSA = calculateBodySurfaceArea(WSettings,weight,height)
%CALCULATEBODYSURFACEAREA calculates the body surface area
%
% BSA = calculateBodySurfaceArea(WSettings,weight,height)
%
%Inputs: 
%   WSettings (structure)    definition of properties used in all
%                   workflow functions see GETDEFAULTWORKFLOWSETTINGS
%  weight (double vector) body weight of individuals [kg]
%  height (double vector) height of individuals [dm]

% Open Systems Pharmacology Suite;  http://forum.open-systems-pharmacology.org
% Date: 31-July-2017

switch WSettings.BSA_calculationmethode

    case 'DuBois'        
        BSA = (0.007184.*(height.*10).^0.725 .*weight.^0.425).*100;
        
    case 'Haycock'
        BSA = (0.024265.* (height.*10).^0.3964.* weight.^0.5378).*100;
    
end

return



