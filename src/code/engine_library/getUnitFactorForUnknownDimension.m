function [unitFactor,success] = getUnitFactorForUnknownDimension(Settings,unit,targetUnit,MW)
% GETUNITFACTORFORUNKNOWNDIMENSION getfactor to convert units, 
%
% Inputs: 
%      - Settings (structure)    definition of properties used in all
%                   workflow functions see GETDEFAULTWORKFLOWSETTINGS
%      - targetUnit (string) unit to convert internal unit
%      - MW (double) moelcular weigth
%  Outputs:
%      - unitFactor (double)  factor to convert interanl unit to target unit
%      - success (boolean) if false, it was not possible to determine unit conversion factor

% Open Systems Pharmacology Suite;  http://forum.open-systems-pharmacology.org
% Date: 26-July-2017

persistent unitList;
persistent unitList_dimensionList;

if isempty(unitList)
    [unitList,unitList_dimensionList]=iniUnitList(0);
end

% initialize return values
unitFactor = nan;
success = true;



% check if taregt unit exist, and get index of dimension
ix_dim_display = find(cellfun(@(x) any(strcmp(x,targetUnit)),{unitList.unit_txt}));
if isempty(ix_dim_display)
    writeToLog(sprintf('ERROR: unit "%s" for seems to be no default OSPSuite unit',targetUnit),Settings.logfile,true,false);
    success = false;
    return
end
    
% get index of dimension for internal unit
ix_dim_internal = find(cellfun(@(x) any(strcmp(x,unit)),{unitList.unit_txt}));

% target unit and internal unit should have a common dimension
ix_dim = intersect(ix_dim_display,ix_dim_internal);
if isempty(ix_dim)
    writeToLog(sprintf('ERROR: For unit "%s", there is no common dimension with display unit "%s"',unit,targetUnit),Settings.logfile,true,false);
    success = false;
    return
end

% get Unitfactor
% for concnetration molweight is als needed
if strcmp(unitList_dimensionList{ix_dim(1)},'Concentration')
    if ~exist('MW','var')
        error('molweight is needed');
    end        
    unitFactor = getUnitFactor(unit,targetUnit,unitList_dimensionList{ix_dim(1)},'MW',MW*1e9);

else    
    unitFactor = getUnitFactor(unit,targetUnit,unitList_dimensionList{ix_dim(1)});
end

return