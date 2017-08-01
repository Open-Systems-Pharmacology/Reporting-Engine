function [unitFactor,success] = getUnitFactorForUnknownDimension(WSettings,unit,targetUnit,MW)
% GETUNITFACTORFORUNKNOWNDIMENSION getfactor to convert units, 
%
% Inputs: 
%      - WSettings (structure)    definition of properties used in all
%                   workflow functions see GETDEFAULTWORKFLOWSETTINGS
%      - targetUnit (string) unit to convert internal unit
%      - MW (double) moelcular weigth
%  Outputs:
%      - unitFactor (double)  factor to convert interanl unit to target unit
%      - success (boolean) if false, it was not possible to determine unit conversion factor

% Open Systems Pharmacology Suite;  http://open-systems-pharmacology.org


persistent unitList;
persistent unitListDimensionList;

if isempty(unitList)
    [unitList,unitListDimensionList] = iniUnitList(0);
end

% initialize return values
unitFactor = nan;
success = true;



% check if taregt unit exist, and get index of dimension
ixDimDisplay = find(cellfun(@(x) any(strcmp(x,targetUnit)),{unitList.unit_txt}));
if isempty(ixDimDisplay)
    writeToLog(sprintf('ERROR: unit "%s" for seems to be no default OSPSuite unit',targetUnit),WSettings.logfile,true,false);
    success = false;
    return
end
    
% get index of dimension for internal unit
ixDimInternal = find(cellfun(@(x) any(strcmp(x,unit)),{unitList.unit_txt}));

% target unit and internal unit should have a common dimension
ixDim = intersect(ixDimDisplay,ixDimInternal);
if isempty(ixDim)    
    
    % merge molar and mass AUC
    jj = ismember({'AUC (molar)','AUC (mass)'},unitListDimensionList([ixDimDisplay,ixDimInternal]));
    if all(jj) && exist('MW','var')
        switch unitListDimensionList{ixDimDisplay}
            case 'AUC (mass)'
                unitFactorMolar = getUnitFactor(unit,'µmol*min/l','AUC (molar)');
                unitFactorMass = getUnitFactor('µg*min/l',targetUnit,'AUC (mass)');
                unitFactor = unitFactorMolar.*MW.*unitFactorMass;
            case 'AUC (molar)'
                unitFactorMass = getUnitFactor(targetUnit,'µg*min/l','AUC (mass)');
                unitFactorMolar = getUnitFactor('µmol*min/l',unit,'AUC (molar)');
                unitFactor = unitFactorMolar./MW.*unitFactorMass;
        end
        return
    end
    
    writeToLog(sprintf('ERROR: For unit "%s", there is no common dimension with display unit "%s"',unit,targetUnit),WSettings.logfile,true,false);
    success = false;
    return
end

% get Unitfactor
% for concnetration molweight is als needed
if strcmp(unitListDimensionList{ixDim(1)},'Concentration')
    if ~exist('MW','var')
        error('molweight is needed');
    end        
    unitFactor = getUnitFactor(unit,targetUnit,unitListDimensionList{ixDim(1)},'MW',MW*1e9);

else    
    unitFactor = getUnitFactor(unit,targetUnit,unitListDimensionList{ixDim(1)});
end

return
