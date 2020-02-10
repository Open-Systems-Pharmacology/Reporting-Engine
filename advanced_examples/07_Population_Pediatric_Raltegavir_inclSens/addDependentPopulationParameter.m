function [parPaths,parValues,success] = addDependentPopulationParameter(WSettings,parPaths,parValues)
%ADDDEPENDENTPOPULATIONPARAMETER add additional population parameter 
%   for plotting which can be caluclated out of the existing ones
%
%  [parPaths,parValues] = addDependentPopulationParameter(WSettings,parPaths,parValues)
%
% Inputs:
%   WSettings   (structure) containing global settings see GETDEFAULTWORKFLOWSETTINGS
%   parPaths (cellarray of strings) pathnames of parameter
%   parValues ( double array) values of parameter
%   success (boolean) if false error has occured.
%
% Outputs:
%   parPaths (cellarray of strings) pathnames of parameter
%   parValues ( double array) values of parameter

% Open Systems Pharmacology Suite;  http://open-systems-pharmacology.org

success = true;

[parPaths,parValues,successTmp] = addBSA(WSettings,parPaths,parValues);
success = success & successTmp;

%% add GFR
[parPaths,parValues,successTmp] = addGFR(parPaths,parValues);
success = success & successTmp;

%% add GFR Norm
[parPaths,parValues,successTmp] = addGFRBSANorm(parPaths,parValues);
success = success & successTmp;

%% add fraction unbounds
[parPaths,parValues,successTmp] = addFU(parPaths,parValues);
success = success & successTmp;


%% add unnecessary ontogeny factors
parPaths{end+1} = 'ABCG2 (BCRP)|Ontogeny factor';
parValues(:,end+1) = zeros(size(parValues,1),1);
parPaths{end+1} = 'ABCG2 (BCRP)|Ontogeny factor GI';
parValues(:,end+1) = zeros(size(parValues,1),1);
parPaths{end+1} = 'ABCB1 (P-gp)|Ontogeny factor';
parValues(:,end+1) = zeros(size(parValues,1),1);
parPaths{end+1} = 'ABCB1 (P-gp)|Ontogeny factor GI';
parValues(:,end+1) = zeros(size(parValues,1),1);

return

function [parPaths,parValues,success] = addFU(parPaths,parValues) 

success = true;
nPar = length(parPaths);
    
try
    % %% fraction unbounds
    [albuminOntgenyFactor,success] = getPopulationParameter(parPaths,parValues,'Organism|Ontogeny factor (albumin)',success);
        
    [agpOntogenyFactor,success] = getPopulationParameter(parPaths,parValues,'Organism|Ontogeny factor (alpha1-acid glycoprotein)',success);

    if success
        % get name of drugs with a positive fraction unbound
        [~,desc] = existsParameter('*Fraction unbound (plasma, reference value)',1,'parametertype','readonly');
        tmpStruct = cell2struct(desc(2:end,:),desc(1,:),2);
        jj = [tmpStruct.Value] > 0;
        for iJ = find(jj)
            fu_ref = tmpStruct(iJ).Value;
            tmp = regexp(tmpStruct(iJ).Path,'\|','split');
            drug = tmp{2};
            bindingPartner = getParameter(strjoin([tmp(1:2),{'Plasma protein binding partner'}],'|'),1,'parametertype','readonly');
            switch bindingPartner
                case 1
                    FU =  1 ./ (1 + albuminOntgenyFactor .* (1 - fu_ref) ./ fu_ref);
                case 2
                    FU =  1 ./ (1 + agpOntogenyFactor .* (1 - fu_ref) ./ fu_ref);
                otherwise
                    FU = fu_ref ;
            end
            
            parPaths{end+1} = strjoin({'DependentParameter',drug,'Fraction unbound (plasma)'},'|'); %#ok<AGROW>
            parValues(:,end+1) = FU.*100; %#ok<AGROW>
            
        end
        writeToReportLog('INFO','fraction unbound was added to population',false);
    end
catch
    success = false;
end

if ~success
    writeToReportLog('WARNING','It was not possible to add fraction unbound to population',false);
    parPaths = parPaths(1:nPar);
    parValues = parValues(:,1:nPar);
end

return

function [value,success] = getPopulationParameter(popPaths,popValues,parPath,success)

% initialize return value
value = nan;

% if it failed before this call return
if ~success
    return
end

% check if it is in the population
jj = strcmp(parPath,popPaths);
if any(jj)
   value = popValues(:,jj);
else
    % get parameter from simulation file
    [ise,desc] = existsParameter(['*|' parPath],1,'parametertype','readonly');
    if ~ise
        success = false;
        writeToReportLog('WARNING',sprintf('Missing parameter %s for calculation of dependent parameter!',parPath),false);
    elseif size(desc,1)>2
        success = false;
        writeToReportLog('WARNING',sprintf('Ambiguous parameter %s for calculation of dependent parameter!',parPath),false);
    else
        value = getParameter(['*|' parPath],1,'parametertype','readonly');
        isFormula = getParameter(['*|' parPath],1,'parametertype','readonly','property','isFormula');
        if isFormula
            writeToReportLog('WARNING',sprintf(['Parameter %s for calculation of dependent parameter is a formula parameter!'...
                ' Be careful may be different for other individuals'],parPath),false);
        end
    end
end

return


function [parPaths,parValues,success] = addBSA(WSettings,parPaths,parValues)

success = true;

%% add BSA if not available
if ~any(strcmp(parPaths,'Organism|BSA'))
    
    nPar = length(parPaths);
    
    try
        [weight,success] = getPopulationParameter(parPaths,parValues,'Organism|Weight',success);
        [height,success] = getPopulationParameter(parPaths,parValues,'Organism|Height',success);
        
        if success
            parPaths{end+1} = 'Organism|BSA';
            parValues(:,end+1) = calculateBodySurfaceArea(WSettings,weight,height);
            writeToReportLog('INFO','Organism|BSA was added to population',false);
        end
    catch
        success = false;
    end

    if ~success
        writeToReportLog('WARNING','It was not possible to add Organism|BSA to population',false);
        parPaths = parPaths(1:nPar);
        parValues = parValues(:,1:nPar);
    end

end

return


function [parPaths,parValues,success] = addGFR(parPaths,parValues)

success = true;
nPar = length(parPaths);

try
    % renal Aging factor
    [Age,success] = getPopulationParameter(parPaths,parValues,'Organism|Age',success);
    
    [AgeOnset,success] = getPopulationParameter(parPaths,parValues,'Organism|Kidney|Age of aging onset',success);
    
    [MaxDecreaseFactor,success] = getPopulationParameter(parPaths,parValues,'Organism|Kidney|Maximal decreasing rate factor',success);% 10.9
    
    [HillAging,success] = getPopulationParameter(parPaths,parValues,'Organism|Kidney|Hill coefficient for aging GFR',success);% 1.5
    
    [AgingHalfTime,success] = getPopulationParameter(parPaths,parValues,'Organism|Kidney|Aging half-time',success);

    if success
        renalAgingScalingFactor = (Age <= AgeOnset).* 1 + ...
            (Age > AgeOnset).* (1 - MaxDecreaseFactor.*(Age-AgeOnset).^HillAging./(AgingHalfTime.^HillAging+(Age-AgeOnset).^HillAging));
    end

    
    % GFRspecific
    
    [GestationalAge,success] = getPopulationParameter(parPaths,parValues,'Organism|Gestational age',success);

        
    [Hill,success] = getPopulationParameter(parPaths,parValues,'Organism|Kidney|Hill coefficient for GFR',success);
    
    [TM50,success] = getPopulationParameter(parPaths,parValues,'Organism|Kidney|TM50 for GFR',success);
    
    [GFRPremat,success] = getPopulationParameter(parPaths,parValues,'Organism|Kidney|fGFRpremat',success);%0.258;
    
    [GFRmat,success] = getPopulationParameter(parPaths,parValues,'Organism|Kidney|GFRmat',success);
    
    [vKidStd,success] = getPopulationParameter(parPaths,parValues,'Organism|Kidney|Volume (standard kidney)',success);%0.44;
    
    
    if success
        PMA = Age* 365.25/7 + GestationalAge;
        GFRspec = (PMA.^Hill./(TM50.^Hill+PMA.^Hill).*(1-GFRPremat)+GFRPremat) .* GFRmat ./ vKidStd .* renalAgingScalingFactor;
    end
    
    
    %GFR
    [KidneyVolume,success] = getPopulationParameter(parPaths,parValues,'Organism|Kidney|Volume',success);

    if success
        GFR = GFRspec.*KidneyVolume;
        
        parPaths{end+1} = 'DependentParameter|GFR';
        parValues(:,end+1) = GFR;
        
        writeToReportLog('INFO','GFR parameter were added to population',false);
    end

catch
    success = false;
end

if ~success
    writeToReportLog('WARNING','It was not possible to add GFR parameter to population',false);
    parPaths = parPaths(1:nPar);
    parValues = parValues(:,1:nPar);
end

return

function [parPaths,parValues,success] = addGFRBSANorm(parPaths,parValues)

success = true;

nPar = length(parPaths);
try
    % GFR BSA Norm
    [BSA,success] = getPopulationParameter(parPaths,parValues,'Organism|BSA',success);
    
    [GFR,success] = getPopulationParameter(parPaths,parValues,'DependentParameter|GFR',success);

    if success
        GFRbsaNorm = GFR ./ BSA * 1.73 .*1000 ; % [mL/min/1.73m2]
        
        parPaths{end+1} = 'DependentParameter|GFRbsaNorm1dot73';
        parValues(:,end+1) = GFRbsaNorm;
        
        writeToReportLog('INFO','GFR norm parameter was added to population',false);
    end

catch
    success = false;
end

if ~success
    writeToReportLog('WARNING','It was not possible to add GFR norm parameter to population',false);
    parPaths = parPaths(1:nPar);
    parValues = parValues(:,1:nPar);
end
return
