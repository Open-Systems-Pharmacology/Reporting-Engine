function calculatePopulationPKParameter(Settings,PopRunSet)
% CALCULATEPOPULATIONPKPARAMETER calculates PK-Parameter for simulations defined within a PopRunSet 
% see GETDEFAULTPOPRUNSET
%
% Inputs:
%   Settings  structure containing global settings see GETDEFAULTWORKFLOWSETTINGS
%   PopRunSet  defines properties and corresponding files of a simulation
%           see GETDEFAULTPOPRUNSET
 

% Open Systems Pharmacology Suite;  http://forum.open-systems-pharmacology.org
% Date: 14-July-2017

writeToLog(sprintf('Calculate PK Parameter of %s',PopRunSet.name),Settings.logfile,true,false);


% load population 
load(fullfile('tmp',PopRunSet.name,'pop.mat'),'parPaths','parValues');

% reduce to demographic
[jj,ix] = ismember({'IndividualId','Organism|Weight','Organism|Height','Organism|BMI'},parPaths); %#ok<NODEF>
parPaths = parPaths(ix(jj));
parValues = parValues(:,ix(jj)); %#ok<NODEF>


% load applicationProtocol
load(fullfile('tmp',PopRunSet.name,'applicationProtocol.mat'),'ApplicationProtocol','isValid');

% load time profiles
SimResult = readPopulationResultfile(PopRunSet.name);

% initialize PK List
PKPList = cell(length(SimResult.outputPathList) ,1);

% decide on the PKParameter function and calculate PKParameter
% if given take the projectspecific function
if ~isempty(PopRunSet.calculatePKParameter_fh)

    for iO = 1:length(SimResult.outputPathList)
        PKPList{iO} = feval(PopRunSet.calculatePKParameter_fh,Settings,SimResult.time,SimResult.values{iO},parPaths,parValues,ApplicationProtocol);
    end

    
elseif isValid
    
    jj_BW = strcmp('Organism|Weight',parPaths);
    weight =  parValues(:,jj_BW);
    
    for iO = 1:length(SimResult.outputPathList)
        PKPList{iO} = calculatePKParameterForApplicationProtocol(Settings,SimResult.time,SimResult.values{iO},weight,ApplicationProtocol);
    end

else            
    error('application protocol could not be interpreted. Please use project specific function for PK Parameter calculation');
end

% export PKParameter
exportPKParameter(PKPList,SimResult);

% save as temporary file
save(fullfile('tmp',PopRunSet.name,'pKPList.mat'),'PKPList');


writeToLog(sprintf('Calculation finished \n'),Settings.logfile,true,false);

return

function exportPKParameter(PKPList,SimResult)
% get name of resultfile
resultfile =fullfile('Simulations',[SimResult.name '-PK-Analyses.csv']); 

% check if siluation directory already exist
if ~exist('Simulations','dir')
    mkdir('Simulations')
end

% construct result cell array
csvResult(1,:) = {'IndividualId','Quantity Path','Parameter','Value','Unit'};
% unit type
csvResult(2,:) = {'double','string','string','double','string'};

% get Dimensions
nInd = length(SimResult.individualIdVector);
nO = length(PKPList);

% get values
for iO = 1:nO
    
    PKParameter = PKPList{iO};
    
    for iPKP = 1:length(PKParameter)
        
        if iO*iPKP ==1
            csvResult{3,1} = SimResult.individualIdVector;
            csvResult{3,2} = repmat(SimResult.outputPathList(iO),nInd,1);
            csvResult{3,3} = repmat({PKParameter(iPKP).name},nInd,1);
            csvResult{3,4} = PKParameter(iPKP).value;
            csvResult{3,5} = repmat({PKParameter(iPKP).unit},nInd,1);
        else
            csvResult{3,1} = [csvResult{3,1}; SimResult.individualIdVector];
            csvResult{3,2} = [csvResult{3,2}; repmat(SimResult.outputPathList(iO),nInd,1)];
            csvResult{3,3} = [csvResult{3,3}; repmat({PKParameter(iPKP).name},nInd,1)];
            csvResult{3,4} = [csvResult{3,4}; PKParameter(iPKP).value];
            csvResult{3,5} = [csvResult{3,5}; repmat({PKParameter(iPKP).unit},nInd,1)];
        end        
        
    end
end

% write result
writetab(resultfile,csvResult,';',0,0,1,0);


return