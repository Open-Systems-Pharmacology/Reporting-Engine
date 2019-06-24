function calculatePopulationPKParameter(WSettings,PopRunSet)
% CALCULATEPOPULATIONPKPARAMETER calculates PK-Parameter for simulations defined within a PopRunSet 
% see GENERATEWORKFLOWINPUTFORPOPULATIONSIMULATION
%
% calculatePopulationPKParameter(WSettings,PopRunSet)
% 
% Inputs:
%   WSettings  structure containing global settings see GETDEFAULTWORKFLOWSETTINGS
%   PopRunSet  defines properties and corresponding files of a simulation
%           see GENERATEWORKFLOWINPUTFORPOPULATIONSIMULATION
 

% Open Systems Pharmacology Suite;  http://open-systems-pharmacology.org

try
    writeToReportLog('INFO',sprintf('Calculate PK Parameter of %s',PopRunSet.name),false);
    
    
    % load population
    load(fullfile('tmp',PopRunSet.name,'pop.mat'),'parPaths','parValues');
    
    % reduce to demographic
    [jj,ix] = ismember(WSettings.demographicParameters,parPaths); %#ok<NODEF>
    parPaths = parPaths(ix(jj));
    parValues = parValues(:,ix(jj)); %#ok<NODEF>
    
    % start the caluclation
    [PKPList,individualIdVector,pathList] = calculatesPKParameterList(WSettings,PopRunSet.name,PopRunSet.calculatePKParameterFh,parPaths,parValues);
    
    % export PKParameter
    exportPKParameter(WSettings,PopRunSet.name,PKPList,individualIdVector,pathList);
    
    % save as temporary file
    save(fullfile('tmp',PopRunSet.name,'pKPList.mat'),'PKPList');
    
    
    writeToReportLog('INFO',sprintf('Calculation finished \n'),false);

catch exception
        
    save(sprintf('exception_%s.mat',datestr(now,'ddmmyy_hhMM')),'exception')
    writeToReportLog('ERROR',exception.message,false,exception);
    writeToReportLog('INFO',sprintf('Calculation finished with error \n'),false);
        
end

return

function exportPKParameter(WSettings,simulationName,PKPList,individualIdVector,outputPathList)

% check if siluation directory already exist
if ~exist(fullfile(cd,'simulations'),'dir')
    mkdir('simulations')
end

% prepare bunches
nInd = length(individualIdVector);
bunches = unique([1:WSettings.nIndPerSimResult:nInd nInd+1]);
nBunch = length(bunches)-1;

for iBunch = 1:nBunch

    
    indVector = bunches(iBunch):(bunches(iBunch+1)-1);

    % get name of resultfile
    if iBunch==1
        resultfile =fullfile('simulations',[simulationName '-PK-Analyses.csv']);
    else
        resultfile = fullfile('simulations',sprintf('%s-%d-PK-Analyses.csv',simulationName,iBunch));
    end
    
    
    % get Dimensions
    nInd = length(indVector);
    nO = length(PKPList);

    % get values
    for iO = 1:nO
        
        PKParameter = PKPList{iO};
        
        for iPKP = 1:length(PKParameter)
            
                IndividualId = individualIdVector(indVector);
                QuantityPath = repmat(outputPathList(iO),nInd,1);
                Parameter = repmat({PKParameter(iPKP).name},nInd,1);
                Value = PKParameter(iPKP).value(indVector)';
                Unit = repmat({PKParameter(iPKP).unit},nInd,1);
            if iO*iPKP ==1
                tTmp = table(IndividualId,QuantityPath,Parameter,Value,Unit);
            else
                tTmp = [tTmp;table(IndividualId,QuantityPath,Parameter,Value,Unit)]; %#ok<AGROW>
            end
            
        end
    end

    % write result
    writetable(tTmp,resultfile,'fileType','text','delimiter',';');
end

return
