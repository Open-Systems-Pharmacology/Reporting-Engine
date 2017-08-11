function runPopulationSimulation(WSettings,PopRunSet)
% RUNPOPULATIONSIMULATION  runs a population simulation defined within a PopRunSet 
%
% runPopulationSimulation(WSettings,PopRunSet)
% 
% Inputs:
%   WSettings  structure containing global settings see GETDEFAULTWORKFLOWSETTINGS
%   PopRunSet  defines properties and corresponding files of a simulation
%           see GENERATEWORKFLOWINPUTFORPOPULATIONSIMULATION


% Open Systems Pharmacology Suite;  http://open-systems-pharmacology.org



writeToLog(sprintf('Start Simulation of %s',PopRunSet.name),WSettings.logfile,true,false);

% read population 
load(fullfile('tmp',PopRunSet.name,'pop.mat'),'parPaths','parValues');

% load OutpulList
load(fullfile('tmp',PopRunSet.name,'outputList.mat'),'OutputList');

% get individual ids
nInd = size(parValues,1); %#ok<NODEF>
individualIdVector = parValues(:,strcmp(parPaths,'IndividualId')); 

% prepare bunches
bunches = unique([1:WSettings.nIndPerSimResult:nInd nInd+1]);
nBunch = length(bunches)-1;

for iBunch = 1:nBunch
    
    indVector = bunches(iBunch):(bunches(iBunch+1)-1);
    
    % do the simulation
    [SimResult] = generateSimResult(WSettings,PopRunSet.xml,nan,parPaths,parValues(indVector,:),...
        {OutputList.pathID},length(indVector),individualIdVector(indVector));
        
    % export results to  new result file
    exportResult(SimResult,PopRunSet.name,iBunch);
    
    % save as temporary file
    save(fullfile('tmp',PopRunSet.name,sprintf('simResult_%d.mat',iBunch)),'SimResult');
    
end
writeToLog(sprintf('Simulation finished \n'),WSettings.logfile,true,false);

return

function exportResult(SimResult,simulationName,iBunch)

% get name of resultfile
if iBunch==1
    resultfile = fullfile('simulations',[simulationName '-Results.csv']);
else
    resultfile = fullfile('simulations',sprintf('%s-%d-Results.csv',simulationName,iBunch));
end
% check if siluation directory already exist
if ~exist('simulations','dir')
    mkdir('simulations')
end


% get dimensions
nT =  length(SimResult.time);
nInd = length(SimResult.individualIdVector);
nO = length(SimResult.outputPathList);


% construct result cell array
csvResult(1,:) = {'IndividualId','Time [min]'};
for iO = 1:nO;
    csvResult{1,iO+2} = sprintf('%s [%s]',SimResult.outputPathList{iO},SimResult.outputUnit{iO});
end
% unit type
csvResult(2,:) = repmat({'double'},1,nO+2);

% values
csvResult{3,1} = reshape(repmat(SimResult.individualIdVector',nT,1),nT*nInd,1);
csvResult{3,2} = repmat(SimResult.time',nInd,1);
for iO = 1:length(SimResult.outputPathList);
    csvResult{3,2+iO} = reshape(SimResult.values{iO},nInd*nT,1);
end

% write result
writetab(resultfile,csvResult,';',0,0,1,0);

return
