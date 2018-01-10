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

try

    writeToReportLog('INFO',sprintf('Start Simulation of %s',PopRunSet.name),false);
    
    % check if siluation directory already exist
    if ~exist(fullfile(cd,'simulations'),'dir')
        mkdir('simulations')
    end
    
    
    % load OutpulList
    load(fullfile('tmp',PopRunSet.name,'outputList.mat'),'OutputList');

    % generate meanModel Result
    resultfile = fullfile('simulations',[PopRunSet.name '-MeanModel.csv']);
    if ~WSettings.restart || ~exist(fullfile(cd,resultfile),'file')
        SimResult = generateSimResult(WSettings,PopRunSet,nan,{},[],{OutputList.pathID},1,-1); 
        save(fullfile('tmp',PopRunSet.name,'meanModel_simResult_1.mat'),'SimResult');
        
        % export results to  new result file
        exportResult(SimResult,resultfile);
    end
    
    % read population
    load(fullfile('tmp',PopRunSet.name,'pop.mat'),'parPaths','parValues');
    
    
    % get individual ids
    nInd = size(parValues,1); %#ok<NODEF>
    individualIdVector = parValues(:,strcmp(parPaths,'IndividualId')); %#ok<IDISVAR>
    
    % prepare bunches
    bunches = unique([1:WSettings.nIndPerSimResult:nInd nInd+1]);
    nBunch = length(bunches)-1;
    
    for iBunch = 1:nBunch
        
        % get name of resultfile
        if iBunch==1
            resultfile = fullfile('simulations',[PopRunSet.name '-Results.csv']);
        else
            resultfile = fullfile('simulations',sprintf('%s-%d-Results.csv',PopRunSet.name,iBunch));
        end
        % simulate only if not restart and file does not exist
        if WSettings.restart && exist(fullfile(cd,resultfile),'file')
            writeToReportLog('INFO',sprintf('simulation for %s part %d already exist \n',PopRunSet.name,iBunch),false);
        else
            
            indVector = bunches(iBunch):(bunches(iBunch+1)-1);
        

            % do the simulation
            [SimResult] = generateSimResult(WSettings,PopRunSet,nan,parPaths,parValues(indVector,:),...
                {OutputList.pathID},length(indVector),individualIdVector(indVector));
            
            % export results to  new result file
            exportResult(SimResult,resultfile);
            
            % save as temporary file
            save(fullfile('tmp',PopRunSet.name,sprintf('simResult_%d.mat',iBunch)),'SimResult');
        end
        
    end
    writeToReportLog('INFO',sprintf('Simulation finished \n'),false);

catch exception
        
    save(sprintf('exception_%s.mat',datestr(now,'ddmmyy_hhMM')),'exception');
    writeToReportLog('ERROR',exception.message,false);
    writeToReportLog('INFO',sprintf('Simulation finished with error \n'),false);
        
end


return

function exportResult(SimResult,resultfile)


% get dimensions
nT =  length(SimResult.time);
nInd = length(SimResult.individualIdVector);
nO = length(SimResult.outputPathList);


% construct result cell array
csvResult(1,:) = {'IndividualId','Time [min]'};
for iO = 1:nO
    csvResult{1,iO+2} = sprintf('%s [%s]',SimResult.outputPathList{iO},SimResult.outputUnit{iO});
end

% values
csvResult(2:(nT*nInd+1),1) = num2cell(reshape(repmat(SimResult.individualIdVector',nT,1),nT*nInd,1));
csvResult(2:(nT*nInd+1),2) = num2cell(repmat(SimResult.time',nInd,1));
for iO = 1:length(SimResult.outputPathList)
    csvResult(2:(nT*nInd+1),2+iO) = num2cell(reshape(SimResult.values{iO},nInd*nT,1));
end

% write result
writeTabCellArray(csvResult,resultfile);

return
