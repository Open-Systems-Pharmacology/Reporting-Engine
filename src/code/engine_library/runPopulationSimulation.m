function SimResult = runPopulationSimulation(Settings,PopRunSet)
% RUNPOPULATIONSIMULATION  runs a population simulation defined within a PopRunSet 
% see GETDEFAULTPOPRUNSET
%
% Inputs:
%   Settings  structure containing global settings see GETDEFAULTWORKFLOWSETTINGS
%   PopRunSet  defines properties and corresponding files of a simulation
%           see GETDEFAULTPOPRUNSET


% Open Systems Pharmacology Suite;  http://forum.open-systems-pharmacology.org
% Date: 14-July-2017


writeToLog(sprintf('Start Simulation of %s',PopRunSet.name),Settings.logfile,true,false);

% read population 
load(fullfile('tmp',PopRunSet.name,'pop.mat'),'parPaths','parValues');

% get individual ids
nInd = size(parValues,1); %#ok<NODEF>
individualIdVector = parValues(:,strcmp(parPaths,'IndividualId')); 

% do the simulation
[SimResult] = doIndividualsimulation(Settings,PopRunSet.name,PopRunSet.xml,parPaths,parValues,{PopRunSet.OutputList.pathID},nInd,individualIdVector);

% export results to  new result file
exportResult(SimResult);

% save as temporary file
save(fullfile('tmp',PopRunSet.name,'simResult.mat'),'SimResult');

writeToLog(sprintf('Simulation finished \n'),Settings.logfile,true,false);

return

function exportResult(SimResult)

% get name of resultfile
resultfile =fullfile('Simulations',[SimResult.name '-Results.csv']); 
% check if siluation directory already exist
if ~exist('Simulations','dir')
    mkdir('Simulations')
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
    csvResult{3,2+iO} = reshape(SimResult.values{iO}',nInd*nT,1);
end

% write result
writetab(resultfile,csvResult,';',0,0,1,0);

return

function SimResult = doIndividualsimulation(Settings,name,xml,parPaths,parValues,outputPathList,nInd,individualIdVector)

% initialize model
initStruct=[];
for iPar=1:length(parPaths)
    initStruct=initParameter(initStruct,[ '*|' parPaths{iPar}],'always','throwWarningIfNotExisting',false);
end
initSimulation(xml,initStruct);

% Only one simulation is intialized so the simulation index is always 1
simulationIndex = 1;

% getSimulationTime For the import all results must have the same
% timepoints, which happens to be not the case for simulations outside
% PK-Sim (???), but the differences are only small % todo @Juri discuss
% this comment with Juri
simTimeNorm = getSimulationTime(simulationIndex);

% initialize rowindex of parameter for speedy setting
rowIndex = [];
% initialize rowindex of output for speedy fetching
rowIndex_S = nan(size(outputPathList));

% initialize result array
SimResult.name = name;
SimResult.time = simTimeNorm;
SimResult.values = cell(length(outputPathList),1);
for iO = 1:length(outputPathList)
    SimResult.values{iO} = nan(nInd,length(simTimeNorm));
end
SimResult.individualIdVector = individualIdVector;
SimResult.outputPathList = outputPathList;
SimResult.outputUnit = cell(size(outputPathList));



% Loop On Individuals
for iInd=1:nInd
    
    % set physiology parameters
    % first time get row index
    if isempty(rowIndex)
        
        rowIndex=nan(size(parPaths));
        for iPar=1:length(parPaths)
            if existsParameter(['*|' parPaths{iPar}],simulationIndex);
                 tmp = setParameter(parValues(iInd,iPar),[ '*|' parPaths{iPar}],simulationIndex);
                 
                rowIndex(iPar) =tmp(1);
                % chek if it is a unique parameter 
                if length(tmp)>1
                    
                    writeToLog(sprintf('WARNING %s is a parameter of the population csv, which correspond to more then one simulation parameter',...
                        parPaths{iPar}),Settings.logfile,true,false);
                    
                    for iTmp = 2:length(tmp)
                        indx = length(rowIndex)+1;
                        rowIndex(indx) = tmp(iTmp);
                        parValues(:,indx) = parValues(:,iPar);
                    end
                end
                
            else
                % check if is a demographic parameter which is not an xml
                % parameter, otherwise,  set a warning
                if ~ismember(parPaths{iPar},{'IndividualId','Gender','RaceIndex','Population Name','Organism|Weight','Organism|BMI'})                
                    writeToLog(sprintf('WARNING %s is parameter of the population csv, but not available as parameter in the xml', ...
                        parPaths{iPar}),Settings.logfile,true,false);
                end
            end                
        end
    else
        % with row index use fast methode
        for iPar=find(~isnan(rowIndex))
            setParameter(parValues(iInd,iPar),'',simulationIndex,...
                'speedy','variable',rowIndex(iPar));
        end
    end
    
    
    success=processSimulation;
    
    if success
        % get result
        for iO=1:length(outputPathList)
            [simTime,tmpValues,rowIndex_S(iO)] = getSimulationResult(['*|' outputPathList{iO}],simulationIndex,'rowindex',rowIndex_S(iO));
            % map outputs to defined time raster
            if ~all(simTime - simTimeNorm==0)
                SimResult.values{iO}(iInd,:) = interp1(simTime,tmpValues,simTimeNorm,'linear','extrap');
            else
                SimResult.values{iO}(iInd,:) = tmpValues;
            end
        end
    else
        writeToLog(sprintf('WARNING solver error, at individual %d',iInd),Settings.logfile,true,false);
    end
    
end

% add Units of outputs to output structure for export
for iO=1:length(outputPathList)
    if existsObserver(['*|' outputPathList{iO}],1)
        SimResult.outputUnit{iO} = getObserverFormula(['*|' outputPathList{iO}],1,'property','Unit');
    else    
        SimResult.outputUnit{iO} = getSpeciesInitialValue(['*|' outputPathList{iO}],1,'parametertype','readonly','property','Unit');
    end
    
end


return

