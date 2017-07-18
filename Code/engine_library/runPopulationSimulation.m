function SimResult = runPopulationSimulation(Settings,PopRunSet)
% RUNPOPULATIONSIMULATION  runs a population simulation defined within a PopRunSet 
% see GETDEFAULTPOPRUNSET
%
% Inputs:
%   Settings  structure containing global settings see GETDEFAULTWORKFLOWSETTINGS
%   PopRunSet  defines properties and corresponding files of a simulation
%           see GETDEFAULTPOPRUNSET


% Open Systems Pharmacology Suite;  support@systems-biology.com
% Date: 14-July-2017


writeToLog(sprintf('Start Simulation of %s',PopRunSet.name),Settings.logfile,true,false);

% read population csv
[parPathes,parValues] = readPopulationCSV(PopRunSet.pop);
nInd = size(parValues,1);
individualIdVector = parValues(:,strcmp(parPathes,'IndividualId'));

% add studyDesign if one is given
if ~isempty(PopRunSet.studyDesign)
    [parPathes,parValues] = readStudyDesign(Settings,PopRunSet.studyDesign,parPathes,parValues);
end

% read outputs
outputList = readOutputSelection(PopRunSet.outputList);

% do the simulation
[SimResult] = doIndividualsimulation(Settings,PopRunSet.name,PopRunSet.xml,parPathes,parValues,outputList,nInd,individualIdVector);

% export results to  new result file
exportResult(SimResult);

writeToLog(sprintf('Simulation finished \n'),Settings.logfile,true,false);

return

function exportResult(SimResult)

% get name of resultfile
resultfile =fullfile('Simulations',[SimResult.name '-Results.csv']); 
% check if siluation directory already exist
if ~exist('Simulations','dir')
    mkdir('Simulations')
end

% create header for Results
header='IndividualId;Time [min]';
for iO=1:length(SimResult.outputList);
    header=[header  ';' sprintf('%s [%s]',SimResult.outputList{iO},SimResult.outputUnit{iO})]; %#ok<AGROW>
end


% construct numeric result array
nT =  length(SimResult.time);
nInd = length(SimResult.individualIdVector);
R(:,1) = reshape(repmat(SimResult.individualIdVector',nT,1),nT*nInd,1);
R(:,2) = repmat(SimResult.time',nInd,1);
for iO = 1:length(SimResult.outputList);
    R(:,2+iO) = reshape(SimResult.values{iO}',nInd*nT,1);
end

% Write header
fid = fopen(resultfile,'w');
fprintf(fid,'%s',header);
fprintf(fid,'\r\n');
% write reults
dlmwrite(resultfile,R,'-append','delimiter',';')
fclose(fid);

return

function SimResult = doIndividualsimulation(Settings,name,xml,parPathes,parValues,outputList,nInd,individualIdVector)

% initialize model
initStruct=[];
for iPar=1:length(parPathes)
    initStruct=initParameter(initStruct,[ '*|' parPathes{iPar}],'always','throwWarningIfNotExisting',false);
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
rowIndex_S = nan(size(outputList));

% initialize result array
SimResult.name = name;
SimResult.time = simTimeNorm;
SimResult.values = cell(length(outputList),1);
for iO = 1:length(outputList)
    SimResult.values{iO} = nan(nInd,length(simTimeNorm));
end
SimResult.individualIdVector = individualIdVector;
SimResult.outputList = outputList;
SimResult.outputUnit = cell(size(outputList));



% Loop On Individuals
for iInd=1:nInd
    
    % set physiology parameters
    % first time get row index
    if isempty(rowIndex)
        
        rowIndex=nan(size(parPathes));
        for iPar=1:length(parPathes)
            if existsParameter(['*|' parPathes{iPar}],simulationIndex);
                rowIndex(iPar)=setParameter(parValues(iInd,iPar),[ '*|' parPathes{iPar}],simulationIndex);
            else
                % check if is a demographic parameter which is not an xml
                % parameter, otherwise,  set a warning
                if ~ismember(parPathes{iPar},{'IndividualId','Gender','RaceIndex','Population Name','Organism|Weight','Organism|BMI'})                
                    writeToLog(sprintf('WARNING %s is parameter of the population csv, but not available as parameter in the xml', parPathes{iPar}),Settings.logfile,true,false);
                end
            end                
        end
    else
        % with row index use fast methode
        for iPar=find(~isnan(rowIndex))
            setParameter(parValues(iInd,iPar),['*|' parPathes{iPar}],simulationIndex,...
                'speedy','variable',rowIndex(iPar));
        end
    end
    
    
    success=processSimulation;
    
    if success
        % get result
        for iO=1:length(outputList)
            [simTime,tmpValues,rowIndex_S(iO)] = getSimulationResult(['*|' outputList{iO}],simulationIndex,'rowindex',rowIndex_S(iO));
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
for iO=1:length(outputList)
    if existsObserver(['*|' outputList{iO}],1)
        SimResult.outputUnit{iO} = getObserverFormula(['*|' outputList{iO}],1,'property','Unit');
    else    
        SimResult.outputUnit{iO} = getSpeciesInitialValue(['*|' outputList{iO}],1,'parametertype','readonly','property','Unit');
    end
    
end


return

