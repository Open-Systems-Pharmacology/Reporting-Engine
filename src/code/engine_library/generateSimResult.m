function SimResult = generateSimResult(WSettings,SimulationSet,simulationIndex,parPaths,parValues,outputPathList,nInd,individualIdVector)
% GENERATESIMRESULT processes given  simulation with differnet parameter values
%
% Inputs
%   WSettings  structure containing global settings see GETDEFAULTWORKFLOWSETTINGS
%   name (string) naem of simulation
%   SimulationSet (structure) structure with name of xmlfile and simulation time may be empty, if model is already
%       initialized, but then simulationIndex is needed.
%   simulationIndex (nan) if xmls is simulationIndex is needed
%   parPaths (cellarray of strings) pathnames of parameter
%   parValues ( double array) values of parameter
%   outputPathList (list of desired output
%   nInd (double) number of individuals
%   individualIdVector  (double vector) vector of individuals ids


% Open Systems Pharmacology Suite;  http://open-systems-pharmacology.org

% initialize model
if ~isempty(SimulationSet)
        
    initStruct=[];
    for iPar=1:length(parPaths)
        initStruct=initParameter(initStruct,[ '*|' parPaths{iPar}],'always','throwWarningIfNotExisting',false);
    end
    initSimulation(SimulationSet.xml,initStruct);
    
    % Only one simulation is intialized so the simulation index is always 1
    simulationIndex = 1;
    
    % set simulaion time if given
    if isfield(SimulationSet,'simulationTime') && ~isempty(PopRunSet.simulationTime)
        setSimulationTime(SimulationSet.simulationTime,simulationIndex)
    end

end


% getSimulationTime For the import all results must have the same
% timepoints, which happens to be not the case for simulations outside PK-Sim 
simTimeNorm = getSimulationTime(simulationIndex);

% initialize rowindex of parameter for speedy setting
rowIndex = [];
% initialize rowindex of output for speedy fetching
rowIndexS = nan(size(outputPathList));

% initialize result array
SimResult.time = simTimeNorm;
SimResult.values = cell(length(outputPathList),1);
for iO = 1:length(outputPathList)
    SimResult.values{iO} = nan(length(simTimeNorm),nInd);
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
                        parPaths{iPar}),WSettings.logfile,true,false);
                    
                    for iTmp = 2:length(tmp)
                        indx = length(rowIndex)+1;
                        rowIndex(indx) = tmp(iTmp);
                        parValues(:,indx) = parValues(:,iPar);
                    end
                end
                
            else
                % check if is a demographic parameter which is not an xml
                % parameter, otherwise,  set a warning
                if ~ismember(parPaths{iPar},{'IndividualId','Gender','RaceIndex','Population Name',...
                        'Organism|Weight','Organism|BMI','Organism|BSA'})                
                    writeToLog(sprintf('WARNING %s is parameter of the population csv, but not available as parameter in the xml', ...
                        parPaths{iPar}),WSettings.logfile,true,false);
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
            [simTime,tmpValues,rowIndexS(iO)] = getSimulationResult(['*|' outputPathList{iO}],simulationIndex,'rowindex',rowIndexS(iO));
            % map outputs to defined time raster
            if ~all(simTime - simTimeNorm==0)
                SimResult.values{iO}(:,iInd) = interp1(simTime,tmpValues,simTimeNorm,'linear','extrap');
            else
                SimResult.values{iO}(:,iInd) = tmpValues;
            end
        end
    else
        writeToLog(sprintf('WARNING solver error, at individual %d',iInd),WSettings.logfile,true,false);
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

