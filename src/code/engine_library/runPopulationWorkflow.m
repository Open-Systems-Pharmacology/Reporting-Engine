function runPopulationWorkflow(TaskList,Settings,PopRunSet,varargin)
% master routine for population workflow
%
% Inputs:
%       TaskList (structure)    list of task which should be executed see GETDEFAULTTASKLISTPOPULATIONWORKFLOW
%       Settings (structure)    definition of properties used in all
%                   workflow functions see GETDEFAULTWORKFLOWSETTINGS
%       PopRunSet (structure)   list of population simulations see GETDEFAULTPOPRUNSET

% Open Systems Pharmacology Suite;  http://forum.open-systems-pharmacology.org
% Date: 14-July-2017


% try
    %% check optional inputs
    [dataTpFile,VPC] = ...
    checkInputOptions(varargin,{...
    'dataTpFile','',{},...
    'VPC','',[],...
    });

    %% initialize workflow
    [Settings] = initializeWorkflow('Population Simulation',Settings);
    
    %% check if necessary inputs for task are available
    writeToLog(sprintf('Start input checks'),Settings.logfile,true,false);
    
    successInputCheck = true;

    % simulation
    if TaskList.simulatePopulation
        for iSet = 1:length(PopRunSet)
            successInputCheck = checkPopulationSimulationInput(Settings,PopRunSet(iSet)) & successInputCheck;
        end
    end

    % PK Parameter
    if TaskList.calculatePKParameter
        if ~ TaskList.simulatePopulation
            for iSet = 1:length(PopRunSet)
                if ~exist(fullfile('Simulations',[PopRunSet(iSet).name '-Results.csv']),'file')
                    writeToLog(sprintf('ERROR: results for "%s" does not exist, please set Task simulatePopulation to true',PopRunSet(iSet).name),Settings.logfile,true,false);
                    successInputCheck = false;
                end
                successInputCheck = checkPopulationSimulationInput(Settings,PopRunSet(iSet)) & successInputCheck;
            end
        end
    end

    % data TP file
    % check if data read in is necessary
    doReadDataTPfile = TaskList.doVPC & ~isempty(dataTpFile);
    
    % chekc file name of data read in
    if doReadDataTPfile
        if ~iscell(dataTpFile) || length(dataTpFile)~=2
            writeToLog(sprintf('ERROR: variable dataTpFile for timeprofile data must be a cellarray, first entry datafile, second entry dictionary'),Settings.logfile,true,false);
            successInputCheck = false;
        else           
            if ~exist(dataTpFile{1},'file')
                writeToLog(sprintf('ERROR: Datafile %s does not exist',dataTpFile{1}),Settings.logfile,true,false);
                successInputCheck = false;
            end
            if ~exist(dataTpFile{2},'file')
                writeToLog(sprintf('ERROR: Dictionary %s for timeprofile datafile does not exist',dataTpFile{2}),Settings.logfile,true,false);
                successInputCheck = false;
            end
        end
    end

    
    % stop if not all inputs are available
    if ~successInputCheck
        return
    end

    
    % evaluate xmls, check for consistency  and create temporary
    % directory with infos like application protocols, and factor for
    % unitconversion 
    for iSet = 1:length(PopRunSet)
        successInputCheck = preparePopulationSimulationConsistency(Settings,PopRunSet(iSet)) & successInputCheck;
    end
    
    % stop if inpt is inconsistent
    if ~successInputCheck
        return
    end
    
    % read datafiles
    if doReadDataTPfile
        successInputCheck = readTimeprofileDataForPopulation(Settings,dataTpFile,PopRunSet);
    end
    
    % stop if inpt is inconsistent
    if ~successInputCheck
        return
    end

    
    writeToLog(sprintf('Input checks were successfully executed. \n'),Settings.logfile,true,false);

    %% Tasks processing
    
    % simulate the poulationssimulations
    if TaskList.simulatePopulation
        
        for iSet = 1:length(PopRunSet)
            runPopulationSimulation(Settings,PopRunSet(iSet));
        end
    end

    % calculate PK Parameter
    if TaskList.calculatePKParameter
        for iSet = 1:length(PopRunSet)
            calculatePopulationPKParameter(Settings,PopRunSet(iSet));
        end
    end

    % generate VPC
    if TaskList.doVPC
        runPopulationVPC(Settings,VPC,PopRunSet);
    end

    
% catch exception
%     
%     save('exception.mat','exception')
%     writeToLog(exception.message,Settings.logfile,true,false);
%     
% end



return


function successInputCheck = preparePopulationSimulationConsistency(Settings,PopRunSet)

% initialize return value
successInputCheck = true;

% create temporary diretory
tmpDir = fullfile('tmp',PopRunSet.name); 
% check if siluation directory already exist
if ~exist(tmpDir,'dir')
    mkdir(tmpDir)
else
    % clear up all temporary results
    delete(fullfile(tmpDir,'*.mat'));
end


% read population csv
[parPaths,parValues] = readPopulationCSV(PopRunSet.popcsv);

% add studyDesign if one is given
if ~isempty(PopRunSet.studyDesign)
    [parPaths,parValues,ixStudyDesign] = readStudyDesign(Settings,PopRunSet.studyDesign,parPaths,parValues);
    
    parPathsStudyDesign = parPaths(ixStudyDesign);
    parValuesStudyDesign = parValues(:,ixStudyDesign);
else
    ixStudyDesign = []; %#ok<NASGU>
    parPathsStudyDesign = [];
    parValuesStudyDesign = [];
end

save(fullfile(tmpDir,'pop.mat'),'parPaths','parValues','ixStudyDesign');


% analyse applicationProtocol
[ApplicationProtocol,isValid] = getApplicationProtocollFromXML(Settings,PopRunSet.xml,parPathsStudyDesign,parValuesStudyDesign);   %#ok<NASGU,ASGLU>
save(fullfile(tmpDir,'applicationProtocol.mat'),'ApplicationProtocol','isValid');

simulationIndex = 1; %Initialised in getApplicationProtocollFromXML


% get defaultUnits for population
unit = cell(length(parPaths),1);
for iPar = 1:length(parPaths)
    [ise,desc] = existsParameter(['*' parPaths{iPar}],simulationIndex,'parametertype','readonly');
    if ise
        unit{iPar} = desc{2,strcmp(desc(1,:),'Unit')};
    else
        switch parPaths{iPar}
            case {'IndividualId','Gender','RaceIndex','Population Name'}
                unit{iPar} = 'none';
            otherwise
                unit{iPar} = 'none';
                writeToLog(sprintf('WARNING: For "%s" no base unit could be identified, may be cause problems if plotted',parPaths{iPar}),Settings.logfile,true,false);
                
        end
    end
end
save(fullfile(tmpDir,'pop.mat'),'-append','unit');


% analyse outputs, chekc if existing and add unitfactor
OutputList = PopRunSet.OutputList;
for iO = 1:length(OutputList)
    
    [unitFactor,success] = getUnitfactorForOutputPath(Settings,OutputList(iO).pathID,OutputList(iO).displayUnit,simulationIndex);
    successInputCheck = success & successInputCheck;
    OutputList(iO).unitFactor = unitFactor; 
end

save(fullfile(tmpDir,'outputList.mat'),'OutputList');

return
    

function successInputCheck = checkPopulationSimulationInput(Settings,PopRunSet)

successInputCheck = true;

% check name on special letters
if any(ismember(PopRunSet.name,'./ ?§$%&()[]{}+~*#'))
    writeToLog(sprintf('ERROR: Popsetname "%s" contains special signs, do not use them.',PopRunSet.name),Settings.logfile,true,false);
    successInputCheck = false;
end

% mandatory inputs
if ~exist(PopRunSet.xml,'file')
    writeToLog(sprintf('ERROR: "%s" does not exist',PopRunSet.xml),Settings.logfile,true,false);
    successInputCheck = false;
end
if ~exist(PopRunSet.popcsv,'file')
    writeToLog(sprintf('ERROR: "%s" does not exist',PopRunSet.popcsv),Settings.logfile,true,false);
    successInputCheck = false;
end
if isempty(PopRunSet.OutputList)
    writeToLog(sprintf('ERROR: outputList is empty'),Settings.logfile,true,false);
    successInputCheck = false;
end
                
% optional inputs
if  ~isempty(PopRunSet.studyDesign) && ~exist(PopRunSet.studyDesign,'file') 
    writeToLog(sprintf('ERROR: "%s" is given, but does not exist',PopRunSet.studyDesign),Settings.logfile,true,false);
    successInputCheck = false;
end
 
return


function successInputCheck = readTimeprofileDataForPopulation(Settings,dataTpFile,PopRunSet) 

successInputCheck = true; %#ok<NASGU>

% collect filters
k=0;
for iSet = 1:length(PopRunSet)
    
    if ~isempty(PopRunSet(iSet).dataTpFilter)
        k=k+1;
        filterList{k} = PopRunSet(iSet).dataTpFilter; %#ok<AGROW>
        PopRunSet(iSet).filterIndex = k;
    else
        PopRunSet(iSet).filterIndex = nan;
    end
    

    for iO = 1:length(PopRunSet(iSet).OutputList);
        
        if ~isempty(PopRunSet(iSet).OutputList(iO).dataTpFilter)
            k=k+1;
            filterList{k,:} = PopRunSet(iSet).OutputList(iO).dataTpFilter; %#ok<AGROW>
            PopRunSet(iSet).OutputList(iO).filterIndex = k;
        else
            PopRunSet(iSet).OutputList(iO).filterIndex = nan;
        end

    end
end
    
    
% read nonmefile uses Settings and dataTpFile
[successInputCheck,X,filter,dict] = readNonmemFile(Settings,dataTpFile,'timeprofile',filterList);

if ~successInputCheck
    return
end



% get fieldnamees of cavariates
jj = strcmp({dict.type},'covariate');
fn_covariate = intersect(fieldnames(X),{dict(jj).matlabID});


% get dataset for each popset
for iSet = 1:length(PopRunSet)

    % get filter of popset
    if ~isnan(PopRunSet(iSet).filterIndex)
        jj_set = filter(:,PopRunSet(iSet).filterIndex);
    else
        jj_set = true(size(X.stud));
    end
    
    % get unique identifier
    uniSTUD = unique(X.stud(jj_set));
        

    for iStud = 1:length(uniSTUD)
        
        jj_STUD = jj_set & X.stud == uniSTUD(iStud);
        
        uniSID = unique(X.sid(jj_STUD));
        
        for iSID = 1:length(uniSID)
            
            jj_SID = jj_STUD & X.sid == uniSID(iSID) ;
            
            % get index of dataset
            if exist('DataTP','var')
                indx = length(DataTP)+1;
            else
                indx = 1;
            end
            
            for iO = 1:length(PopRunSet(iSet).OutputList);

                % get filter of popset
                if ~isnan(PopRunSet(iSet).OutputList(iO).filterIndex)
                    jj_O = jj_SID & filter(:,PopRunSet(iSet).OutputList(iO).filterIndex);
                else
                    jj_O = jj_SID;
                end
                
                if any(jj_O)
                    
                    % add to data structure
                    DataTP(indx).stud = uniSTUD(iStud); %#ok<AGROW>
                    DataTP(indx).sid = uniSID(iSID) ; %#ok<AGROW>
                    DataTP(indx).Tp(iO).time = X.time(jj_O); %#ok<AGROW>
                    DataTP(indx).Tp(iO).dv = X.dv(jj_O); %#ok<AGROW>
                    if isfield(X,'tad')
                        DataTP(indx).Tp(iO).tad = X.tad(jj_O); %#ok<AGROW>
                    end
                    DataTP(indx).Tp(iO).lloq = nan(size( DataTP(indx).Tp(iO).time)); %#ok<AGROW>
                    if isfield(X,'lloq')
                        lloq =X.lloq(jj_O);
                        jj_lloq = DataTP(indx).Tp(iO).dv < lloq;
                        DataTP(indx).Tp(iO).lloq(jj_lloq) = lloq(jj_lloq); %#ok<AGROW>
                    end
                
                    % covariates
                    for iFn = 1:length(fn_covariate)
                        DataTP(indx).(fn_covariate{iFn}) = X.(fn_covariate{iFn})(find(jj_SID,1)); %#ok<AGROW>
                
                    end
                end
            end
        end
    end
    
    % get temporary directory
    tmpDir = fullfile('tmp',PopRunSet(iSet).name);
    % save data
    save(fullfile(tmpDir,'dataTp.mat'),'DataTP','dict');
end

return
