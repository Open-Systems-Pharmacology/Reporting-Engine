function runPopulationWorkflow(WSettings,TaskList,PopRunSet,VPC,Datafiles)
% master routine for population workflow
%
% runPopulationWorkflow(TaskList,WSettings,PopRunSet,varargin)
% 
% Inputs:
%       WSettings (structure)    definition of properties used in all
%                   workflow functions see GETDEFAULTWORKFLOWSETTINGS
%       TaskList (structure)    list of task which should be executed see GENERATEWORKFLOWINPUTFORPOPULATIONSIMULATION
%       PopRunSet (structure)   list of population simulations see GENERATEWORKFLOWINPUTFORPOPULATIONSIMULATION
%       VPC  (structure) contains information which plots should be
%               generated  see GETDEFAULTVPCPOPULATIONSETTINGS

% Open Systems Pharmacology Suite;  http://open-systems-pharmacology.org


% try
    %% check optional inputs
    if ~exist('Datafiles','var')
        Datafiles = {};
    end

    %% initialize workflow
    [WSettings] = initializeWorkflow('Population Simulation',WSettings);
    
    %% check if necessary inputs for task are available
    writeToLog(sprintf('Start input checks'),WSettings.logfile,true,false);
    
    successInputCheck = true;

    % simulation
    if TaskList.simulatePopulation
        for iSet = 1:length(PopRunSet)
            successInputCheck = checkPopulationSimulationInput(WSettings,PopRunSet(iSet)) & successInputCheck;
        end
    end

    % PK Parameter
    if TaskList.calculatePKParameter
        if ~ TaskList.simulatePopulation
            for iSet = 1:length(PopRunSet)
                if ~exist(fullfile('simulations',[PopRunSet(iSet).name '-Results.csv']),'file')
                    writeToLog(sprintf('ERROR: results for "%s" does not exist, please set Task simulatePopulation to true',PopRunSet(iSet).name),WSettings.logfile,true,false);
                    successInputCheck = false;
                end
                successInputCheck = checkPopulationSimulationInput(WSettings,PopRunSet(iSet)) & successInputCheck;
            end
        end
    end

    % data TP file
    % check if data read in is necessary
    doReadDataTPfile = TaskList.doVPC & ~isempty(Datafiles);
    
    % check file name of data read in
    if doReadDataTPfile
        if ~iscell(Datafiles) || size(Datafiles,2)~=3
            writeToLog(sprintf(['ERROR: variable Datafiles must be a cellarray size n x 3, \n',...
                'first entry datafile, second entry dictionary, third entry type of data']),WSettings.logfile,true,false);
            successInputCheck = false;
        else          
            for iData = 1:size(Datafiles,1)
                if ~exist(Datafiles{iData,1},'file')
                    writeToLog(sprintf('ERROR: Datafile %s does not exist',Datafiles{iData,1}),WSettings.logfile,true,false);
                    successInputCheck = false;
                end
                if ~exist(Datafiles{iData,2},'file')
                    writeToLog(sprintf('ERROR: Dictionary %s for timeprofile datafile does not exist',Datafiles{iData,2}),WSettings.logfile,true,false);
                    successInputCheck = false;
                end
                if ~ismember(Datafiles{iData,3},{'timeprofile'})
                    writeToLog(sprintf('ERROR: Datatype %s is unknown.',Datafiles{iData,3}),WSettings.logfile,true,false);
                    successInputCheck = false;
                end
            end
        end
    end

    % PK Parameter
    if TaskList.doVPC
        if isempty(VPC)
            writeToLog(sprintf('ERROR: VPC definition is missing!'),WSettings.logfile,true,false);
            successInputCheck = false;
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
        [success,VPC] = preparePopulationSimulation(WSettings,PopRunSet(iSet),VPC);
        successInputCheck = success & successInputCheck;
    end
    % ontogeny factors are added now, get rid of keay word
    if ~isempty(VPC)
        jj = strcmp(VPC.PhysProperties.yList(:,1),'<addOntogenyFactor>');
        if any(jj)
            VPC.PhysProperties.yList = VPC.PhysProperties.yList(~jj,:);
        end
    end
        

    
    % stop if inpt is inconsistent
    if ~successInputCheck
        return
    end
    
    % read datafiles
    if doReadDataTPfile
        for iData = 1:size(Datafiles,1)
            switch Datafiles{iData,3}
                case 'timeprofile'
                    successInputCheck = readTimeprofileDataForPopulation(WSettings,Datafiles(iData,[1 2]),PopRunSet);
                otherwise
                    error ('unknown datatype')
            end
        end
    end              
    
    % stop if inpt is inconsistent
    if ~successInputCheck
        return
    end

    
    writeToLog(sprintf('Input checks were successfully executed. \n'),WSettings.logfile,true,false);

    %% Tasks processing
    
    % simulate the poulationssimulations
    if TaskList.simulatePopulation
        
        for iSet = 1:length(PopRunSet)
            runPopulationSimulation(WSettings,PopRunSet(iSet));
        end
    end

    % calculate PK Parameter
    if TaskList.calculatePKParameter
        for iSet = 1:length(PopRunSet)
            calculatePopulationPKParameter(WSettings,PopRunSet(iSet));
        end
    end

    % generate VPC
    if TaskList.doVPC
        runPopulationVPC(WSettings,VPC,PopRunSet);
    end

    
% catch exception
%     
%     save('exception.mat','exception')
%     writeToLog(exception.message,WSettings.logfile,true,false);
%     
% end



return


function [successInputCheck,VPC] = preparePopulationSimulation(WSettings,PopRunSet,VPC)

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
    [parPaths,parValues,ixStudyDesign] = readStudyDesign(WSettings,PopRunSet.studyDesign,parPaths,parValues);
    
    parPathsStudyDesign = parPaths(ixStudyDesign);
    parValuesStudyDesign = parValues(:,ixStudyDesign);
else
    ixStudyDesign = []; %#ok<NASGU>
    parPathsStudyDesign = [];
    parValuesStudyDesign = [];
end

save(fullfile(tmpDir,'pop.mat'),'parPaths','parValues','ixStudyDesign');


% analyse applicationProtocol
[ApplicationProtocol,isValid] = getApplicationProtocollFromXML(WSettings,PopRunSet.xml,parPathsStudyDesign,parValuesStudyDesign);   %#ok<NASGU>
save(fullfile(tmpDir,'applicationProtocol.mat'),'ApplicationProtocol','isValid');

simulationIndex = 1; %Initialised in getApplicationProtocollFromXML

if isempty(PopRunSet.calculatePKParameterFh)
    PKParameterTemplate = calculatePKParameterForApplicationProtocol(WSettings,ApplicationProtocol);
else
    PKParameterTemplate = feval(PopRunSet.calculatePKParameterFh,WSettings,ApplicationProtocol);
end

save(fullfile(tmpDir,'applicationProtocol.mat'),'-append','PKParameterTemplate');


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
                writeToLog(sprintf('WARNING: For "%s" no base unit could be identified, may be cause problems if plotted',parPaths{iPar}),WSettings.logfile,true,false);
                
        end
    end
end
save(fullfile(tmpDir,'pop.mat'),'-append','unit');


% analyse outputs, check if existing and add unitfactor
OutputList = readOutputCsv(PopRunSet.outputcsv);
for iO = 1:length(OutputList)
    
    [unitFactor,success,MW] = getUnitfactorForOutputPath(WSettings,OutputList(iO).pathID,OutputList(iO).displayUnit,simulationIndex);
    successInputCheck = success & successInputCheck;
    OutputList(iO).unitFactor = unitFactor; 
    
    % check PK Parameter
    for iPK = 1:size(OutputList(iO).pKParameterList,2)
        
        jj = strcmp(OutputList(iO).pKParameterList{1,iPK},{PKParameterTemplate.name});
        if ~any(jj)
            writeToLog(sprintf('ERROR: Output %s has listed PK Parameter "%s" which is not calculated by the PK parameter function.',...
                OutputList(iO).pathID,OutputList(iO).pKParameterList{1,iPK}),WSettings.logfile,true,false);
            successInputCheck = false;
        else
            
            [unitFactor,success] = getUnitFactorForUnknownDimension(WSettings,PKParameterTemplate(jj).unit,OutputList(iO).pKParameterList{2,iPK},MW);
            
            successInputCheck = success & successInputCheck;
           OutputList(iO).pKParameterList{3,iPK} = unitFactor;
        end
    end
        
    
end

save(fullfile(tmpDir,'outputList.mat'),'OutputList');

% check if ontogenyFactors should be added to the physology plotList
if ~isempty(VPC) && ...
 any(strcmp(VPC.PhysProperties.yList(:,1),'<addOntogenyFactor>'))
    ontogenyIDs = getParameter('*|Ontogeny factor*',1,'parametertype','readonly','property','ID');
    jj = getParameter(ontogenyIDs,1,'parametertype','readonly','property','isFormula');
    
    ontogenyPaths = getParameter(ontogenyIDs(~jj),1,'parametertype','readonly','property','Path');
    
    for iO = 1:length(ontogenyPaths)
        tmp = regexp(ontogenyPaths{iO},'\|','split');
        pathID = strjoin(tmp(2:end),'|');
        if ~ismember(pathID,VPC.PhysProperties.yList(:,1))
            if ismember(tmp{end},{'Ontogeny factor','Ontogeny factor GI'})
                ontogenyName = strjoin(tmp([end 2:end-1]),' ');
            else
                ontogenyName = tmp{end};
            end
            VPC.PhysProperties.yList(end+1,:) = {pathID,ontogenyName,'',{}};
        end
    end
end


return
    

function successInputCheck = checkPopulationSimulationInput(WSettings,PopRunSet)

successInputCheck = true;

% check name on special letters
if any(ismember(PopRunSet.name,'./ ?§$%&()[]{}+~*#'))
    writeToLog(sprintf('ERROR: Popsetname "%s" contains special signs, do not use them.',PopRunSet.name),WSettings.logfile,true,false);
    successInputCheck = false;
end

% mandatory inputs
if ~exist(PopRunSet.xml,'file')
    writeToLog(sprintf('ERROR: "%s" does not exist',PopRunSet.xml),WSettings.logfile,true,false);
    successInputCheck = false;
end
if ~exist(PopRunSet.popcsv,'file')
    writeToLog(sprintf('ERROR: "%s" does not exist',PopRunSet.popcsv),WSettings.logfile,true,false);
    successInputCheck = false;
end
if ~exist(PopRunSet.outputcsv,'file')
    writeToLog(sprintf('ERROR: "%s" does not exist',PopRunSet.outputcsv),WSettings.logfile,true,false);
    successInputCheck = false;
end
                
% optional inputs
if  ~isempty(PopRunSet.studyDesign) && ~exist(PopRunSet.studyDesign,'file') 
    writeToLog(sprintf('ERROR: "%s" is given, but does not exist',PopRunSet.studyDesign),WSettings.logfile,true,false);
    successInputCheck = false;
end
 
return


function successInputCheck = readTimeprofileDataForPopulation(WSettings,dataTpFile,PopRunSet) 

successInputCheck = true; %#ok<NASGU>

% collect filters
k=0;
filterList = {};
for iSet = 1:length(PopRunSet)
    
    if ~isempty(PopRunSet(iSet).dataTpFilter)
        k=k+1;
        filterList{k} = PopRunSet(iSet).dataTpFilter; %#ok<AGROW>
        filterIndex(iSet,1) = k; %#ok<AGROW>
    else
        filterIndex(iSet,1) = nan; %#ok<AGROW>
    end
    
    load(fullfile('tmp',PopRunSet(iSet).name,'outputList.mat'),'OutputList');
    
    for iO = 1:length(OutputList);
        
        if ~isempty(OutputList(iO).dataTpFilter)
            k=k+1;
            filterList{k,:} = OutputList(iO).dataTpFilter; %#ok<AGROW>
            filterIndex(iSet,iO+1) = k; %#ok<AGROW>
        else
            filterIndex(iSet,iO+1) = nan;%#ok<AGROW>
        end

    end
end
    
    
% read nonmefile uses WSettings and dataTpFile
[successInputCheck,X,filter,Dict] = readNonmemFile(WSettings,dataTpFile,'timeprofile',filterList);

if ~successInputCheck
    return
end



% get fieldnamees of cavariates
jj = strcmp({Dict.type},'covariate');
fn_covariate = intersect(fieldnames(X),{Dict(jj).matlabID});


% get dataset for each popset
 
for iSet = 1:length(PopRunSet)
    
    % get temporary directory
    tmpDir = fullfile('tmp',PopRunSet(iSet).name);

    
    % initialize TP
    for iO = 1:length(OutputList)
        TP{iO} = []; %#ok<AGROW>
    end

    
    % get filter of popset
    if ~isnan(filterIndex(iSet,1))
        jj_set = filter(:,filterIndex(iSet,1));
    else
        jj_set = false(size(X.stud));
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
            
            
            for iO = 1:length(OutputList)

                % get filter of popset
                if ~isnan(filterIndex(iSet,1+iO))
                    jj_O = jj_SID & filter(:,filterIndex(iSet,1+iO));
                else
                    jj_O = false(size(X.stud));
                end
                
                if any(jj_O)
                    
                    % add to data structure
                    DataTP(indx).stud = uniSTUD(iStud); %#ok<AGROW>
                    DataTP(indx).sid = uniSID(iSID) ; %#ok<AGROW>
                    
                    % timeprofile
                    TP{iO}(indx).time = X.time(jj_O); %#ok<AGROW>
                    TP{iO}(indx).dv = X.dv(jj_O); %#ok<AGROW>
                    if isfield(X,'tad')
                        TP{iO}(indx).tad = X.tad(jj_O); %#ok<AGROW>
                    end
                    TP{iO}(indx).isLloq = false(size( TP{iO}(indx).time)); %#ok<AGROW>
                    TP{iO}(indx).lloq = nan; %#ok<AGROW>
                    if isfield(X,'lloq')
                        lloq =X.lloq(jj_O);
                        jj_lloq = TP{iO}(indx).dv < lloq;
                        if any(jj_lloq)
                            TP{iO}(indx).isLloq(jj_lloq) = true; %#ok<AGROW>
                            tmp =  unique(lloq(jj_lloq));
                            if length(tmp) > 1
                                writeToLog(sprintf('WARNING: non unique lloq for STUD %d SID %d', ...
                                    DataTP(indx).stud ,DataTP(indx).sid),WSettings.logfile,true,false);
                            end
                            TP{iO}(indx).lloq = min(tmp); %#ok<AGROW>
                        end

                    end
                
                    % covariates
                    for iFn = 1:length(fn_covariate)
                        DataTP(indx).(fn_covariate{iFn}) = X.(fn_covariate{iFn})(find(jj_SID,1)); %#ok<AGROW>
                
                    end
                end
            end
        end
    end
    
    % save data
    if exist('DataTP','var')
        save(fullfile(tmpDir,'dataTp.mat'),'DataTP','TP','Dict');
    end
end

return
