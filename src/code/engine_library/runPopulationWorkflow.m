function runPopulationWorkflow(WSettings,TaskList,PopRunSet,workflowType,VPC,Datafiles,sensParameterList)
% master routine for population workflow
%
% runPopulationWorkflow(TaskList,WSettings,PopRunSet,varargin)
% 
% Inputs:
%       WSettings (structure)    definition of properties used in all
%                   workflow functions see GETDEFAULTWORKFLOWSETTINGS
%       TaskList (structure)    list of task which should be executed see GENERATEWORKFLOWINPUTFORPOPULATIONSIMULATION
%       PopRunSet (structure)   list of population simulations see GENERATEWORKFLOWINPUTFORPOPULATIONSIMULATION
%       workflowType (string)  defines type of population, default = 'paralellComparison'
%                               other possibilities are 'pediatric', ratioComparison;
%       VPC  (structure) contains information which plots should be
%               generated  see GETDEFAULTVPCPOPULATIONSETTINGS
%       sensParameterList  (cellarry) 1. column pathid of parameter,
%                       2. number of steps
%                       3. variation range
%                       4. minimal value
%                       5. maximal value

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
                tmp = dir(fullfile('simulations',[PopRunSet(iSet).name '*-Results.csv']));
                if isempty(tmp)
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
        successInputCheck = checkInputDatafiles(WSettings,Datafiles) & successInputCheck;
    end

    % PK Parameter
    if TaskList.doVPC
        if isempty(VPC)
            writeToLog(sprintf('ERROR: VPC definition is missing!'),WSettings.logfile,true,false);
            successInputCheck = false;
        end
    end
    
    % PK Parameter
    if TaskList.doSensitivityAnalysis
        if isempty(sensParameterList)
            writeToLog(sprintf('ERROR: sensitivity definition is missing!'),WSettings.logfile,true,false);
            successInputCheck = false;
        end
        % for sensitivit population workflow same output is needed
        success = checkOutputConsistencyForSensitivity(WSettings,PopRunSet,workflowType);
        successInputCheck = successInputCheck & success;
        
        % PK Parameter are needed
        if ~ TaskList.calculatePKParameter
            for iSet = 1:length(PopRunSet)
                tmp = dir(fullfile('simulations',[PopRunSet(iSet).name '*PK-Analyses.csv']));
                if isempty(tmp)
                    writeToLog(sprintf('ERROR: PK Parameter results for "%s" does not exist, please set Task calculatePKParameter to true',...
                        PopRunSet(iSet).name),WSettings.logfile,true,false);
                    successInputCheck = false;
                end
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
        [success,VPC] = preparePopulationSimulation(WSettings,PopRunSet(iSet),VPC,sensParameterList);
        successInputCheck = success & successInputCheck;
    end
    % ontogeny factors are added now, get rid of key word
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
                    successInputCheck = readTimeprofileDataForSimulation(WSettings,Datafiles(iData,[1 2]),PopRunSet);
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

    % generate sesnitivity analysis
    if TaskList.doSensitivityAnalysis
        runPopulationSensitivity(WSettings,PopRunSet,workflowType);
    end

    
% catch exception
%     
%     save('exception.mat','exception')
%     writeToLog(exception.message,WSettings.logfile,true,false);
%     
% end
% 


return


function [successInputCheck,VPC] = preparePopulationSimulation(WSettings,PopRunSet,VPC,sensParameterList)

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


% add BSA if not available
if ~any(strcmp(parPaths,'Organism|BSA'))
    parPaths{end+1} = 'Organism|BSA';
    
    jj = strcmp('Organism|Weight',parPaths);
    weight =  parValues(:,jj);
    jj = strcmp('Organism|Height',parPaths);
    height =  parValues(:,jj);
    parValues(:,end+1) = calculateBodySurfaceArea(WSettings,weight,height);
end


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

% initialize PK Parameter
if isempty(PopRunSet.calculatePKParameterFh)
    PKParameterTemplate = calculatePKParameterForApplicationProtocol(WSettings,ApplicationProtocol);
else
    PKParameterTemplate = feval(PopRunSet.calculatePKParameterFh,WSettings,ApplicationProtocol);
end

save(fullfile(tmpDir,'applicationProtocol.mat'),'-append','PKParameterTemplate');


% get defaultUnits for population parameter
unit = cell(length(parPaths),1);
for iPar = 1:length(parPaths)
    [ise,desc] = existsParameter(['*' parPaths{iPar}],simulationIndex,'parametertype','readonly');
    if ise
        unit{iPar} = desc{2,strcmp(desc(1,:),'Unit')};
    else
        switch parPaths{iPar}
            case {'IndividualId','Gender','RaceIndex','Population Name'}
                unit{iPar} = 'none';
            case 'Organism|BSA'
                 unit{iPar} = 'dm²';
            otherwise
                unit{iPar} = 'none';
                writeToLog(sprintf('WARNING: For "%s" no base unit could be identified, may be cause problems if plotted',parPaths{iPar}),WSettings.logfile,true,false);
                
        end
    end
end
save(fullfile(tmpDir,'pop.mat'),'-append','unit');


% analyse outputs, check if existing and add unitfactor
OutputList = PopRunSet.OutputList;
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

% check sensitivity definition
if ~isempty(sensParameterList)
    sensParameterList = checkSensitivityParameter(WSettings,sensParameterList,simulationIndex);     %#ok<NASGU>
    save(fullfile(tmpDir,'senssitivity.mat'),'sensParameterList');
end


% check if ontogenyFactors should be added to the physiology plotList
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
                
% optional inputs
if  ~isempty(PopRunSet.studyDesign) && ~exist(PopRunSet.studyDesign,'file') 
    writeToLog(sprintf('ERROR: "%s" is given, but does not exist',PopRunSet.studyDesign),WSettings.logfile,true,false);
    successInputCheck = false;
end
 
return


function  successInputCheck = checkOutputConsistencyForSensitivity(WSettings,PopRunSet,workflowType)

successInputCheck = true;

if strcmp(workflowType,'pediatric')
    jjRef = [PopRunSet.isReference];
    setList = find(~jjRef);
else
    setList = 1:length(PopRunSet);
end

%  get outputs and PKPList
baseOutputList = PopRunSet(setList(1)).OutputList;

for iSet = setList(2:end)

    OutputList = PopRunSet(iSet).OutputList;
    
    if ~(length(OutputList) == length(baseOutputList)) || ...
        ~all(strcmp({OutputList.pathID},{baseOutputList.pathID})) 
            writeToLog(sprintf('ERROR: Outputs must be the same for a sensitivity analysis!'),WSettings.logfile,true,false);
            successInputCheck = false;
    end
    
    if successInputCheck
        for iO = 1:length(baseOutputList)
            if ~(length(OutputList(iO).pKParameterList(1,:)) == length(baseOutputList(iO).pKParameterList(1,:))) || ...
                    ~all(strcmp(OutputList(iO).pKParameterList(1,:),baseOutputList(iO).pKParameterList(1,:)))
                writeToLog(sprintf('ERROR: PK parameter of Outputs must be the same for a sensitivity analysis!'),WSettings.logfile,true,false);
                successInputCheck = false;                
            end
        end
    end
    
end
