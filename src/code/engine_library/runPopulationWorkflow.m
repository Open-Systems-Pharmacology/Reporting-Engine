function runPopulationWorkflow(WSettings,TaskList,PopRunSet,VPC,Datafiles,sensParameterList)
% master routine for population workflow
%
% runPopulationWorkflow(WSettings,TaskList,PopRunSet,VPC,Datafiles,sensParameterList)
% 
% Inputs:
%       WSettings (structure)    definition of properties used in all
%                   workflow functions see GETDEFAULTWORKFLOWSETTINGS
%       TaskList (structure)    list of task which should be executed see GENERATEWORKFLOWINPUTFORPOPULATIONSIMULATION
%       PopRunSet (structure)   list of population simulations see GENERATEWORKFLOWINPUTFORPOPULATIONSIMULATION
%       VPC  (structure) contains information which plots should be
%               generated  see GETDEFAULTVPCPOPULATIONSETTINGS
%       sensParameterList  (cellarry) 1. column pathid of parameter,
%                       2. number of steps
%                       3. variation range
%                       4. minimal value
%                       5. maximal value

% Open Systems Pharmacology Suite;  http://open-systems-pharmacology.org


try
    %% check optional inputs
    if ~exist('Datafiles','var')
        Datafiles = {};
    end

    %% initialize workflow
    
    [WSettings] = initializeWorkflow(WSettings);
    
    %% check if necessary inputs for task are available
    [successInputCheck,VPC] = inputCheck(WSettings,TaskList,PopRunSet,VPC,Datafiles,sensParameterList);
    if ~successInputCheck
        return
    end
   
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
        runPopulationSensitivity(WSettings,PopRunSet);
    end
    
    writeToReportLog('INFO',sprintf('Population workflow finalized \n'),false);

    
catch exception
    
    save(sprintf('exception_%s.mat',datestr(now,'ddmmyy_hhMM')),'exception');
    writeToReportLog('ERROR',exception.message,false);
    writeToReportLog('INFO',sprintf('Population workflow finished with error \n'),false);
  
end



return

function  [successInputCheck,VPC] = inputCheck(WSettings,TaskList,PopRunSet,VPC,Datafiles,sensParameterList)
    writeToReportLog('INFO',sprintf('Start input checks'),false);
    
    successInputCheck = true;

    % simulation
    if TaskList.simulatePopulation
        for iSet = 1:length(PopRunSet)
            successInputCheck = checkPopulationSimulationInput(PopRunSet(iSet)) & successInputCheck;
        end
    end

    % PK Parameter
    if TaskList.calculatePKParameter
        if ~ TaskList.simulatePopulation
            for iSet = 1:length(PopRunSet)
                tmp = dir(fullfile('simulations',[PopRunSet(iSet).name '*-Results.csv']));
                if isempty(tmp)
                    writeToReportLog('ERROR',sprintf('results for "%s" does not exist, please set Task simulatePopulation to true',...
                        PopRunSet(iSet).name),false);
                    successInputCheck = false;
                end
                successInputCheck = checkPopulationSimulationInput(PopRunSet(iSet)) & successInputCheck;
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
            writeToReportLog('ERROR',sprintf('VPC definition is missing!'),false);
            successInputCheck = false;
        end
    end
    
    % PK Parameter
    if TaskList.doSensitivityAnalysis
        if isempty(sensParameterList)
            writeToReportLog('ERROR','sensitivity definition is missing!',false);
            successInputCheck = false;
        end
        % for sensitivit population workflow same output is needed
        success = checkOutputConsistencyForSensitivity(WSettings,PopRunSet);
        successInputCheck = successInputCheck & success;
        
        % PK Parameter are needed
        if ~ TaskList.calculatePKParameter
            for iSet = 1:length(PopRunSet)
                tmp = dir(fullfile('simulations',[PopRunSet(iSet).name '*PK-Analyses.csv']));
                if isempty(tmp)
                    writeToReportLog('ERROR',sprintf('PK Parameter results for "%s" does not exist, please set Task calculatePKParameter to true',...
                        PopRunSet(iSet).name),false);
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
    if ~isempty(VPC) && ~isempty(VPC.PhysProperties)
        for iD = 1:length(VPC.PhysProperties)
            jj = strcmp(VPC.PhysProperties(iD).yList(:,1),'<addOntogenyFactor>');
            if any(jj)
                VPC.PhysProperties(iD).yList = VPC.PhysProperties(iD).yList(~jj,:);
            end
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

    
    writeToReportLog('INFO',sprintf('Input checks were successfully executed. \n)'),false);

return 

function [successInputCheck,VPC] = preparePopulationSimulation(WSettings,PopRunSet,VPC,sensParameterList)

writeToReportLog('INFO',sprintf('Prepare population: %s! \n',PopRunSet.name),false);


% initialize return value
successInputCheck = true;

% create temporary diretory
tmpDir = fullfile('tmp',PopRunSet.name); 
% check if siluation directory already exist
if ~exist(tmpDir,'dir')
    mkdir(tmpDir)
else
    % clear up all temporary results
    if ~WSettings.restart
        delete(fullfile(tmpDir,'*.mat'));
    else
        writeToReportLog('WARNING',sprintf('As script is started in restart mode, temporary results are not deleted! \n'),false);
    end

end

% initialize simulation file
initSimulation(PopRunSet.xml,'none');
simulationIndex=1;


% read population csv
[parPaths,parValues] = readPopulationCSV(PopRunSet.popcsv);

% add aditional populations
if ~isfield(PopRunSet,'addDependentPopulationParameterFh') || isempty(PopRunSet.addDependentPopulationParameterFh)
    addDependentPopulationParameterFh = 'addDependentPopulationParameter';
else
    addDependentPopulationParameterFh = PopRunSet.addDependentPopulationParameterFh;
end
[parPaths,parValues,successTmp] = feval(addDependentPopulationParameterFh,WSettings,parPaths,parValues);
successInputCheck = successInputCheck & successTmp;

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
[ApplicationProtocol,isValid] = getApplicationProtocollFromXML(WSettings,simulationIndex,parPathsStudyDesign,parValuesStudyDesign);     %#ok<ASGLU>
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
                writeToReportLog('WARNING',sprintf('For "%s" no base unit could be identified, may be cause problems if plotted',parPaths{iPar}),...
                    false);
                
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
            writeToReportLog('ERROR',sprintf('Output %s has listed PK Parameter "%s" which is not calculated by the PK parameter function.',...
                OutputList(iO).pathID,OutputList(iO).pKParameterList{1,iPK}),false);
            successInputCheck = false;
        else
            
            [unitFactor,success] = getUnitFactorForUnknownDimension(PKParameterTemplate(jj).unit,OutputList(iO).pKParameterList{2,iPK},MW);
            
            successInputCheck = success & successInputCheck;
           OutputList(iO).pKParameterList{3,iPK} = unitFactor;
        end
    end
        
    
end

save(fullfile(tmpDir,'outputList.mat'),'OutputList');

% check sensitivity definition
if ~isempty(sensParameterList)
    sensParameterList = checkSensitivityParameter(sensParameterList,simulationIndex);     %#ok<NASGU>
    save(fullfile(tmpDir,'sensitivity.mat'),'sensParameterList');
end
    
% check if ontogenyFactors should be added to the physiology plotList
if ~isempty(VPC) && ~isempty(VPC.PhysProperties) 
    
    for iD = 1:length(VPC.PhysProperties)
        if any(strcmp(VPC.PhysProperties(iD).yList(:,1),'<addOntogenyFactor>'))
            ontogenyIDs = getParameter('*|Ontogeny factor*',1,'parametertype','readonly','property','ID');
            jj = getParameter(ontogenyIDs,1,'parametertype','readonly','property','isFormula');
    
            ontogenyPaths = getParameter(ontogenyIDs(~jj),1,'parametertype','readonly','property','Path');
            
            for iO = 1:length(ontogenyPaths)
                tmp = regexp(ontogenyPaths{iO},'\|','split');
                pathID = strjoin(tmp(2:end),'|');
                if ~ismember(pathID,VPC.PhysProperties(iD).yList(:,1))
                    if ismember(tmp{end},{'Ontogeny factor','Ontogeny factor GI'})
                        ontogenyName = strjoin(tmp([end 2:end-1]),' ');
                    else
                        ontogenyName = tmp{end};
                    end
                    VPC.PhysProperties(iD).yList(end+1,:) = {pathID,ontogenyName,'',{},ontogenyName};
                end
            end
        end
    end
end
    
return
    

function successInputCheck = checkPopulationSimulationInput(PopRunSet)

successInputCheck = true;

% check name on special letters
if any(ismember(PopRunSet.name,'./ ?§$%&()[]{}+~*#'))
    writeToReportLog('ERROR',sprintf('Popsetname "%s" contains special signs, do not use them.',PopRunSet.name),false);
    successInputCheck = false;
end

% mandatory inputs
if ~exist(PopRunSet.xml,'file')
    writeToReportLog('ERROR',sprintf('"%s" does not exist',PopRunSet.xml),false);
    successInputCheck = false;
end
if ~exist(PopRunSet.popcsv,'file')
   writeToReportLog('ERROR',sprintf('"%s" does not exist',PopRunSet.popcsv),false);
    successInputCheck = false;
end
                
% optional inputs
if  ~isempty(PopRunSet.studyDesign) && ~exist(PopRunSet.studyDesign,'file') 
    writeToReportLog('ERROR',sprintf('"%s" is given, but does not exist',PopRunSet.studyDesign),false);
    successInputCheck = false;
end
 
return


function  successInputCheck = checkOutputConsistencyForSensitivity(WSettings,PopRunSet)

successInputCheck = true;

if strcmp(WSettings.workflowMode,'pediatric')
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
            writeToReportLog('ERROR',sprintf('Outputs must be the same for a sensitivity analysis!'),false);
            successInputCheck = false;
    end
    
    if successInputCheck
        for iO = 1:length(baseOutputList)
            if ~(length(OutputList(iO).pKParameterList(1,:)) == length(baseOutputList(iO).pKParameterList(1,:))) || ...
                    ~all(strcmp(OutputList(iO).pKParameterList(1,:),baseOutputList(iO).pKParameterList(1,:)))
                writeToReportLog('ERROR',sprintf('PK parameter of Outputs must be the same for a sensitivity analysis!'),false);
                successInputCheck = false;                
            end
        end
    end
    
end
