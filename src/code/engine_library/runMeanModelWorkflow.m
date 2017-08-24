function runMeanModelWorkflow(WSettings,TaskList,MeanModelSet,VPC,Datafiles,sensParameterList,MBS)
%RUNMEANMODELWORKFLOW master routine for Mean Model workflow
%
% runMeanModelWorkflow(TaskList,WSettings,PopRunSet,varargin)
% 
% Inputs:
%       WSettings (structure)    definition of properties used in all
%                   workflow functions see GETDEFAULTWORKFLOWSETTINGS
%       TaskList (structure)    list of task which should be executed see GENERATEWORKFLOWINPUTFORMEANMODELSIMULATION
%       RUNMEANMODELWORKFLOW (structure)   list of population simulations see GENERATEWORKFLOWINPUTFORMEANMODELSIMULATION
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
    writeToReportLog('INFO','Start input checks',false);
    
    successInputCheck = true;

    % simulation
    for iSet = 1:length(MeanModelSet)
        successInputCheck = checkMeanModelInput(WSettings,MeanModelSet(iSet)) & successInputCheck;
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
           writeToReportLog('ERROR',sprintf('sensitivity definition is missing!'),false);
            successInputCheck = false;
        end
    end

    % PK Parameter
    if TaskList.checkMassbalance
        if isempty(MBS)
            writeToReportLog('ERROR',sprintf('ERROR: MBS definition is missing!'),false);
            successInputCheck = false;
        end
    end
    
    % stop if not all inputs are available
    if ~successInputCheck
        return
    end

    
    % process xmls, check for consistency  and create temporary
    % directory with infos like application protocols, and factor for
    % unitconversion, and timeprofiles of desired outputs 

    for iSet = 1:length(MeanModelSet)
        [success,VPC,simulationIndex] = prepareMeanModelSimulation(WSettings,MeanModelSet(iSet),VPC,sensParameterList);
        successInputCheck = success & successInputCheck;
        
        % processSimulation
        if successInputCheck
            doSim(WSettings,MeanModelSet(iSet),VPC,simulationIndex);
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
                    successInputCheck = readTimeprofileDataForSimulation(WSettings,Datafiles(iData,[1 2]),MeanModelSet);
                otherwise
                    error ('unknown datatype')
            end
        end
    end              
    
    % stop if input is inconsistent
    if ~successInputCheck
        return
    end

    
    writeToReportLog('INFO',sprintf('Input checks were successfully executed. \n'),false);

    %% Tasks processing
    
    % generate VPC
    if TaskList.doVPC
        runMeanModelVPC(WSettings,VPC,MeanModelSet);
    end

    % generate sesnitivity analysis
    if TaskList.doSensitivityAnalysis
        runMeanModelSensitivity(WSettings,MeanModelSet);
    end

    % generate absorption plots
    if TaskList.doAbsorptionPlots
        runMeanModelAbsorption(WSettings,MeanModelSet);
    end

    % checkMassbalance
    if TaskList.checkMassbalance
        runMeanModelCheckMassbalance(WSettings,MeanModelSet,MBS);
    end
    
    writeToReportLog('INFO',sprintf('Mean model workflow finalized \n'),false);

    
catch exception
    
    save(sprintf('exception_%s.mat',datestr(now,'ddmmyy_hhMM')),'exception');
    writeToReportLog('ERROR',exception.message,false);
    writeToReportLog('INFO',sprintf('Mean model workflow finished with error \n'),false);
  
end


return


function [successInputCheck,VPC,simulationIndex] = prepareMeanModelSimulation(WSettings,MeanModelSet,VPC,sensParameterList)

% initialize return value
successInputCheck = true;

% create temporary diretory
tmpDir = fullfile('tmp',MeanModelSet.name); 
% check if siluation directory already exist
if ~exist(tmpDir,'dir')
    mkdir(tmpDir)
else
    % clear up all temporary results
    if ~WSettings.restart
        delete(fullfile(tmpDir,'*.mat'));
    else
        writeToReportLog('WARNING',sprintf('As script is started in restart mode, temporaray results are nor deleted! \n'),false);
    end
end


% analyse applicationProtocol
[ApplicationProtocol,isValid] = getApplicationProtocollFromXML(WSettings,MeanModelSet.xml);   %#ok<NASGU>
save(fullfile(tmpDir,'applicationProtocol.mat'),'ApplicationProtocol','isValid');

simulationIndex = 1; %Initialised in getApplicationProtocollFromXML

if isempty(MeanModelSet.calculatePKParameterFh)
    PKParameterTemplate = calculatePKParameterForApplicationProtocol(WSettings,ApplicationProtocol);
else
    PKParameterTemplate = feval(MeanModelSet.calculatePKParameterFh,WSettings,ApplicationProtocol);
end

save(fullfile(tmpDir,'applicationProtocol.mat'),'-append','PKParameterTemplate');


% analyse outputs, check if existing and add unitfactor
OutputList = MeanModelSet.OutputList;
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


% do the simulation
if ~successInputCheck
    return
end

% if no outputs are defined, return now
if isempty(OutputList)
    return
end

return
    

function successInputCheck = checkMeanModelInput(~,MeanModelSet)

successInputCheck = true;

% check name on special letters
if any(ismember(MeanModelSet.name,'./ ?§$%&()[]{}+~*#'))
    writeToReportLog('ERROR',sprintf('Name "%s" contains special signs, do not use them.',MeanModelSet.name),false);
    successInputCheck = false;
end

% mandatory inputs
if ~exist(MeanModelSet.xml,'file')
    writeToReportLog('ERROR',sprintf('"%s" does not exist',MeanModelSet.xml),false);
    successInputCheck = false;
end
                 
return

function doSim(WSettings,MeanModelSet,VPC,simulationIndex)

% do the simualtion of the model
tmpDir = fullfile('tmp',MeanModelSet.name); 
load(fullfile(tmpDir,'outputList.mat'),'OutputList');

if ~isempty(VPC) && ~isempty(VPC.optimizedParameters)
    initStruct = [];
    for iPar = 1:length(VPC.optimizedParameters)
        initStruct = initParameter(initStruct,['*|' VPC.optimizedParameters{iPar}],'always','calculateSensitivity',true);
    end
    
    initSimulation(MeanModelSet.xml,initStruct);
    
    simulationIndex=1;
    
end

SimResult = generateSimResult(WSettings,'',simulationIndex,{},[],{OutputList.pathID},1,0); 

if ~isempty(VPC) && ~isempty(VPC.optimizedParameters)

    for iO = 1:length(OutputList)
        for iPar = 1:length(VPC.optimizedParameters)
            [~,v(:,iPar)] = getSimulationSensitivityResult(['*|' OutputList.pathID],['*|' VPC.optimizedParameters{iPar}],simulationIndex); %#ok<AGROW>
        end
        SimResult.sensitivity{iO} = v;
    end
    
end
% save as temporary file
save(fullfile(tmpDir,'simResult_1.mat'),'SimResult');

% start the PKParameter caluclation
parPaths = {'Organism|Weight','Organism|Height'};
parValues(1,1) = getParameter('*Organism|Weight',1,'parametertype','readonly');
parValues(1,2) = getParameter('*Organism|Height',1,'parametertype','readonly');
PKPList = calculatesPKParameterList(WSettings,MeanModelSet.name,MeanModelSet.calculatePKParameterFh,parPaths,parValues); %#ok<NASGU>

% save as temporary file
save(fullfile(tmpDir,'pKPList.mat'),'PKPList');


return