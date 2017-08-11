function runMeanModelWorkflow(WSettings,TaskList,MeanModelSet,VPC,Datafiles)
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

% Open Systems Pharmacology Suite;  http://open-systems-pharmacology.org


% try
    %% check optional inputs
    if ~exist('Datafiles','var')
        Datafiles = {};
    end

    %% initialize workflow
    [WSettings] = initializeWorkflow('Mean Model Simulation',WSettings);
    
    %% check if necessary inputs for task are available
    writeToLog(sprintf('Start input checks'),WSettings.logfile,true,false);
    
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
            writeToLog(sprintf('ERROR: VPC definition is missing!'),WSettings.logfile,true,false);
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
        [success,VPC] = prepareMeanModelSimulation(WSettings,MeanModelSet(iSet),VPC);
        successInputCheck = success & successInputCheck;
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

    
    writeToLog(sprintf('Input checks were successfully executed. \n'),WSettings.logfile,true,false);

    %% Tasks processing
    
    % generate VPC
    if TaskList.doVPC
        runMeanModelVPC(WSettings,VPC,MeanModelSet);
    end

    
% catch exception
%     
%     save('exception.mat','exception')
%     writeToLog(exception.message,WSettings.logfile,true,false);
%     
% end



return


function [successInputCheck,VPC] = prepareMeanModelSimulation(WSettings,MeanModelSet,VPC)

% initialize return value
successInputCheck = true;

% create temporary diretory
tmpDir = fullfile('tmp',MeanModelSet.name); 
% check if siluation directory already exist
if ~exist(tmpDir,'dir')
    mkdir(tmpDir)
else
    % clear up all temporary results
    delete(fullfile(tmpDir,'*.mat'));
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

% do the simulation
if ~successInputCheck
    return
end

SimResult = generateSimResult(WSettings,'',simulationIndex,{},[],{OutputList.pathID},1,0); %#ok<NASGU>

% save as temporary file
save(fullfile(tmpDir,'simResult_1.mat'),'SimResult');

% start the caluclation
parPaths = {'Organism|Weight','Organism|Height'};
parValues(1,1) = getParameter('*Organism|Weight',1,'parametertype','readonly');
parValues(1,2) = getParameter('*Organism|Height',1,'parametertype','readonly');
PKPList = calculatesPKParameterList(WSettings,MeanModelSet.name,MeanModelSet.calculatePKParameterFh,parPaths,parValues); %#ok<NASGU>

% save as temporary file
save(fullfile(tmpDir,'pKPList.mat'),'PKPList');
  

return
    

function successInputCheck = checkMeanModelInput(WSettings,MeanModelSet)

successInputCheck = true;

% check name on special letters
if any(ismember(MeanModelSet.name,'./ ?§$%&()[]{}+~*#'))
    writeToLog(sprintf('ERROR: Name "%s" contains special signs, do not use them.',MeanModelSet.name),WSettings.logfile,true,false);
    successInputCheck = false;
end

% mandatory inputs
if ~exist(MeanModelSet.xml,'file')
    writeToLog(sprintf('ERROR: "%s" does not exist',MeanModelSet.xml),WSettings.logfile,true,false);
    successInputCheck = false;
end
                 
return

