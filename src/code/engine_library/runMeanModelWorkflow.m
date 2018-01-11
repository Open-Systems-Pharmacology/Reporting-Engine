function runMeanModelWorkflow(WSettings,TaskList,MeanModelSet,VPC,Datafiles,sensParameterList,MBS)
%RUNMEANMODELWORKFLOW master routine for Mean Model workflow
%
% runMeanModelWorkflow(WSettings,TaskList,MeanModelSet,VPC,Datafiles,sensParameterList,MBS)
% 
% Inputs:
%       WSettings (structure)    definition of properties used in all
%                   workflow functions see GETDEFAULTWORKFLOWSETTINGS
%       TaskList (structure)    list of task which should be executed
%           TaskList.doVPC (boolean) if true outputs are plotted
%           TaskList.doSensitivityAnalysis (boolean) if true sensitivity analysis is executed
%           TaskList.doAbsorptionPlots (boolean) if true absorption plots are generated
%           TaskList.checkMassbalance (boolean) if truemassbalance plots are generated
%       MeanModelSet (structure)   list of MeanModel simulations
%           - name:  (string) identifier: it is used as directory or file name for the results, please avoid special signs like \
%           - reportName:  (string) used in report figures and tables as description and legend entry
%           - xml: (string) name of the xml file
%           - calculatePKParameterFh: (string) name of the function to calculate the PK Parameter, 
%                   if it is empty the default function, which calculates default OSP Suite PK Parameter is used.
%           - dataTpFilter: (string)  filter for nonmem data file. It must be matlab readable code, 
%                   using the nonmem headers as variable names.  If it is  empty no data are filtered, 
%                   so if you want to use all, include a filter including all data (e.g. SID>0)
%           - dataReportName:  (string) name used in figures for the filtered data
%           - OutputList:  (structure) which describes the properties of the outputs to generate
%               - pathID:  (string) pathname of the output in the mode, e.g. 
%                   “Organism|PeripheralVenousBlood|Theophylline|Plasma (Peripheral Venous Blood)”
%               - reportName: (string) name used in the report for the output
%               - displayUnit:  (string) unit to display the output. 
%                   It must be a unit known to the OSPSuite (case insensitive)
%               - dataTpFilter: (string) filter for nonmem data file. It must be matlab readable code, 
%                   using the nonmem headers as variable names.  If it is  empty no data are filtered, 
%                   so if you want to use all, include a filter including all data (e.g. SID>0)
%               - residualScale: (string) scale of residual calculation, for lin residuals are calculated as data - simulation 
%                   for log as log(data) - log(simulation)
%               - unitFactor: (double)  initialvalue nan, it is calculated internally to 
%                   transfer the internal units of the outputs to the display unit
%               - pKParameterList: (2xn cellarray) first row identifier to select the method 
%                   within the calculate PKParameter function, second row, display unit
%   	Datafiles (cellarray) % set to {}, if no data available, first column name of nonmem file, 
%                   second column dictionary file , third column datatype, so far only implemented type is 'timeprofile' 
%       VPC  (structure) contains information which plots should be
%               generated  see GETDEFAULTVPCSETTINGS
%       sensParameterList  (cellarry) columns: 1. pathid of parameter,
%                       2. number of steps
%                       3. variation range
%                       4. name used in report 
%       MBS (structure) contains information how to generate the
%                   Massbalance plots see GETDEFAULTMASSBALANCESETTINGS


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
        writeToReportLog('WARNING',sprintf('Script is started in restart mode, temporary results are nor deleted! \n'),false);
    end
end


% analyse applicationProtocol
[ApplicationProtocol,isValid] = getApplicationProtocollFromXML(WSettings,MeanModelSet.xml);    %#ok<ASGLU>
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
        
        simValues = SimResult.values{iO};
        jj = simValues~=0;

        
        for iPar = 1:length(VPC.optimizedParameters)
            [~,v(:,iPar)] = getSimulationSensitivityResult(['*|' OutputList.pathID],['*|' VPC.optimizedParameters{iPar}],simulationIndex); %#ok<AGROW>
            % ( d log(y(p))/dp = 1/y(p) * dy(p)/dp
            if strcmp(OutputList(iO).residualScale,'log')
                v(jj,iPar) = v(jj,iPar)./simValues(jj); %#ok<AGROW>
            end
        end
        
        

        
        
        SimResult.sensitivity{iO} = v;
    end
    
end
% save as temporary file
save(fullfile(tmpDir,'simResult_1.mat'),'SimResult');

% start the PKParameter caluclation
parPaths = {'Organism|Weight','Organism|Height'};
if existsParameter('*Organism|Weight',1,'parametertype','readonly')
    parValues(1,1) = getParameter('*Organism|Weight',1,'parametertype','readonly');
else
    parValues(1,1) = nan;
end
if existsParameter('*Organism|Height',1,'parametertype','readonly')
    parValues(1,2) = getParameter('*Organism|Height',1,'parametertype','readonly');
else
    parValues(1,2) = nan;
end
PKPList = calculatesPKParameterList(WSettings,MeanModelSet.name,MeanModelSet.calculatePKParameterFh,parPaths,parValues); %#ok<NASGU>

% save as temporary file
save(fullfile(tmpDir,'pKPList.mat'),'PKPList');


return
