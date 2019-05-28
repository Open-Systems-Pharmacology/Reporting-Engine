function [WSettings, ConfigurationPlan, TaskList, ObservedDataSets] = initializeQualificationWorkflow(jsonFile, REInput_path, REOutput_path)
% function that initialise the workflow, sets up the logfile, and sets global settings
%
% [WSettings, ConfigurationPlan, TaskList, ObservedDataSets] = initializeQualificationWorkflow(jsonFile)
% 
% Inputs:
%       jsonFile (string)  name of configuration .json file
%       REInput_path (string) name of the input RE folder
%       REOutput_path (string) name of the output RE folder
% Outputs:
%       WSettings (structure)  see GETDEFAULTWORKFLOWSETTINGS
%       ConfigurationPlan (structure): Meta structure with all the plots and
%       mappings information
%       TaskList (cells): list of tasks to be performed
%       ObservedDataSets (structure): extracted observations
%

% Open Systems Pharmacology Suite;  http://open-systems-pharmacology.org

% Create output directory
if ~isempty(dir(REOutput_path))
    % Override RE Output if already created
    rmdir(REOutput_path, 's');
    writeToReportLog('WARNING', sprintf('In initializeQualificationWorkflow, Reset Output RE folder "%s" \n', ...
                REOutput_path), 'true');
end
mkdir(REOutput_path);

% Get complete for REOutput and REInput paths
REOutput_structure = what(REOutput_path);
REInput_structure = what(REInput_path);

REOutput_path = REOutput_structure.path;
REInput_path = REInput_structure.path;

% Read the .json file
ConfigurationPlan = getConfigurationPlan(fullfile(REInput_path, jsonFile));

% Set Reporting Engine Input Path 
ConfigurationPlan.REInput_path=REInput_path;

% Linearize the Sections input and create the folder structure
ConfigurationPlan.Sections = generateOutputFolders(ConfigurationPlan.Sections, REInput_path, REOutput_path);

% --------------------------------------------------------------
% Copy Intro and Input content into right section
if isfield(ConfigurationPlan, 'Intro')
    for i=1:length(ConfigurationPlan.Intro)
        try
            copyfile(fullfile(REInput_path, ConfigurationPlan.Intro(i).Path), fullfile(REOutput_path));
        catch exception
            writeToReportLog('ERROR', sprintf('Error in copying Intro: \n %s \n Intro index: %d \n Input path: %s \n Output path: %s \n', ...
                exception.message, i, fullfile(REInput_path, ConfigurationPlan.Intro(i).Path), fullfile(REOutput_path)), 'true', exception);
            rethrow(exception);
        end
    end
end
for i=1:length(ConfigurationPlan.Inputs)
    [SectionPath, indexed_item] = getSection(ConfigurationPlan.Sections, ConfigurationPlan.Inputs(i).SectionId);
    try
        copyfile(fullfile(REInput_path, ConfigurationPlan.Inputs(i).Path), SectionPath);
    catch exception
        writeToReportLog('ERROR', sprintf('Error in copying Inputs: \n %s \n Inputs index: %d \n Input path: %s \n Output path: %s \n', ...
                exception.message, i, fullfile(REInput_path, ConfigurationPlan.Inputs(i).Path), SectionPath), 'true', exception);
            rethrow(exception);
    end
end

% --------------------------------------------------------------
% Global settings
% Not developed yet: uses getDefaultWorkflowSettings from Reporting Engine
WSettings = getDefaultWorkflowSettings('Qualification','default'); close;

% --------------------------------------------------------------
% Load the Observations Data Sets from csv into structures
% Convert the input into cell array
if ~iscell(ConfigurationPlan.ObservedDataSets)
    ConfigurationPlan.ObservedDataSets=num2cell(ConfigurationPlan.ObservedDataSets);
end
if isempty(ConfigurationPlan.ObservedDataSets)
    % If empty input, output ObservedDataSets is an empty variable
    ObservedDataSets=[];
else
    for i=1:length(ConfigurationPlan.ObservedDataSets)
        % For DDI and PK Ratio plots, Observed Data Type is different
        % ObservedDataSets format to be determined
        % Temporary format: structure where the observed data is a table y
        if isfield(ConfigurationPlan.ObservedDataSets{i}, 'Type')
            ObsTable = readtable(fullfile(REInput_path, ConfigurationPlan.ObservedDataSets{i}.Path), 'Encoding', 'UTF-8');
            ObservedDataSets(i).Id=ConfigurationPlan.ObservedDataSets{i}.Id;
            ObservedDataSets(i).Type=ConfigurationPlan.ObservedDataSets{i}.Type;
            ObservedDataSets(i).name=fullfile(REInput_path, ConfigurationPlan.ObservedDataSets{i}.Path);
            ObservedDataSets(i).y=ObsTable;
        else
            ObservedDataSets(i) = loadObservationscsv(fullfile(REInput_path, ConfigurationPlan.ObservedDataSets{i}.Path),...
                ConfigurationPlan.ObservedDataSets{i}.Id);
        end
    end
end

% --------------------------------------------------------------
% Task list: fields in ConfigurationPlan.Plots
TaskList = fields(ConfigurationPlan.Plots);

% --------------------------------------------------------------
% Add the units to unitList.mat if not present
setAmountMassUnitList