function [WSettings, ConfigurationPlan, TaskList, ObservedDataSets] = initializeQualificationWorkflow(jsonFile, REOutput_path)
% function that initialise the workflow, sets up the logfile, and sets global settings
%
% [WSettings, ConfigurationPlan, TaskList, ObservedDataSets] = initializeQualificationWorkflow(jsonFile)
% 
% Inputs:
%       jsonFile (string)  name of configuration .json file
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
mkdir(REOutput_path);

% Read the .json file
ConfigurationPlan = getConfigurationPlan(jsonFile);

% Linearize the Sections input and create the folder structure
ConfigurationPlan.Sections = generateOutputFolders(ConfigurationPlan.Sections, REOutput_path);

% --------------------------------------------------------------
% Copy Intro and Input content into right section
if isfield(ConfigurationPlan, 'Intro')
    for i=1:length(ConfigurationPlan.Intro)
        try
            copyfile(ConfigurationPlan.Intro(i).Path, fullfile(REOutput_path));
        catch exception
            writeToReportLog('ERROR', sprintf('Error in copying Intro: \n %s \n Intro index: %d \n Input path: %s \n Output path: %s \n', ...
                exception.message, i, ConfigurationPlan.Intro(i).Path, fullfile(REOutput_path)), 'true', exception);
            rethrow(exception);
        end
    end
end
for i=1:length(ConfigurationPlan.Inputs)
    [SectionPath, indexed_item] = getSection(ConfigurationPlan.Sections, ConfigurationPlan.Inputs(i).SectionId);
    try
        copyfile(ConfigurationPlan.Inputs(i).Path, SectionPath);
    catch exception
        writeToReportLog('ERROR', sprintf('Error in copying Inputs: \n %s \n Inputs index: %d \n Input path: %s \n Output path: %s \n', ...
                exception.message, i, ConfigurationPlan.Inputs(i).Path, SectionPath), 'true', exception);
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
            ObsTable = readtable(ConfigurationPlan.ObservedDataSets{i}.Path, 'Encoding', 'UTF-8');
            ObservedDataSets(i).Id=ConfigurationPlan.ObservedDataSets{i}.Id;
            ObservedDataSets(i).Type=ConfigurationPlan.ObservedDataSets{i}.Type;
            ObservedDataSets(i).name=ConfigurationPlan.ObservedDataSets{i}.Path;
            ObservedDataSets(i).y=ObsTable;
        else
            ObservedDataSets(i) = loadObservationscsv(ConfigurationPlan.ObservedDataSets{i}.Path,...
                ConfigurationPlan.ObservedDataSets{i}.Id);
        end
    end
end

% --------------------------------------------------------------
% Task list: fields in ConfigurationPlan.Plots
TaskList = fields(ConfigurationPlan.Plots);
