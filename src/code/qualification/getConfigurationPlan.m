function ConfigurationPlan = getConfigurationPlan(jsonFile)
%GETCONFIGURATIONPLAN Support function:
% Read a .json ConfigurationPlan and output the .json structure and create
% the folder structure for the report
%
%   getConfigurationPlan(jsonFile)
%       jsonFile (string) name .json file where the information of the
%       configuration plan is stored
%

% Open Systems Pharmacology Suite;  http://open-systems-pharmacology.org

%% --------------------------------------------------% 
% Read .json file and output its structure

try
jsonFileContent = fileread(jsonFile);
ConfigurationPlan = jsondecode(jsonFileContent);
writeToReportLog('INFO',[jsonFile ' was extracted successfully'],'true');
catch
writeToReportLog('ERROR',[jsonFile ' could not be extracted'],'true');
end

%% Check if every important parts are within the .json content
% So far mandatory, could be split with Optional
mandatoryFields = {'SimulationMappings', 'ObservedDataSets', 'Plots', 'Inputs', 'Sections'};
for i=1:length(mandatoryFields)
if ~isfield(ConfigurationPlan, mandatoryFields)
    writeToReportLog('WARNING',[mandatoryFields(i) ' is missing from ' jsonFile],'false');
end
end

%% Generate Reporting Engine Output folder structure
% So far, will create a Reporting Engine Output folder in the parent folder
% This can be tuned or can use an input from .json
REOutputFolder = '..\reporting engine output';
mkdir(REOutputFolder);

ConfigurationPlan.Sections = generateOutputFolders(ConfigurationPlan.Sections, REOutputFolder);