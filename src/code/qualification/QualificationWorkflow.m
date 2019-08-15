% Script to perform a Qualification Plan workflow
% Qualification Plan Workflow developed with Matlab 2017b
% --------------------------------------------------------------

close all
tic

% --------------------------------------------------------------
% replace REInput_path, REOutput_path and jsonFile with your paths
%
% - REInput_path: input path for the Reporting engine 
%                 (corresponds to the output path defined in the Qualification Runner)
%
% - REOutput_path: outputs of the Reporting Engine will be created here
%                  CAUTION: if the folder is not empty, its contents will be deleted
%
% - jsonFile: report configuration file defined in the Qualification Runner
%             (default is 'report-configuration-plan.json')
REInput_path = 'C:\re_input';
REOutput_path = 'C:\re_output';
jsonFile = 'report-configuration-plan.json';

% --------------------------------------------------------------
% Get the Configuration Plan Settings
[WSettings, ConfigurationPlan, TaskList, ObservedDataSets] = initializeQualificationWorkflow(jsonFile, REInput_path, REOutput_path);

%OPTIONAL: set watermark. If set, it will appear in all generated plots
%WSettings.Watermark = 'Preliminary';

% --------------------------------------------------------------
% run the Worklfow tasklist of ConfigurationPlan
runQualificationWorkflow(WSettings, ConfigurationPlan, TaskList, ObservedDataSets);

QualificationWorkflowTime = toc/60;
fprintf('\n Qualification Workflow Duration: %0.1f minutes \n', QualificationWorkflowTime);

% --------------------------------------------------------------
%OPTIONAL: call MarkdownJoiner to combine Reporting Engine output into the final report
% (this can be also done later from Windows command line)
%
% replace ReportOutput_path and MarkdownJoiner_path with your paths
%
% - ReportOutput_path: final report will be generated here
%
% - MarkdownJoiner_path: location of markdown-joiner.exe

ReportOutput_path='C:\report';
MarkdownJoiner_path='C:\MD\markdown-joiner.exe';

% alternative #1: ReportOutput_path must be empty. If not, report generation will fail
system([MarkdownJoiner_path ' -i ' REOutput_path ' -o ' ReportOutput_path]);

% alternative #2: (CAUTION) ReportOutput_path will be cleared first
%system([MarkdownJoiner_path ' -i ' REOutput_path ' -o ' ReportOutput_path ' -f']);

mergeQualificationMarkdown([ReportOutput_path filesep 'markdown_for_pdf'], [ReportOutput_path filesep 'report_merged.md']);
