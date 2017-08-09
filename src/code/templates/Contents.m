% TEMPLATES functions to create support workflow scripts creations
%
% master functions to create workflow script
%   prepareMeanModelWorkflow  - creates out of xls table workflow script
%   preparePopulationWorkflow - creates out of xls table workflow script
%
% auxiliary functions   
%   readOutputXls             - reads the output csv and translate it to an array of structures
%   readWorkflowInput         - converts workflow xls to STructures
%   writeDataFiles            - write datalist to workflow.m
%   writeTaskList             - write Tasklist to workflow.m
%   writeWorkflowHeader       - write header of a wrokflow script
%
% xls files with templates
%   WorkflowInput.xlsx   - template for input for mastr script       
%   Output.xlsx          - template for output.xslx linked in workflow.xlsx       
%   NonmemFieldDicitonary.xlsx  - template for nonmem dictionary
%   StudyDesign.xlsx      - template to create a studydesign.csv
