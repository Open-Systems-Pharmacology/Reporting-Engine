function [parPathes,parValues] = readStudyDesign(Settings,studyDesignCsv,parPathes,parValues)
% read the study design description and adds new pathes
%
% Inputs:
%   - Settings       (structure) containing global settings see GETDEFAULTWORKFLOWSETTINGS
%   - studyDesignCsv (string)     name of csv, which contains the information of the study
%                           design (ToDo documentations or templates)
%   - parPathes (cellarray):  pathnames of the population parameters
%   - parValues (double matrix):  values of the population parameters and individuals
%
% Outputs:
%   - parPathes (cellarray):  pathnames of the population parameters,
%                   changed acoording to study desing
%   - parValues (double matrix):  values of the population parameters and
%                   individuals, changed acoording to study desing

% Open Systems Pharmacology Suite;  support@systems-biology.com
% Date: 18-July-2017


% chekc if file exist
if ~exist(studyDesignCsv,'file')
    error('This file does no exist: %s')
end
    
% read non numeric content of file to variable C
fid = fopen(studyDesignCsv);
tmp=textscan(fid,'%s','delimiter','\n');
C=tmp{1};
fclose(fid);

% read numeric content of file to variable M
M=dlmread(studyDesignCsv,';',3,0);

% get function handle
if isempty(strfind(C{1},'functionHandle'))
    error('first line  of study design csv should begin with: "functionHandle = @"');
end
eval(C{1});
% get targetparameterList
if isempty(strfind(C{2},'targetParameterList'))
    error('second line  of study design csv should begin with: "targetParameterList = {"');
end
eval(C{2});


% Header
header = regexp(C{3},';','split');
SplitCondition = cell2struct(num2cell(M),header,2);

[parPathes,parValues]  = feval(functionHandle,Settings,targetParameterList,parPathes,parValues,SplitCondition);

return
