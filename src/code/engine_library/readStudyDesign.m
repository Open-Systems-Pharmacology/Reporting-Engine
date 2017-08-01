function [parPaths,parValues,ixChanged] = readStudyDesign(WSettings,studyDesignCsv,parPaths,parValues)
%READSTUDYDESIGN read the study design description and adds this to the population
%
% [parPaths,parValues,ixChanged] = readStudyDesign(WSettings,studyDesignCsv,parPaths,parValues)
% 
% Inputs:
%   - WSettings       (structure) containing global settings see GETDEFAULTWORKFLOWSETTINGS
%   - studyDesignCsv (string)     name of csv, which contains the information of the study design 
%   - parPaths (cellarray):  pathnames of the population parameters
%   - parValues (double matrix):  values of the population parameters and individuals
%
% Outputs:
%   - parPaths (cellarray):  pathnames of the population parameters,
%                   changed acoording to study desing
%   - parValues (double matrix):  values of the population parameters and
%                   individuals, changed acoording to study desing
%   - ixChanged (double vector) index of changed or added parameter

% Open Systems Pharmacology Suite;  http://open-systems-pharmacology.org



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

% get Values for targetParameters
targetParameterValues  = feval(functionHandle,WSettings,targetParameterList,parPaths,parValues,SplitCondition);

% add new path to parValues
% check if path already exist
[jj_target,ix_target] = ismember(targetParameterList,parPaths);
% get indices
ixChanged = [ix_target(jj_target) length(parPaths)+[1:sum(~jj_target)]];

% set values
parPaths(ixChanged) = targetParameterList; 
parValues(:,ixChanged) = targetParameterValues;


return
