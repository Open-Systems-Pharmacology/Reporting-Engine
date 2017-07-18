function [ PopRunSet ] = getDefaultPopRunSet(name,xml,pop,outputList)
% GETDEFAULTPOPRUNSET defines a  combinations of xml,population, outputs and study design to simulate
%
% Inputs:  
%       - name  (string)    this name is used as directory or file name for the outputs, please avoid special signs like /
%       - xml   (string)    simulation file as base for the simulation
%       - outputList (string)  csv to define outputs  
% Outputs:
%       - PopRunSet (structure) with following fields
%           name (string) as input variable name
%           xml (string) as input variable xml
%           outputList (string) as input variable outputList
%           studyDesign (string) csv which defines study design like dose
%                   table this is an optional input, it may be empty
%           calculatePKParameter_fh (function handle) ToDo decribe format
%           of function handle
%               this is an optional input, it may be empty

% Open Systems Pharmacology Suite;  support@systems-biology.com
% Date: 14-July-2017


%this name is used as directory or file name for the outputs, please avoid special signs like /
PopRunSet.name = name; 
% simulation file as base for the simulation
PopRunSet.xml = xml;
% population csv
PopRunSet.pop = pop;
% csv to define outputs
PopRunSet.outputList = outputList;

% study design, like a dose table or application protocol
% this is an optional input, it may be empty
PopRunSet.studyDesign = '';

% for PK Parameter calculation a function handle can be added
% this is an optional input, it may be empty
PopRunSet.calculatePKParameter_fh = '';
