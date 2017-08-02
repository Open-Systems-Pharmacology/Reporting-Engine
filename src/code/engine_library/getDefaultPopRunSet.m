function [ PopRunSet ] = getDefaultPopRunSet(name,xml,pop,OutputList)
% GETDEFAULTPOPRUNSET defines a  combinations of xml,population, outputs and study design to simulate
%
% Inputs:  
%    - name  (string)    this name is used as directory or file name for the outputs, please avoid special signs like /
%    - xml   (string)    simulation file as base for the simulation
%    - pop   (cell array of strings)    first entry name of csv file, second entry description for reporting
%    - outputList (structure)  list of Outputs see GETDEFAULTOUTPUT
%           
% Outputs:
%    - PopRunSet (structure) with following fields
%       - name (string) as input variable name
%       - xml (string) as input variable xml
%       - OutputList (struture)  List of outputs see GETDEFAULTOUTPUT
%       - studyDesign (string) csv which defines study design like dose
%                   table this is an optional input, it may be empty
%       - dataTpFilter (string) this is an optional input, it may be empty
%               if timeprofile data are available, and a filter per popRunSet is necessary 
%               please enter here the filter conditions. 
%               use nonmem column headers as variable name.
%               e.g. 'AGE>4 & STUD==11111'
%       - isReference (boolean)  marks if it is a refernce population 
%                   (e.g. adult in pediatric or healthy in disease comparison, 
%                   which shall be treated differently in plotting routines
%                   default false;
%       - calculatePKParameter_fh (function handle) 
%               ToDo decribe format of function handle
%               this is an optional input, it may be empty

% Open Systems Pharmacology Suite;  http://forum.open-systems-pharmacology.org
% Date: 14-July-2017


%this name is used as directory or file name for the outputs, please avoid special signs like /
PopRunSet.name = name; 

% simulation file as base for the simulation
PopRunSet.xml = xml;

% population csv
PopRunSet.popcsv = pop{1};
PopRunSet.popReportName = pop{2};

% structure which defines outputs
PopRunSet.OutputList = OutputList;

% optional fields

% study design, like a dose table or application protocol
% this is an optional input, it may be empty
PopRunSet.studyDesign = '';

% if timeprofile data are available, and a filter per popRunSet is necessary 
% please enter here the filter conditions. 
% use nonmem column headers as variable name.
% e.g. 'AGE>4 & STUD==11111'
PopRunSet.dataTpFilter='';

% marks if it is a reference simulation, which is compared in the plots to
% other simulations
PopRunSet.isReference = false;

% for PK Parameter calculation a function handle can be added
% this is an optional input, it may be empty
PopRunSet.calculatePKParameter_fh = '';
