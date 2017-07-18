function addDosetablePerWeight(Settings,targetParameterList,parPathes,parValues,SplitCondition)
% add new lines to the population, which sets
%
% Inputs:
%   - Settings       (structure) containing global settings see GETDEFAULTWORKFLOWSETTINGS
%   - targetParameterList (cellarray of string)   List of
%                       applicationparameter, which should be added whit adjusted dose to the
%                       population values
%   - parPathes (cellarray):  pathnames of the population parameters
%   - parValues (double matrix):  values of the population parameters and individuals
%   - SplitCondition structure with following fields
%           - BWmin	 (double vector) [kg] vector which defines lower limit of one body weight bin
%           - BWmax	 (double vector) [kg] vector which defines upper limit of one body weight bin
%           - targetParameter (double vector) values, which the target parameter will be in the corresponding body weight bin
%
% Outputs:
%   - parPathes (cellarray):  pathnames of the population parameters,
%                   changed acoording to study desing
%   - parValues (double matrix):  values of the population parameters and
%                   individuals, changed acoording to study desing

% Open Systems Pharmacology Suite;  support@systems-biology.com
% Date: 18-July-2017


% get Body weight of poulation
jj_BW = strcmp('Organism|Weight',parPathes);
weight =  parValues(:,jj_BW);

%    