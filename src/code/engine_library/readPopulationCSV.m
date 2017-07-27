function [parPaths,parValues] = readPopulationCSV(csvfile)
%READPOPULATIONCSV reads a poulation csv file 
%
% Inputs 
%   - csvfile (string) name of the exported file
% Outputs 
%   - parPaths (cellarray):  pathnames of the exported parameters
%   - parValues (double matrix):  values of the exported parameters and individuals

% Open Systems Pharmacology Suite;  http://forum.open-systems-pharmacology.org
% Date: 14-July-2017

if ~exist(csvfile,'file')
    error('This file does no exist: %s')
end

% read data
data = readtab(csvfile,';',0,0,1,0);
parPaths = data(1,:);
parValues = nan(length(data{3,1}),size(parPaths,1));

% get numeric columns
jj = strcmp(data(2,:),'double');
parValues(:,jj) = cell2mat(data(3,jj));


return
