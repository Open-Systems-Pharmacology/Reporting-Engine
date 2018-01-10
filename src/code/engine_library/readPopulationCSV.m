function [parPaths,parValues] = readPopulationCSV(csvfile)
%READPOPULATIONCSV reads a poulation csv file
%
% Inputs
%   - csvfile (string) name of the exported file
% Outputs
%   - parPaths (cellarray):  pathnames of the exported parameters
%   - parValues (double matrix):  values of the exported parameters and individuals

% Open Systems Pharmacology Suite;  http://open-systems-pharmacology.org


if ~exist(csvfile,'file')
    error('This file does no exist: %s',csvfile)
end

    
% read table
t = readtable(csvfile);

parPaths = strrep(strrep(t.Properties.VariableDescriptions,'Original column heading: ',''),'''','');
jj = cellfun(@isempty,parPaths);
parPaths(jj) = t.Properties.VariableNames(jj);

% get numeric values
C =  table2cell(t);
jj = cellfun(@(x) isnumeric(x),C);
jjCol = all(jj);
parValues = nan(size(C));
parValues(:,jjCol) = cell2mat(C(:,jjCol));

return
