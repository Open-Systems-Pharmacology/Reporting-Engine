function OutputList = readOutputCsv(outputcsv)
% READOUTPUTCSV reads the output csv and translate it to an array of structures
% see GETDEFAULTOUTPUT
%
% OutputList = readOutputCsv(outputcsv)
%
% Inputs:
%   outputcsv (string) name of csv file
% Outputs
%   OutputList (array of structure)  see GETDEFAULTOUTPUT

% Open Systems Pharmacology Suite;  http://forum.open-systems-pharmacology.org
% Date: 31-July-2017


% read data
data = readtab(outputcsv,';',0,0,0,0);

% get List of PK Parameter
pKidentifier = data{3,1}(4:end);

for iOx = 2:size(data,2)
    
    iO = iOx-1;
    
    % initialize output
    OutputList(iO) = getDefaultOutput(data{3,iOx}{1},data{3,iOx}{2},data{3,iOx}{3}); %#ok<AGROW>
    
    if ~isempty(pKidentifier)
        % get List of PK Parameter
        tmp = data{3,iOx}(4:end);
        
        % chekc for exclusion citeria -1
        jj = ~strcmp(tmp,'-1');
        
        if any(jj)
            OutputList(iO).pKParameterList = [pKidentifier(jj),tmp(jj)]'; %#ok<AGROW>
        end
    end
end
   
return