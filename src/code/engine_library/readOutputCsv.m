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

% Open Systems Pharmacology Suite;  http://open-systems-pharmacology.org


% read data
data = readtab(outputcsv,';',0,0,1,0);

% get rows via identifier;
jjPK = strcmp(data{3,1},'pkparameter');
jjPathID = strcmp(data{3,1},'pathID');
jjReportName = strcmp(data{3,1},'reportName');
jjDisplayUnit = strcmp(data{3,1},'displayUnit');
jjDataTpFilter = strcmp(data{3,1},'dataTpFilter');


% get List of PK Parameter
pKidentifier = data{3,2}(jjPK);

for iOx = 3:size(data,2)
    
    iO = iOx-2;
    
    % initialize output
    OutputList(iO) = getDefaultOutput(data{3,iOx}{jjPathID},data{3,iOx}{jjReportName},data{3,iOx}{jjDisplayUnit}); %#ok<AGROW>
    
    % add data timeprofile filter
    if any(jjDataTpFilter)
        OutputList(iO).dataTpFilter = data{3,iOx}{jjDataTpFilter}; %#ok<AGROW>
    end
    
    if ~isempty(pKidentifier)
        % get List of PK Parameter
        tmp = data{3,iOx}(jjPK);
        
        % check for exclusion citeria -1
        jj = ~strcmp(tmp,'-1');
        
        if any(jj)
            OutputList(iO).pKParameterList = [pKidentifier(jj),tmp(jj)]'; %#ok<AGROW>
        end
    end
end
   
return