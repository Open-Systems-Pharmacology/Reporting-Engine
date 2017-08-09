function OutputList = readOutputXls(outputXls,outputsheet)
%READOUTPUTXLS reads the output csv and translate it to an array of structures
% see GETDEFAULTOUTPUT
%
% OutputList = readOutputXls(outputcsv)
%
% Inputs:
%   outputXls (string) name of xls file
%   outputsheet (string) name of sheet in slxs file, if empty first sheet is taken
% Outputs
%   OutputList (array of structure)  see GETDEFAULTOUTPUT

% Open Systems Pharmacology Suite;  http://open-systems-pharmacology.org


if ~exist('outputsheet','var') || isempty(outputsheet)
    outputsheet = 1;
end

% read data
[~,~,xlsdata] = xlsread(outputXls,outputsheet);

xlsIdentifier = xlsdata(2:end,1);
xlsPKIdentifier = xlsdata(2:end,2);
xlsInfo = xlsdata(2:end,3:end);
% remove empty rows
jj = ~cellfun(@isnumeric,xlsIdentifier) & ~cellfun(@isempty,xlsIdentifier);
xlsIdentifier = xlsIdentifier(jj);
xlsPKIdentifier = xlsPKIdentifier(jj);
xlsInfo = xlsInfo(jj,:);


% get rows via identifier;
jjPK = strcmp(xlsIdentifier,'pkparameter');
jjPathID = strcmp(xlsIdentifier,'pathID');
jjReportName = strcmp(xlsIdentifier,'reportName');
jjDisplayUnit = strcmp(xlsIdentifier,'displayUnit');
jjDataTpFilter = strcmp(xlsIdentifier,'dataTpFilter');
jjResidualScale = strcmp(xlsIdentifier,'residualScale');


% get cols for output info
jjCol =  ~cellfun(@isempty,xlsInfo(strcmp(xlsIdentifier,'pathID'),:)) & ...
     ~cellfun(@isnumeric,xlsInfo(strcmp(xlsIdentifier,'pathID'),:));


% get List of PK Parameter
pKidentifier = xlsPKIdentifier(jjPK);

iO = 0;
for iOx = find(jjCol)
    
    iO = iO+1;
    
    tmp = xlsInfo(:,iOx);
    
    % set nan to empty
    ij = find(cellfun(@isnumeric,tmp));
    jj = cellfun(@isnan,tmp(ij));
    tmp(ij(jj)) = {''};
    
    % initialize output
    OutputList(iO) = getDefaultOutput(tmp{jjPathID},tmp{jjReportName},tmp{jjDisplayUnit}); %#ok<AGROW>
    
    % add data timeprofile filter
    if any(jjDataTpFilter)
        OutputList(iO).dataTpFilter = tmp{jjDataTpFilter}; %#ok<AGROW>
    end
    if any(jjResidualScale)
        OutputList(iO).residualScale = tmp{jjResidualScale}; %#ok<AGROW>
    else
        OutputList(iO).residualScale ='log'; %#ok<AGROW>
    end
    
    if ~isempty(pKidentifier)
        % get List of PK Parameter
        tmpPK = tmp(jjPK);
        
        % check for exclusion citeria -1
        ij = find(cellfun(@isnumeric,tmpPK));
        jj = cellfun(@(x) x==-1,tmpPK(ij));
        keep = setdiff(1:length(pKidentifier),ij(jj));
        
        if ~isempty(keep)
            OutputList(iO).pKParameterList = [pKidentifier(keep),tmpPK(keep)]'; %#ok<AGROW>
        end
    end
end
   
return
