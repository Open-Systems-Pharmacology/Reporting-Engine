function [success,X,filter,Dict] = readNonmemFile(WSettings,dataFile,dataType,filterList)  %#ok<INUSL>
%READNONMEMFILE reads nonmemfile and executes filter,
%
% [success,X,filter,dict] = readNonmemFile(WSettings,dataFile,dataType,filterList) 
% 
% Inputs:
%       WSettings (structure)    definition of properties used in all
%                   workflow functions see GETDEFAULTWORKFLOWSETTINGS
%       dataFile (string) name of nonmemfile
%       dataType (string) type of data, known are 'timeprofile', % toDo:
%       'aggregatedTimeProfiles', 'PKParameter', 'aggregatedPKParameter';
% 
% Outputs:
%       success (boolean) if false read in did not work

% Open Systems Pharmacology Suite;  http://open-systems-pharmacology.org


% initialize return values
X=[];
filter = [];
Dict = []; %#ok<NASGU>
success = true;

% readtable
warning('OFF', 'MATLAB:table:ModifiedAndSavedVarnames');
t = readtable(dataFile{1},'delimiter',' ','FileType','text');
warning('ON', 'MATLAB:table:ModifiedAndSavedVarnames');
header = t.Properties.VariableNames;

% convert content array to column variables
nRow = height(t);
nCol = width(t);
for iCol = 1:nCol
    eval(sprintf('%s = t.%s;',header{iCol},header{iCol}));
end
clear header;
clear t;

% get dictionary
Dict = readDictionary(dataFile{2});

% get mandatory fields and field list for datatype
manadatoryfields = getMandatoryFieldsbyType(Dict,dataType);

% check if necessary filed are available
for iM = 1:size(manadatoryfields,2)
    if ~exist(manadatoryfields{2,iM},'var')
        writeToReportLog('ERROR',sprintf('%s column is missing in %s',manadatoryfields{2,iM},dataFile{1}),false);
        success = false;
    end
end

if ~success 
    return
end


% execute filter
filter = false(nRow,length(filterList));
for iF = 1:length(filterList)
    filter(:,iF) = eval(filterList{iF});
end

% get rid of lines not needed
jj = any(filter,2);
filter = filter(jj,:);

    
% transfer needed variables to structure
for iF = 1:length(Dict)

    if ~exist(Dict(iF).nonmenColumn,'var')
        writeToReportLog('WARNING',sprintf('%s column is missing in %s, but wanted by Dictionary',Dict(iF).nonmenColumn,dataFile{1}),false);
    else
        eval(sprintf('X.(Dict(iF).matlabID) =  %s(jj);',Dict(iF).nonmenColumn));
    end
end

return


function [manadatoryfields] = getMandatoryFieldsbyType(dict,dataType)



switch dataType
    
    case 'timeprofile'
        
        % set madatory fields
        manadatoryfields = {'stud','sid','time','dv'};
        
    otherwise
        error('not implemented yet');
end

[jj,ix] = ismember({dict.matlabID},manadatoryfields);

manadatoryfields(2,:) = {dict(ix(jj)).nonmenColumn};

return


function dict = readDictionary(dataFile)

warning('OFF', 'MATLAB:table:ModifiedAndSavedVarnames');
tTmp = readtable(dataFile,'FileType','text','Delimiter',';');
warning('ON', 'MATLAB:table:ModifiedAndSavedVarnames');
tTmp = tTmp(:,{'matlabID',	'type',	'nonmenColumn',	'nonmemUnit',	'reportName',	'pathID'});

dict = table2struct(tTmp);

return

