function [success,X,filter,dict] = readNonmemFile(Settings,dataFile,dataType,filterList) 
%READNONMEMFILE reads nonmemfile and executes filter,
%
% Inputs:
%       Settings (structure)    definition of properties used in all
%                   workflow functions see GETDEFAULTWORKFLOWSETTINGS
%       dataFile (string) name of nonmemfile
%       dataType (string) type of data, known are 'timeprofile', % toDo:
%       'aggregatedTimeProfiles', 'PKParameter', 'aggregatedPKParameter';
% 
% Outputs:
%       success (boolean) if false read in did not work

% Open Systems Pharmacology Suite;  http://forum.open-systems-pharmacology.org
% Date: 21-July-2017

% initialize return values
X=[];
filter = [];
dict = {};

% readtable
[content,header,success] = readtable(Settings,dataFile{1});
if ~success
    return
end

% convert content array to column variables
nRow = size(content,1);
nCol = size(content,2);
for iCol = 1:nCol
    % check if column is numeric
    contentnum=cellfun(@str2double,content(:,iCol));
    jj = isnan(contentnum);
    
    % less string then numrich values
    if sum(jj) < nRow/2
        tmp = contentnum;
    else
        tmp = content(:,iCol); %#ok<NASGU>
    end
    
    % transfer to variable
    eval(sprintf('%s = tmp;',header{iCol}));
end
clear content;
clear tmp;

% get dictionary
dict = readDictionary(dataFile{2});

% get mandatory fields and field list for datatype
manadatoryfields = getMandatoryFieldsbyType(dict,dataType);

% check if necessary filed are available
for iM = 1:size(manadatoryfields,2)
    if ~exist(manadatoryfields{2,iM},'var')
        writeToLog(sprintf('%s column is missing in %s',manadatoryfields{2,iM},dataFile{1}),Settings.logfile,true,false);
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
for iF = 1:length(dict)

    if ~exist(dict(iF).nonmenColumn,'var')
        writeToLog(sprintf('%s column is missing in %s, but wanted by Dictionary',dict(iF).nonmenColumn,dataFile),Settings.logfile,true,false);
    else
        eval(sprintf('X.(dict(iF).matlabID) =  %s(jj);',dict(iF).nonmenColumn));
    end
end

return


function [content,header,success] = readtable(Settings,dataFile)

success = true;

% read table 
fid = fopen(dataFile);

% get number of columns
lcount=1;
line1=strtrim(fgetl(fid));
delim = ' ';
ncol=length(strfind(line1,delim))+1;

% get number of lines to allocate space, which saves time
while true
    line=fgetl(fid);
    if line==-1
        break
    end
    lcount=lcount+1;
end


% sets the file position indicator to the beginning of the file
frewind(fid);

% allocate space
content=cell(lcount-1,ncol);

% read in
for  iLine =1:lcount
    
    if mod(iLine,1000)==0;
        disp(sprintf('%d lines of %d read in of file: %s',iLine,lcount,dataFile));
    end
    
    % get the text of fthe line
    line=fgetl(fid);
    if line==-1
        break;
    end
    line=strtrim(line);
    
    % split line by delimiter
    tmp = regexp(line,delim,'split');
    
    ncoltmp=length(tmp);
    if isnan(ncol)
        ncol=ncoltmp;
    else
        if ncoltmp~=ncol
            writeToLog(sprintf('Line %d: inconsistent number of columns in %s',iLine,dataFile),Settings.logfile,true,false);
            success = false;
            return
        end
    end
    
    if iLine == 1
        % save header and get rid of beginnign #
        header = strrep(tmp,'#','');
        
    else
        content(iLine-1,:)=tmp;
    end
end

fclose(fid);

writeToLog(sprintf('Read  file %s with %d rows',dataFile,lcount),Settings.logfile,true,false);

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

data  = readtab(dataFile,';',0,0,1,0);

jj = ismember(data(1,:),{'matlabID',	'type',	'nonmenColumn',	'nonmemUnit',	'reportName',	'pathID'});

tmp = cell(length(data{3,1})+1,sum(jj));

k=0;
for iCol = find(jj)
    k=k+1;
    tmp(:,k) = [data(1,iCol); data{3,iCol}];
end

dict = cell2struct(tmp(2:end,:),tmp(1,:),2);

return

