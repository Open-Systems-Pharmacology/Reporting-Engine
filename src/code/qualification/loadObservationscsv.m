function Data = loadObservationscsv(csvfile, Id)
% LOADOBSERVATIONSCSV load the observation datasets of corresponding Id
%
% Data = loadObservationscsv(csvfile, Id)
%
%
%       Data (structure) contains all relevant information of
%       the Dataset
%       csvfile (string) name of the csv file to be read
%       Id (string) of the Observation DataSet

% Open Systems Pharmacology Suite;  http://open-systems-pharmacology.org

% read data
warning('OFF', 'MATLAB:table:ModifiedAndSavedVarnames');
t = readtable(csvfile, 'Encoding', 'UTF-8');
warning('ON', 'MATLAB:table:ModifiedAndSavedVarnames');

% outputs
varNames = strrep(strrep(t.Properties.VariableDescriptions,'Original column heading: ',''),'''','');
jj = cellfun(@isempty,varNames);
varNames(jj) = t.Properties.VariableNames(jj);
outputPathList = varNames(2:end);
timePathList = varNames{1};

% Rename time
t.Properties.VariableNames{strContains(t.Properties.VariableNames,'Time')} = 'Time';

% split path and unit for variables
for iP=1:length(outputPathList)
    
    tmp=outputPathList{iP};
    ji=strfind(tmp,'[');
    
    if ~isempty(ji)
        outputPathList{iP}=strtrim(tmp(1:ji-1));
        outputUnit{iP}=strtrim(tmp(ji:end)); %#ok<AGROW>
    else
        outputPathList{iP}=tmp;
        outputUnit{iP}=''; %#ok<AGROW>
    end
end

% For time units
timeUnit = strtrim(timePathList(strfind(timePathList,'['):end));

IndividualId = zeros(size(t.Time));
individualIdVector = unique(IndividualId);

% get values for each output

for iO=2:size(t,2)
    values{iO-1}=reshape(t{:,iO},length(t.Time),length(individualIdVector)); %#ok<AGROW>
end

% collect infos in Structure
Data.Id = Id;
Data.name = csvfile;
Data.time = t.Time';
Data.timeUnit = timeUnit;
Data.y = values';
Data.individualIdVector = individualIdVector;
Data.outputPathList = outputPathList;
Data.outputUnit = outputUnit;

