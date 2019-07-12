function Data = loadObservationscsv(csvfile, Id)
% LOADOBSERVATIONSCSV load the observation datasets of corresponding Id
%
% Data = loadObservationscsv(csvfile, Id)
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

% Rename time variable as Time in Observed Dataset,
success_index = (strContains(t.Properties.VariableNames,'Time') | strContains(t.Properties.VariableNames,'time'));
try
    t.Properties.VariableNames{success_index} = 'Time';
catch 
    t.Properties.VariableNames{1} = 'Time';
end

[timePath, timeUnit] = strimPath(timePathList);

for iP=1:length(outputPathList)
    % Get unit and dimension from string within brackets
    [outputPathList{iP}, outputUnit{iP}, outputDimension{iP}] = strimPath(outputPathList{iP});
    
end

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
Data.outputDimension = outputDimension;

% Check and add LLOQ to Data
LLOQcsvfile = replace(csvfile, '.csv', '_LLOQ.csv');
if ~isfile(LLOQcsvfile)
    Data.LLOQ = [];
    Data.LLOQUnit = [];
else
    tLLOQ = readtable(LLOQcsvfile, 'Encoding', 'UTF-8');
    % outputs
    Data.LLOQ = tLLOQ{1};
    Data.LLOQUnit = tLLOQ{2};
end

