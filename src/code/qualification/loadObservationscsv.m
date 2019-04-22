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

% Rename time
% In Observed Dataset,
success_index = (strContains(t.Properties.VariableNames,'Time') | strContains(t.Properties.VariableNames,'time'));
t.Properties.VariableNames{success_index} = 'Time';

% split path and unit for variables
% Get all the units for concentration and amount
allUnits=[getUnitsForDimension('mass') getUnitsForDimension('amount') getUnitsForDimension('concentration')];
% To get the best match bewteen header and unit
for jii=1:length(allUnits)
    LallUnits(jii)=length(allUnits{jii});
end

for iP=1:length(outputPathList)
    
    tmp=outputPathList{iP};
    ji=strfind(tmp,'[');
    
    if ~isempty(ji)
        outputPathList{iP}=strtrim(tmp(1:ji-1));
        % Get time unit with brackets
        bracketedUnit=strtrim(tmp(ji:end)); %#ok<AGROW>
        % Get time unit without brackets
        unbracketedUnits=allUnits(strContains(bracketedUnit, allUnits));
        [Lmax, Imax]=max(LallUnits(strContains(bracketedUnit, allUnits)));
        outputUnit{iP}=unbracketedUnits{Imax};
        
    else
        outputPathList{iP}=tmp;
        outputUnit{iP}=[]; %#ok<AGROW>
    end
end

% Get time unit
allTimeUnits=getUnitsForDimension('time');
% Get time unit with brackets
timeUnit = strtrim(timePathList(strfind(timePathList,'['):end));

% Get actual time unit
timeUnit=allTimeUnits{strContains(timeUnit, allTimeUnits)};

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

