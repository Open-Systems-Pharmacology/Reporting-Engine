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
catch exception
writeToReportLog('WARNING', sprintf('Error in TimeProfile plot %d. \n %s \n', csvfile, exception.message), 'true', exception);
warning('Variable Time not found in first column header of observed dataset, %s. \n First column assumed as Time variable \n', csvfile, exception.message);
t.Properties.VariableNames{1} = 'Time';
end

% Split path and unit for variables
[unitList,unitList_dimensionList]=iniUnitList(0);

for iP=1:length(outputPathList)
    
    % Get unit and dimension from string within brackets
    tmp=outputPathList{iP};
    ji=strfind(tmp,'[');
    ij=strfind(tmp,']');
    
    if ~isempty(ji)
        % In case more than one '[' ']' are found, takes the first
        ji=ji(1); ij=ij(1);
        
        % Split between header and unit
        outputPathList{iP}=strtrim(tmp(1:ji-1));
        
        % Get unit within brackets
        tmpUnit=strtrim(tmp(ji+1:ij-1)); %#ok<AGROW>
        
        for jiUnit=1:length(unitList)
            iUnit = find(strcmpi(unitList(jiUnit).unit_txt,tmpUnit));
            if ~isempty(iUnit)
                % Get the unit and dimension as cell
                outputUnit(iP)=unitList(jiUnit).unit_txt(iUnit);
                outputDimension(iP)=unitList_dimensionList(jiUnit);
                break
            end
        end
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
Data.outputDimension = outputDimension;

