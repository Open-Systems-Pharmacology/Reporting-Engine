function SimResult = loadSimResultcsv(csvfile, Simulation)
% LOADSIMRESULTCSV load the simulation data of corresponding Simulation
%
% SimResult = loadSimResultcsv(csvfile, Simulation)
%
%
%       SimResult  (structure) contains all relevant information of
%       the Simulation Datasets
%       csvfile (string) name of the csv file to be read
%       Simulation (string) Id of the Simulation DataSet

% Open Systems Pharmacology Suite;  http://open-systems-pharmacology.org

% read data
warning('OFF', 'MATLAB:table:ModifiedAndSavedVarnames');
t = readtable(csvfile, 'Encoding', 'UTF-8');
warning('ON', 'MATLAB:table:ModifiedAndSavedVarnames');

% rename individual ID
t.Properties.VariableNames{strContains(t.Properties.VariableNames,'IndividualId')} = 'IndividualId';

% outputs
varNames = strrep(strrep(t.Properties.VariableDescriptions,'Original column heading: ',''),'''','');
jj = cellfun(@isempty,varNames);
varNames(jj) = t.Properties.VariableNames(jj);
outputPathList = varNames(3:end);
timePathList = varNames{2};

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

% read numeric data and remove duplicate points
individualIdVector = unique(t.IndividualId);
time = unique(t{:,2});

% get values for each output

for iO=3:size(t,2)
    values{iO-2}=reshape(t{:,iO},length(time),length(individualIdVector)); %#ok<AGROW>
end

% collect infos in Structure
SimResult.Id = Simulation;
SimResult.name = csvfile;
SimResult.time = time';
SimResult.timeUnit = timeUnit;
SimResult.y = values';
SimResult.individualIdVector = individualIdVector;
SimResult.outputPathList = outputPathList;
SimResult.outputUnit = outputUnit;
SimResult.outputDimension = outputDimension;

SimResult.timeRef=[];
SimResult.yRef=[];
SimResult.yMeanModel=[];
