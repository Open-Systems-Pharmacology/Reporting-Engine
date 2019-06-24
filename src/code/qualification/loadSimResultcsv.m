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

[timePath, timeUnit] = strimPath(timePathList);

for iP=1:length(outputPathList)
    % Get unit and dimension from string within brackets
    [outputPathList{iP}, outputUnit{iP}, outputDimension{iP}] = strimPath(outputPathList{iP});
    
end

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
