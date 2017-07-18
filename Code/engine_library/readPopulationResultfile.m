function SimResult = readPopulationResultfile(simulationName)
%READPOPULATIONRESULTFILE reads a simulation result file exported by PK-SIM
%
% Inputs: 
%   - simulationName (string) name of the simulation 
% Outputs: 
%   - SimResult structure with following fields
%       name (string)   name of simulation
%       time (doublevector)     timevector
%       values (cellarray):  the ith entry contains a double matrix with the
%           concentration values of the ith quantity ( modeloutput)
%       individualIdVector (double vector):  vector with the individual ids
%       outputList (cellarray):  the ith entry contains a string with the
%          pathname of the ith quantity ( modeloutput)
%       outputUnit (cellarray):  the ith entry contains a string with the
%          unit of the ith quantity ( modeloutput)

% Open Systems Pharmacology Suite;  support@systems-biology.com
% Date: 14-July-2017

% getName of the csvfile
csvfile = fullfile('Simulations',[simulationName '-Results.csv']); 

% read header
fid = fopen(csvfile);
tmp=textscan(fid,'%s','delimiter','\n');
C=tmp{1};
fclose(fid);

% outputs
tmp=textscan(C{1},'%s','delimiter',';');
outputList=tmp{1}(3:end);

% split path and unit
for iP=1:length(outputList)
    tmp=outputList{iP};
    ji=strfind(tmp,'[');
    if ~isempty(ji)
        outputList{iP}=strtrim(tmp(1:ji-1)); 
        outputUnit{iP}=strtrim(tmp(ji:end)); %#ok<AGROW>
    else
         outputList{iP}=tmp; 
         outputUnit{iP}=''; %#ok<AGROW>
    end
end

% read numeric data
M=dlmread(csvfile,';',1,0);
individualIdVector = unique(M(:,1),'stable');
time = unique(M(:,2));

% get values for each output
for iO=3:size(M,2)
    values{iO-2}=reshape(M(:,iO)',length(time),length(individualIdVector)); %#ok<AGROW>
end


% collect infos in Structure
SimResult.name = simulationName;
SimResult.time = time;
SimResult.values = values;
SimResult.individualIdVector = individualIdVector;
SimResult.outputList = outputList;
SimResult.outputUnit = outputUnit;


return
