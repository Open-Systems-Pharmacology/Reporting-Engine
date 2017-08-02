function PKPList = readPKParameter(simulationName)
%READPKPARAMETER reads a PK Analysis result file 
%
% Inputs: 
%   - simulationName (string) name of the simulation 
% Outputs: 
%   - PKPList cellarray of structures with following fields
%       name (string)   name of pkparameter
%       unit (string)   unit of pkparameter
%       value (double vector): individual values of PK Parameter 

% Open Systems Pharmacology Suite;  http://forum.open-systems-pharmacology.org
% Date: 14-July-2017

% getName of the csvfile
csvfile = fullfile('simulations',[simulationName '-PK-Analyses.csv']); 

load(fullfile('tmp',simulationName,'outputList.mat'),'OutputList');

% read data
data = readtab(csvfile,';',0,0,1,0);

outputPathList = data{3,2};
PKParameter = data{3,3};
value = data{3,4};
unit = data{3,5};

for iO = 1:length(OutputList)
    
    jjO = strcmp(OutputList(iO).pathID,outputPathList);
    
    tmp = [];
    for iPK = 1:length(OutputList(iO).pKParameterList)
       
        jjPK = jjO & strcmp(OutputList(iO).pKParameterList{1,iPK},PKParameter);
        
        if any(jjPK)
            ij = find(jjPK,1);
            
            tmp(iPK).name = PKParameter{ij};
            tmp(iPK).unit = unit{ij};
            tmp(iPK).value = value(jjPK)';
        end
    end

    PKPList{iO} = tmp; %#ok<AGROW>
end

% save as temporary file
save(fullfile('tmp',simulationName,'pKPList.mat'),'PKPList');
return