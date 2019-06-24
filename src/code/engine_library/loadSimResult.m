function [simTime,simValues,pathID,outputUnit,individualIdVector,sensitivities] = loadSimResult(simulationName,iO,pathID,simResultPrefix)
% LOADSIMRESULT load the simulated result for one sepcific output
%
% [simTime,simValues,pathID,outputUnit,individualIdVector] = loadSimResult(simulationName,iO,pathID)
%
%
%       simTime (doublevector)     timevector
%       simValues (double matrix nInd x nTimePOints):   timeprofile  of the modeloutput)
%       individualIdVector (double vector):  vector with the individual ids
%       pathID (string):  pathname of the modeloutput
%       outputUnit (string):   unit of the modeloutput
%       sensitivities (double matrix nPar x nTimePOints): sensitivty vs optimized parameters
%

% Open Systems Pharmacology Suite;  http://open-systems-pharmacology.org

% set dfeualt vlaues
if ~exist('simResultPrefix','var')
    simResultPrefix = '';
end

% initialze return values
simTime = [];
simValues = [];
outputUnit = '';
individualIdVector = [];

% get simulation results as temporary varaibles
fileList = dir(fullfile('tmp',simulationName,sprintf('%ssimResult*.mat',simResultPrefix)));

% convert csv file to temporary variables
nBunch = length(fileList);
if nBunch == 0
    if isempty(simResultPrefix) % defualt population simulation
        nBunch = readPopulationResultfile(simulationName,'Results','');
    else
        nBunch = readPopulationResultfile(simulationName,[upper(simResultPrefix(1)),simResultPrefix(2:end-1)],simResultPrefix);
    end
end

% if no simulation result is found, somethin weird happend
if nBunch == 0
   error('Simulation results expected but not found');
end

for iBunch = 1:nBunch
   
    load(fullfile('tmp',simulationName,sprintf('%ssimResult_%d.mat',simResultPrefix,iBunch)),'SimResult');

    % get all retunr variables from first buch
    if iBunch ==1
        % get Output index
        if ~exist('pathID','var') || isempty(pathID)
            pathID = SimResult.outputPathList{iO};
        else
            iO = find(strcmp(SimResult.outputPathList,pathID));
            if isempty(iO)
                % desierd output does not exist
                return
            end
        end
        
        simTime = SimResult.time;
        individualIdVector = SimResult.individualIdVector;
        outputUnit = SimResult.outputUnit{iO};
        simValues = SimResult.values{iO};
        
        if isfield(SimResult,'sensitivity')
            sensitivities = SimResult.sensitivity{iO};
        else
            sensitivities = [];
        end
        
    else
    % add variables, which changes between bunches
        individualIdVector = [individualIdVector;SimResult.individualIdVector];  %#ok<AGROW>
        simValues = [simValues,SimResult.values{iO}];  %#ok<AGROW>

    end
    
end

clear SimResult

return
