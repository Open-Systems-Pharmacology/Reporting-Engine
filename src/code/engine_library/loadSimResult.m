function [simTime,simValues,pathID,outputUnit,individualIdVector] = loadSimResult(simulationName,iO,pathID)
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
%

% Open Systems Pharmacology Suite;  http://open-systems-pharmacology.org


% initialze return values
simTime = [];
simValues = [];
outputUnit = '';
individualIdVector = [];

% get simulation results as temporary varaibles
fileList = dir(fullfile('tmp',simulationName,sprintf('simResult*.mat')));

% convert csv file to temporary variables
nBunch = length(fileList);
if nBunch == 0
    nBunch = readPopulationResultfile(simulationName);
end

% if no simulation result is found, somethin weird happend
if nBunch == 0
   error('Simulation results expected but not found');
end

for iBunch = 1:nBunch
   
    load(fullfile('tmp',simulationName,sprintf('simResult_%d.mat',iBunch)),'SimResult');

    % get all retunr variables from first buch
    if iBunch ==1
        % get Output index
        if ~exist('pathID','var')
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
        
    else
    % add variables, which changes between bunches
        individualIdVector = [individualIdVector;SimResult.individualIdVector];  %#ok<AGROW>
        simValues = [simValues,SimResult.values{iO}];  %#ok<AGROW>
    end
    
end

clear SimResult

return