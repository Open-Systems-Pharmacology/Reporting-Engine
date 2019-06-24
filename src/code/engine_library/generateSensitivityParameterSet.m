function  [parPathsNew,parValuesNew,SensPointer] = generateSensitivityParameterSet(WSettings,sensParameterList,parPaths,parValues)  %#ok<INUSL>
% GENERATESENSITIVITYPARAMETERSET generates parameter for sensitivity calulation
%
% [parPathsNew,parValuesNew] = generateSensitivityParameterSet(WSettings,sensParameterList,parPaths,parValues)
%
% Inputs:
%       WSettings (structure)    definition of properties used in all
%                   workflow functions see GETDEFAULTWORKFLOWSETTINGS
%       sensParameterList  (cellarry) first column pathid of parameter,
%                       second column number of steps
%                       third column variation range
%                       5. column default value from xml
%   parPaths (cellarray of strings)  Paths of population parameter
%   parValues (double matrix)  Values of population parameter for one
%   individual
%
% Outputs;
%   parPathsNew (cellarray of strings)  Paths of population parameter with
%           added sensitivity parameter
%   parValuesNew (double matrix)  Values of population parameter with added
%           sensitivity values
%   SensPointer (structure) with one entry for each sensitivity parameter follwoing fields
%              indxInd (double vector) get first and last individual for this parameter
%              relValuesPar ( double vector) relValues of Parameter


% Open Systems Pharmacology Suite;  http://open-systems-pharmacology.org

% combine sensitivity parameters with popualtion parameter

[jj,ij] = ismember(sensParameterList(:,1),parPaths);
parPathsNew = parPaths;
% find new ones
ixPar = length(parPaths) + 1:sum(~jj);
parPathsNew(ixPar) = sensParameterList(~jj,1);
parValues(1,ixPar) = cell2mat(sensParameterList(~jj,5));

% combine old and new
ix = nan(size(jj));
ix(jj) = ij(jj);
ix(~jj) = ixPar;

% initialize new Values
parValuesNew = nan(2.*sum(cell2mat(sensParameterList(:,2))),length(parPathsNew));
k=0;

SensPointer= struct('indxInd',[],'relValuesPar',[]);

for iPar=1:size(sensParameterList,1)

    % get new parameter values
    nSteps = sensParameterList{iPar,2};
    vRange = sensParameterList{iPar,3};
    
    relValues =[1./(1+vRange* [1:nSteps]./nSteps) (1+vRange*[1:nSteps]./nSteps)]; %#ok<NBRAK>
    
    newValues = parValues(1,ix(iPar)).*relValues;
    
    % check if they are within allowed ranges
    jj = true(1,length(newValues));
    if ~isnan(sensParameterList{iPar,6})
        jj = jj & newValues >= sensParameterList{iPar,6};
    end
    if ~isnan(sensParameterList{iPar,7})
        jj = jj & newValues <= sensParameterList{iPar,7};
    end
    if any(jj)
        
        nInd = sum(jj);
        
        parValuesNew(k + [1:nInd],:) = repmat(parValues,nInd,1); %#ok<NBRAK>
        parValuesNew(k + [1:nInd],ix(iPar)) = newValues(jj); %#ok<NBRAK>
        
        SensPointer(iPar).indxInd = [k+1,k+nInd];
        SensPointer(iPar).relValuesPar = relValues(jj);
        k = k + nInd;
    else
         writeToReportLog('WARNING',sprintf('Parameter %s is outside the allowed definition range', sensParameterList{iPar,1}),false);
    end
end

% get rid of nan Columns
parValuesNew = parValuesNew(1:k,:);

return



