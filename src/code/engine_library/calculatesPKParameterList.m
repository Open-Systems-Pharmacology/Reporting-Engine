function [PKPList,individualIdVector,pathList] = calculatesPKParameterList(WSettings,simulationName,calculatePKParameterFh,parPaths,parValues,simResultPrefix)
%CALCULATESPKPARAMETERLIST calculates the Pkparameter
%
% [PKPList,individualIdVector,pathList] = calculatesPKParameterList(WSettings,simulationName,calculatePKParameterFh,parPaths,parValues)
%
% Inputs:
%   WSettings   (structure) containing global settings see GETDEFAULTWORKFLOWSETTINGS
%   simulationName      (string) name of poulation or meanmodel set
%   calculatePKParameterFh  ( functionhandle) function to calculate PK
%           Parameter, if empty default function is taken
%   parPaths (cellarray of strings) pathnames of parameter
%   parValues ( double array) values of parameter
%   simResultPrefix (string) if simResult is not defualt result, e.g. during sensitivity
%           this is added to matfile, if empty it is stnadard simulationname
%
% Outputs:
%   PKPList (cellarray) each cell correpsonds to one output, contains PK Parameters structures
%   individualIdVector (double) vector of individuals ids
%   pathList (cellarray of string) outputpathes

% Open Systems Pharmacology Suite;  http://open-systems-pharmacology.org


if ~exist('simResultPrefix','var')
    simResultPrefix = '';
end


% load applicationProtocol
load(fullfile('tmp',simulationName,'applicationProtocol.mat'),'ApplicationProtocol','isValid','PKParameterTemplate');
load(fullfile('tmp',simulationName,'outputList.mat'),'OutputList');
pathList = {OutputList.pathID}; %#ok<NODEF>

% initialize PK List
PKPList = cell(length(OutputList) ,1);

if all(cellfun(@isempty,{OutputList(:).pKParameterList}))
    % no PKParameter defined, so nothing to do
    return
end

% decide on the PKParameter function and calculate PKParameter
% if given take the projectspecific function
if ~isempty(calculatePKParameterFh)
    
    for iO = 1:length(OutputList)
        if ~isempty(OutputList(iO).pKParameterList)
            
            [simTime,simValues,~,~,individualIdVector]  = loadSimResult(simulationName,iO,'',simResultPrefix);
            
            tmp = feval(calculatePKParameterFh,WSettings,ApplicationProtocol,simTime,simValues,parPaths,parValues);
            [jj,ix] = ismember(OutputList(iO).pKParameterList{1,iPK},{PKParameterTemplate.name});
            PKPList{iO} = tmp(ix(jj));
        end
        
    end
    
elseif isValid
    
    jj = strcmp('Organism|Weight',parPaths);
    weight =  parValues(:,jj);
    jj = strcmp('Organism|Height',parPaths);
    height =  parValues(:,jj);
    
    for iO = 1:length(OutputList)

        if ~isempty(OutputList(iO).pKParameterList)
            
             [simTime,simValues,~,~,individualIdVector]  = loadSimResult(simulationName,iO,'',simResultPrefix);

            
            tmp = calculatePKParameterForApplicationProtocol(WSettings,ApplicationProtocol,simTime,simValues,weight,height);

            [jj,ix] = ismember(OutputList(iO).pKParameterList(1,:),{PKParameterTemplate.name});
            PKPList{iO} = tmp(ix(jj));
        else
            PKPList{iO} = [];
        end
        
    end
else
    
    error('application protocol could not be interpreted. Please use project specific function for PK Parameter calculation');
    
end
