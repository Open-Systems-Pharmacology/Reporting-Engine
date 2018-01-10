function [Data,TP,Dict,dataReportName,dataPopIndex] = loadMergedData(WSettings,RunSet,popRunSetMerger) %#ok<INUSL>
% LOADMERGEDDATA load data for different simualtion and merge them
% 
% Data = loadMergedData(WSettings,{PopRunSet(ix).name});
% 
% Inputs:
%       WSettings (structure)    definition of properties used in all
%                   workflow functions see GETDEFAULTWORKFLOWSETTINGS
%       listOfname (cell array of strings)   list defines the selection of
%           populations or mean models
% 
% Outputs:
%       Data (structure) contains demographic an physiological data per
%       individuals
%       TP (cellarray of structures) contains timeprofiles per individual,
%       eache cell entry corresponds to one output
%       Dict (structure) Dictionayr for data interpretation
%       dataReportName (name for legend and caption texts

% Open Systems Pharmacology Suite;  http://open-systems-pharmacology.org


% initialize output
Data = [];
TP = []; 
Dict = [];
dataPopIndex = [];

% get first Name
if ~exist('popRunSetMerger','var')
    popRunSetMerger = 'serial';
end
    
switch popRunSetMerger
    case 'merged'
        dataReportName = {RunSet.dataReportName};
    case {'serial','totalMerge'}
        dataReportName = {RunSet(1).dataReportName};
    otherwise
        error('unknonw merger');
end
% loop on simulations
for iName = 1:length(RunSet)

    fname = fullfile('tmp',RunSet(iName).name,'dataTP.mat');
    if exist(fname,'file')
        TP = {};
        load(fname);  %#ok<LOAD>
    
        if isempty(Data)
            Data = DataTP;
            TPnew = TP;  
            dataPopIndex = ones(size(Data))*iName;
        else
            % get data entries, already existing
            jj =  ismember([DataTP.sid],[Data.sid]) & ismember([DataTP.stud],[Data.stud]);
            % add the other ones
            Data = [Data,DataTP(~jj)]; %#ok<AGROW>
            dataPopIndex(end+1:length(Data)) = ones(1,sum(~jj))*iName;
            for iO = 1:length(TP)
                if ~isempty(TP{iO}) 
                    TPnew{iO} = [TPnew{iO},TP{iO}(~jj)];  
                end
            end
            
        end
    end
end
    
if isempty(Dict)
    return
end

% return data in MoBi Units (min) of Time
jj = strcmp({Dict.matlabID},'time');
dataTimeUnitFactor = getUnitFactor(Dict(jj).nonmemUnit,'min','Time');

TP = TPnew;
for iO = 1:length(TP)
    for iInd = 1:length( TP{iO})
        TP{iO}(iInd).time = TP{iO}(iInd).time.*dataTimeUnitFactor; 
    end
end


return
        
