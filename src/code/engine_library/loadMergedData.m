function [Data,TPnew,Dict] = loadMergedData(WSettings,listOfname) %#ok<INUSL>
% LOADMERGEDDATA load data for different simualtion and merge them
% 
% Data = loadMergedData(WSettings,{PopRunSet(ix).name});
% 
% Inputs:
%       WSettings (structure)    definition of properties used in all
%                   workflow functions see GETDEFAULTWORKFLOWSETTINGS
%       listOfname (cell array of strings)   list dfines the selection of  populations
%       parPathSelection (cell array of strings 2 x n) first row list the selected paths
%                               second row corresponding display units

% Open Systems Pharmacology Suite;  http://open-systems-pharmacology.org


% initialize output
Data = [];
TPnew = [];
Dict = [];

% loop on simulations
for iName = 1:length(listOfname)

    fname = fullfile('tmp',listOfname{iName},'dataTP.mat');
    if exist(fname,'file')
        load(fname);
    
        if isempty(Data)
            Data = DataTP;
            TPnew = TP;  %#ok<USENS>
        else
            % get data entries, already existing
            jj = [Data.sid] == [DataTP.sid] & [Data.stud] == [DataTP.stud];
            % add the other ones
            Data = [Data,DataTP(~jj)]; %#ok<AGROW>
            for iO = 1:length(TP)
                TPnew{iO} = [TPnew{iO},TP{iO}(~jj)];  %#ok<AGROW>
            end
        end
    end
end


return
        
