 function successInputCheck = checkInputDatafiles(WSettings,Datafiles) %#ok<INUSL>
 % CHECKINPUTDATAFILES check if varaibale Datafiles is correctly given
%
% successInputCheck = checkInputDatafiles(WSettings,Datafiles)
%
% Inputs:
%       WSettings (structure)    definition of properties used in all
%                   workflow functions see GETDEFAULTWORKFLOWSETTINGS
%       Datafiles (cellarray of strings) list of data files with dictionary and type of data 
%
% Outputs:
%       successInputCheck (boolean) ture, if varaible is correctly designed

% Open Systems Pharmacology Suite;  http://open-systems-pharmacology.org

successInputCheck = true;
 
 if ~iscell(Datafiles) || size(Datafiles,2)~=3
     writeToReportLog('ERROR',sprintf(['variable Datafiles must be a cellarray size n x 3, \n',...
         'first entry datafile, second entry dictionary, third entry type of data']),false);
     successInputCheck = false;
 else
     for iData = 1:size(Datafiles,1)
         if ~exist(Datafiles{iData,1},'file')
             writeToReportLog('ERROR',sprintf('Datafile %s does not exist',Datafiles{iData,1}),false);
             successInputCheck = false;
         end
         if ~exist(Datafiles{iData,2},'file')
             writeToReportLog('ERROR',sprintf('Dictionary %s for timeprofile datafile does not exist',Datafiles{iData,2}),false);
             successInputCheck = false;
         end
         if ~ismember(Datafiles{iData,3},{'timeprofile'})
             writeToReportLog('ERROR',sprintf('Datatype %s is unknown.',Datafiles{iData,3}),false);
             successInputCheck = false;
         end
     end
 end