 function successInputCheck = checkInputDatafiles(WSettings,Datafiles)
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
     writeToLog(sprintf(['ERROR: variable Datafiles must be a cellarray size n x 3, \n',...
         'first entry datafile, second entry dictionary, third entry type of data']),WSettings.logfile,true,false);
     successInputCheck = false;
 else
     for iData = 1:size(Datafiles,1)
         if ~exist(Datafiles{iData,1},'file')
             writeToLog(sprintf('ERROR: Datafile %s does not exist',Datafiles{iData,1}),WSettings.logfile,true,false);
             successInputCheck = false;
         end
         if ~exist(Datafiles{iData,2},'file')
             writeToLog(sprintf('ERROR: Dictionary %s for timeprofile datafile does not exist',Datafiles{iData,2}),WSettings.logfile,true,false);
             successInputCheck = false;
         end
         if ~ismember(Datafiles{iData,3},{'timeprofile'})
             writeToLog(sprintf('ERROR: Datatype %s is unknown.',Datafiles{iData,3}),WSettings.logfile,true,false);
             successInputCheck = false;
         end
     end
 end