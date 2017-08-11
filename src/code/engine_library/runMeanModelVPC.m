function runMeanModelVPC(WSettings,VPC,MeanModelSet)
%RUNMEANMODELVPC runs the visual predictive check for VPC
%
% runMeanModelVPC(WSettings,VPC,MeanModelSet)
%
% Inputs
%       WSettings (structure)    definition of properties used in all
%                   workflow functions see GETDEFAULTWORKFLOWSETTINGS
%       VPC (structure) contains information which plots should be
%               generated  see GETDEFAULTVPCMEANMODELSETTINGS
%       MeanModelSet (structure)   list of population simulations see GENERATEWORKFLOWINPUTFORMEANMODELSIMULATION


% Open Systems Pharmacology Suite;  http://open-systems-pharmacology.org

writeToLog(sprintf('Start Mean Model VPC'),WSettings.logfile,true,false);

% time profiles
% Initialize figureDir
FP = ReportFigurePrint(fullfile('figures','timeprofile'),WSettings.printFormatList);

for iD = 1:length(VPC.Timeprofile)
    FP = plotLoopVPCTimeprofiles(WSettings,VPC.textFunctionHandle,VPC.Timeprofile(iD),MeanModelSet,FP);
end


writeToLog(sprintf('Finalized Mean Model VPC \n'),WSettings.logfile,true,false);

return


