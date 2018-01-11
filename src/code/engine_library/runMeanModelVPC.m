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

try
    writeToReportLog('INFO',sprintf('Start Mean Model VPC'),false);
    
    % time profiles
    % Initialize figureDir
    FP = ReportFigurePrint(fullfile(WSettings.figures,'timeprofile'),WSettings.printFormatList);
    
    for iD = 1:length(VPC.Timeprofile)
        FP = plotLoopVPCTimeprofiles(WSettings,VPC.textFunctionHandle,VPC.Timeprofile(iD),MeanModelSet,FP);
    end
    
    
    writeToReportLog('INFO',sprintf('Finalized Mean Model VPC \n'),false);
catch exception
        
    save(sprintf('exception_%s.mat',datestr(now,'ddmmyy_hhMM')),'exception');
    writeToReportLog('ERROR',exception.message,false);
    writeToReportLog('INFO',sprintf('Mean Model  VPC  finished with error \n'),false);
        
end

return


