function runPopulationVPC(WSettings,VPC,PopRunSet)
%RUNPOPULATIONVPC runs the visual predictive check for VPC
%
% runPopulationVPC(WSettings,VPC,PopRunSet)
%
% Inputs
%       WSettings (structure)    definition of properties used in all
%                   workflow functions see GETDEFAULTWORKFLOWSETTINGS
%       VPC (structure) contains information which plots should be
%               generated  see GETDEFAULTVPCPOPULATIONSETTINGS
%       PopRunSet (structure)   list of population simulations see GENERATEWORKFLOWINPUTFORPOPULATIONSIMULATION


% Open Systems Pharmacology Suite;  http://open-systems-pharmacology.org

try
    writeToReportLog('INFO',sprintf('Start Population VPC'),false);
    
    
    % demographic
    % loop on demographic definitions
    % Initialize figureDir
    FP = ReportFigurePrint(fullfile(WSettings.figures,'physiology'),WSettings.printFormatList);
    
    for iD = 1:length(VPC.PhysProperties)
        FP = plotLoopPopulationVPCforPhysiology(WSettings,VPC.textFunctionHandle,VPC.PhysProperties(iD),PopRunSet,FP);
    end
    
    % time profiles
    % Initialize figureDir
    FP = ReportFigurePrint(fullfile(WSettings.figures,'timeprofile'),WSettings.printFormatList);
    
    for iD = 1:length(VPC.Timeprofile)
        FP = plotLoopVPCTimeprofiles(WSettings,VPC(iD).textFunctionHandle,VPC(iD).Timeprofile(iD),PopRunSet,FP);
    end
    
    
    % pk parameter
    % Initialize figureDir
    FP = ReportFigurePrint(fullfile(WSettings.figures,'pKParameter'),WSettings.printFormatList);
    
    for iD = 1:length(VPC.PKParameter)
        FP = plotLoopPopulationVPCPKparameter(WSettings,VPC.textFunctionHandle,VPC.PKParameter(iD),PopRunSet,FP);
    end
    
    
    writeToReportLog('INFO',sprintf('Finalized Population VPC \n'),false);

catch exception
        
    save(sprintf('exception_%s.mat',datestr(now,'ddmmyy_hhMM')),'exception');
    writeToReportLog('ERROR',exception.message,false,exception);
    writeToReportLog('INFO',sprintf('Population VPC  finished with error \n'),false);
        
end
    
return


