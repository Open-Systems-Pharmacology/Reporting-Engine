function workflowModeAsText = workflowModeToText(workflowType,workflowMode)
% WORKFLOWMODETOTEXT  generates text for logfile out of vartableworkflowMode
%
%   workflowModeAsText = workflowModeToText(workflowMode)
%
%   Inputs
%       workflowMode (string) chararcetrises type of workflow
%       workflowMode (string) chararcetrises mode of workflow
%

% Open Systems Pharmacology Suite;  http://open-systems-pharmacology.org


switch workflowType
    case 'popModel'
        if ismember(workflowMode,{'paralellComparison','pediatric', 'ratioComparison'})
            workflowModeAsText = sprintf('Population workflow, mode %s',workflowMode);
        else
            workflowModeAsText = sprintf('Population workflow, with project specific mode %s',workflowMode);
        end
    case 'meanModel'
        if ismember(workflowMode,{'default'})
            workflowModeAsText = sprintf('MeanModel workflow, mode %s',workflowMode);
        else
            workflowModeAsText = sprintf('MeanModel workflow, with project specific mode %s',workflowMode);
        end        
    otherwise
        workflowModeAsText = sprintf('Project specified workflowTye %s mode %s',workflowType,workflowMode);
end
   
return
        
