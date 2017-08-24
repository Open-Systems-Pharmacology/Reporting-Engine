function [R] = getVPCRange(y,sensitivities,residuals,Jacob)
% GETVPCRANGE calcualtes ranges for VPC plots
%
% getVPCRange(y,sensitivities,residuals)
%
% Inputs:
%   y  (double matrix nT x 1)  property which should be plotted, nT = number of Timesteps
%   sensitivities  (double matrix nT x nP)  property which should be plotted, nP = number of optimized parameters
%   residuals  (double matrix nRes x 2)  first row resiudals, second row corresponding index of time
%   Jacob (doubple matric nP X nRes) (Jacobian) sesnitivity at the place of the residuals
%
% Outputs 
%   R structure with following fields
%   yMin (double vector 3 x m)  lower limit of range
%   yMax (double vector 3 x m)  max line of range
%   legendTextMean (string) description how the mean is calculated
%   rangeTxt (string) description how upper and lower limits are calculated
%   csvHeader (cellarray of string) dscription of all outputs {min,mean,max}




% degrees of freedom
ndgf = (size(Jacob,1)-size(Jacob,2));

% covariancematrix
dpdpMat = pinv(Jacob) * eye(length(residuals))*sum(residuals(1,:).^2)/ndgf * pinv(Jacob)';

% students t
studentsT90 = tinv(0.05,ndgf);



VPCflag = {'parameterUncertainty','dataUncertainty','parameterAndDataUncertainty'};

for iFlag = 1:length(VPCflag)
    % differnec to mean line
    diff = nan(size(y));
    switch VPCflag{iFlag}
        case 'parameterUncertainty'
            for iT = 1:length(y)
                diff(iT) = studentsT90.*sqrt(sensitivities(iT,:)*dpdpMat*sensitivities(iT,:)');
            end
        case 'dataUncertainty'
            for iT = 1:length(y)
                diff(iT) = studentsT90.*sqrt(sum(residuals(1,:).^2)/ndgf);
            end
        case 'parameterAndDataUncertainty'
            for iT = 1:length(y)
                diff(iT) = studentsT90.*sqrt(sensitivities(iT,:)*dpdpMat*sensitivities(iT,:)' + sum(residuals(1,:).^2)/ndgf);
            end
            
        otherwise
            
            error ('unknown VPCFlag %s,VPCflag');
    end
    
    % get line and rnage limits
    R.yMin(:,iFlag) = y - diff;
    R.yMax(:,iFlag) = y + diff;

end


% % get texts
R.legendTextMean = 'simulated timeprofile';

% get text for range legend
R.rangeTxt = sprintf('range of uncertainty');

R.csvHeader = {'Lower limit of uncertatinty',...
    R.legendTextMean,'Upper limit of uncertainty'};

return