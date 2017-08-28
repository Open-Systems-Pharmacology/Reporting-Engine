function [R] = getVPCRange(y,sensitivities,residuals,Jacob,scale)
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

load test;

R = [];

% degrees of freedom
ndgf = (size(Jacob,1)-size(Jacob,2));

if ndgf <=0
    return;
end


% covariancematrix
%  covariance matrix = fisher_inv* Jacob'*eye(length(resid))*sum(resid.^2)/(size(Jacob,1)-size(Jacob,2))*Jacob*fisher_inv;  
% pinv(Jacob) = fisher_inv* Jacob;
M = eye(length(residuals))*sum(residuals.^2)/ndgf;
dpdpMat = pinv(Jacob)*M*pinv(Jacob)';

% students t
studentsT90 = tinv(0.975,ndgf);



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
                diff = studentsT90.*sqrt(sum(residuals.^2)/ndgf)*ones(size(y));
            end
        case 'parameterAndDataUncertainty'
            for iT = 1:length(y)
                diff(iT) = studentsT90.*sqrt(sensitivities(iT,:)*dpdpMat*sensitivities(iT,:)' + sum(residuals.^2)/ndgf);
            end
            
        otherwise
            
            error ('unknown VPCFlag %s,VPCflag');
    end
    
    % get line and rnage limits
    switch scale
        case 'lin'
            R.yMin(:,iFlag) = y-diff;
            R.yMax(:,iFlag) = y+diff;
        case 'log'
            R.yMin(:,iFlag) = y.*(1+diff);
            R.yMax(:,iFlag) = y./(1+diff);
        otherwise
            eror('unknown scale');
    end

end


% % get texts
R.legendTextMean = 'simulated timeprofile';

% get text for range legend
R.rangeTxt = sprintf('range of uncertainty');

R.csvHeader = {'Lower limit of uncertatinty',...
    R.legendTextMean,'Upper limit of uncertainty'};

return