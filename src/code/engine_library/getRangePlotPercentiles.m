function [yMin,yMean,yMax,legendTextMean,rangeTxt,csvHeader] = getRangePlotPercentiles(WSettings,y)
% GETRANGEPLOTPERCENTILES calculates the upper mean and lower line for ranges plots
%   dependening on the workflow settings
%
% [yMin,yMean,yMax,legendTextMean,rangeTxt,csvHeader] = getRangePlotPercentiles(WSettings,y)
%
% Inputs:
%   y  (double matrix nInd x m)  property which should be plotted 
% Outputs
%   yMin (double vector 1 x m)  lower limit of range
%   yMean (double vector 1 x m)  mid line of range
%   yMax (double vector 1 x m)  max line of range
%   legendTextMean (string) description how the mean is calculated
%   rangeTxt (string) description how upper and lower limits are calculated
%   csvHeader (cellarray of string) dscription of all outputs {min,mean,max}

if size(y,1) == 2 % only one individual
    yMin = nan;
    yMax = nan;
    yMean = y;
    legendTextMean = '';
    
    csvHeader = {};
else
    % exclude non processable individuals
    jjNonan = all(~isnan(y),2);
    yMin=prctile(y(jjNonan,:), WSettings.displayPercentiles(1),1);
    yMax=prctile(y(jjNonan,:), WSettings.displayPercentiles(end),1);
    switch WSettings.shadedAreaMeanType
        case 'median'
            yMean=median(y,1);
            legendTextMean='median';
        case 'geomean'
            yMean=geomean(y,1);
            legendTextMean='geo. mean';
        case 'mean'
            yMean=mean(y,1);
            legendTextMean='arith. mean';
    end

    % get text for range legend
    rangeTxt = sprintf('%s-%s percentile',...
        getPercentilePotenzText(WSettings.displayPercentiles(1)),...
        getPercentilePotenzText(WSettings.displayPercentiles(end)));

    csvHeader = {sprintf('%s percentile',getPercentilePotenzText(WSettings.displayPercentiles(1))),...
        legendTextMean,sprintf('%s percentile',getPercentilePotenzText(WSettings.displayPercentiles(end)))};
        
end

return