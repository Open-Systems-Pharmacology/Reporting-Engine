function [allSensValues,sortSumSensIx,iCut,CsortSumSensPop] = getListOfBestSensitivities(WSettings,sens)
% GETLISTOFBESTSENSITIVITIES  get all sensitivity above cutoff value
%
% [allSensValues,sortSumSensIx,iCut,CsortSumSensPop] = getListOfBestSensitivities(WSettings,sensPop)
%
% Inputs:
%       WSettings (structure)    definition of properties used in all
%                   workflow functions see GETDEFAULTWORKFLOWSETTINGS
%       sens (cellarray of structures) with sensitivity informtaion
%               cell array correpsonds to OutputList, structur entries to
%               PK Parameter list
%
% Outputs:
%   allSensValues (double matrix nPar x nPop*nPrc,3) information of cell is vonverted to matrix for easy plot
%   sortSumSensIx (doubleVector) index vector of sorted sensitivities
%   iCut (double) number of parameters above cut
%   CsortSumSensPop (structure)  cumulated sesnitivity


% Open Systems Pharmacology Suite;  http://open-systems-pharmacology.org

nPop = size(sens,2);
nPrc = size(sens,1);
nPar = length(sens{1,1});

% initialize array 
allSensValues=nan(nPar,nPop*nPrc,3);

k=0;
for iPop=1:nPop
    for iPrc=1:nPrc
        k=k+1;
        allSensValues(:,k,1) = [sens{iPrc,iPop}.slope];
        allSensValues(:,k,2) = [sens{iPrc,iPop}.slopeCILower];
        allSensValues(:,k,3) = [sens{iPrc,iPop}.slopeCIUpper];
    end
end

% sort senitivity
[~,sortSumSensIx]=sort(max(abs(allSensValues(:,:,1)),'',2),'descend');
% [~,sortSumSens_ix]=sort(abs(median(allSensValues,2)),'descend');
% [~,sortSumSens_ix]=sort(abs(mean(allSensValues,2)),'descend');

% get cumlated sensitivities above cutoff
iCut=0;
CsortSumSensPop=struct;
for iPop=1:nPop
    for iPrc=1:nPrc      
        tmp=sqrt(cumsum([sens{iPrc,iPop}(sortSumSensIx(end:-1:1)).slope].^2));
        if ~isempty(tmp) &&  ~isnan( WSettings.sensitivityCutoff)
            tmp = tmp(end:-1:1)./tmp(end);
            jj=find(tmp < (1 - WSettings.sensitivityCutoff),1);
            iCut=max(jj-1,iCut);
        else
            iCut = length(tmp);
        end
        CsortSumSensPop(iPrc,iPop).slope=tmp;
    end
end
if isempty(iCut)
    iCut=1;
end
return

