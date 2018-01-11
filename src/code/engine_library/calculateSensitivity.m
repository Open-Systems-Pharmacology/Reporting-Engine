function sens = calculateSensitivity(WSettings,PKPList,PKPListSens,SensPointer,ixBaseInd) %#ok<INUSL>
% CALCULATESENSITIVITY calculate sensitivity for PK Parameter
%
% sens = calculateSensitivity(WSettings,PKPList,PKPListSens,SensPointer,ixBaseInd)
%
% Inputs:
%       WSettings (structure)    definition of properties used in all
%                   workflow functions see GETDEFAULTWORKFLOWSETTINGS
%   PKPList (cellarray) PK Parameter of base individual,each cell correpsonds to one output, contains PK Parameters structures
%   PKPListSens (cellarray) PKParameter of sensitivity indiivduals,each cell correpsonds to one output, contains PK Parameters structures
%   SensPointer (structure) with one entry for each sensitivity parameter follwoing fields
%              indxInd (double vector) get first and last individual for this parameter
%              relValuesPar ( double vector) relValues of Parameter
% 
% Outputs:
%       sens (cellarray of structure) {iOutput,iPKParameter}, structure entries correpondend to sensitivity parameter
%               structure has following fields
%               slope (double) sensitivity
%               slopeCILower (double) lower confidence limit of 95%
%               slopeCIUpper (double) upper confidence limit of 95%
%               rSquare (double) R-square statistic
%               pValue (double) p Value of full mofel


% Open Systems Pharmacology Suite;  http://open-systems-pharmacology.org

% calculate sensitivity an write values to table array
for iO=1:length(PKPList)
    for iPK=1:length(PKPList{iO})
        
        for iP = 1:length(SensPointer)
        
            relValuesPar =  [SensPointer(iP).relValuesPar];
            relValuesPK = PKPListSens{iO}(iPK).value(SensPointer(iP).indxInd(1):SensPointer(iP).indxInd(end))./PKPList{iO}(iPK).value(ixBaseInd);
            
            
            [b,bint,~,~,stats]  = regress([relValuesPK,1]',[[relValuesPar,1]',ones(length(relValuesPar)+1,1)]);
            % calculate sesnitivity
            sens{iO,iPK}(iP).slope = b(1); %#ok<AGROW>
            sens{iO,iPK}(iP).slopeCILower = bint(1,1); %#ok<AGROW>
            sens{iO,iPK}(iP).slopeCIUpper = bint(1,2); %#ok<AGROW>
            sens{iO,iPK}(iP).rSquare = stats(1); %#ok<AGROW>
            sens{iO,iPK}(iP).pValue = stats(3); %#ok<AGROW>
            sens{iO,iPK}(iP).values = relValuesPK;  %#ok<AGROW>
        end        
    end
end

return
