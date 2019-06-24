function exportSensitivityForPopulation(WSettings,analysisName,sens,sensParameterList,individualVector,prctileSelection,simulationName) 
% EXPORTSENSITIVITYFORPOPULATION export list of all caclualted sensitivities
%
% exportSensitivityForPopulation(WSettings,analysisName,sens,sensParameterList,individualVector,prctileSelection,simulationName)
%
% Inputs:
%       WSettings (structure)    definition of properties used in all
%                   workflow functions see GETDEFAULTWORKFLOWSETTINGS
%       analysisName (string) name used to creat file name
%       sens (cellarray of structures) with sensitivity informtaion
%               cell array correpsonds to OutputList, structur entries to
%               PK Parameter list
%       OutputList (structure) with output properties
%       sensParameterList  (cellarry) 1. column pathid of parameter,
%                       2. number of steps
%                       3. variation range
%                       4. minimal value
%                       5. maximal value
%                       6. column default value from xml
%       iPKList (double vector) index of PK Parameter

% Open Systems Pharmacology Suite;  http://open-systems-pharmacology.org



% initialize export structure
data(1,:) = {'simulationName','Individual','Parameter','sensitivity value','lower CI 95','upper CI 95','R-square','pValue'};

% fill structure
for iPop = 1:size(sens,2)
    
    for iInd =  1: size(sens,1)
        
        individualname = sprintf('%s-Individual ID %d',getPercentilePotenzText(prctileSelection(iInd)),individualVector(iInd));
        
        if ~isempty(sens{iInd,iPop})
        
            offset = size(data,1);
            nPar = size(sensParameterList,1);
            
            data(1,offset + (1:nPar)) = simulationName(iPop);
            data(2,offset + (1:nPar)) = {individualname};
            data(3,offset + (1:nPar)) = num2cell(sensParameterList(:,4));
            data(4,offset + (1:nPar)) = num2cell([sens{iInd,iPop}.slope]);
            data(5,offset + (1:nPar)) = num2cell([sens{iInd,iPop}.slopeCILower]);
            data(6,offset + (1:nPar)) = num2cell([sens{iInd,iPop}.slopeCIUpper]);
            data(7,offset + (1:nPar)) = num2cell([sens{iInd,iPop}.rSquare]);
            data(8,offset + (1:nPar)) = num2cell([sens{iInd,iPop}.pValue]);
        end
        
    end
    
end

% write file
fname = fullfile(WSettings.figures,'sensitivity',[analysisName,'.csv']);
writeTabCellArray(data,fname)

return
