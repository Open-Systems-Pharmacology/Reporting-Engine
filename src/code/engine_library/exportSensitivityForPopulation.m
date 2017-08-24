function exportSensitivityForPopulation(WSettings,analysisName,sens,sensParameterList,individualVector,prctileSelection,simulationName) %#ok<INUSL>
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
data(2,:) = {'string','string','string','double','double','double','double','double'};
for iCol = 1:size(data,2)
    switch data{2,iCol}
        case 'string'
            data{3,iCol} = {};
        case 'double'
            data{3,iCol} = [];
    end
end

% fill structure
for iPop = 1:size(sens,2)
    
    for iInd =  1: size(sens,1)
        
        individualname = sprintf('%s-Individual ID %d',getPercentilePotenzText(prctileSelection(iInd)),individualVector(iInd));
        
        if ~isempty(sens{iInd,iPop})
        
            offset = length(data{3,1});
            nPar = size(sensParameterList,1);
            
            data{3,1}(offset+[1:nPar],1) = simulationName(iPop);
            data{3,2}(offset+[1:nPar],1) = {individualname};
            data{3,3}(offset+[1:nPar],1) = sensParameterList(:,4);
            data{3,4}(offset+[1:nPar],1) = [sens{iInd,iPop}.slope];
            data{3,5}(offset+[1:nPar],1) = [sens{iInd,iPop}.slopeCILower];
            data{3,6}(offset+[1:nPar],1) = [sens{iInd,iPop}.slopeCIUpper];
            data{3,7}(offset+[1:nPar],1) = [sens{iInd,iPop}.rSquare];
            data{3,8}(offset+[1:nPar],1) = [sens{iInd,iPop}.pValue];
        end
        
    end
    
end

% write file
fname = fullfile('figures','sensitivity',[analysisName,'.csv']);

writetab(fname, data, ';', 0, 0, 1,0);
return