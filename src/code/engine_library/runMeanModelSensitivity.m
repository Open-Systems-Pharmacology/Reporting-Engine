function runMeanModelSensitivity(WSettings,MeanModelSet)
% RUNMEANMODELSENSITIVITY runs the sensitivity analysis for all mean models defined in set
%
% runMeanModelSensitivity(WSettings,MeanModelSet)
%
% Inputs
%       WSettings (structure)    definition of properties used in all
%                   workflow functions see GETDEFAULTWORKFLOWSETTINGS
%       MeanModelSet (structure)   list of population simulations see GENERATEWORKFLOWINPUTFORMEANMODELSIMULATION


% Open Systems Pharmacology Suite;  http://open-systems-pharmacology.org

try
    writeToReportLog('INFO','start sensitivity analysis',false);
    
    % create resultdirector
    % Initialize figureDir
    FP = ReportFigurePrint(fullfile('figures','sensitivity'),WSettings.printFormatList);
    
    for iSet = 1:length(MeanModelSet)

        % get directory of temporyr varaibles
        tmpDir = fullfile('tmp',MeanModelSet(iSet).name);
        
        % load OutputList of population
        load(fullfile(tmpDir,'outputList.mat'),'OutputList');
        load(fullfile(tmpDir,'sensitivity.mat'),'sensParameterList');
        load(fullfile(tmpDir,'applicationProtocol.mat'),'PKParameterTemplate');

        % in restart mode use previously produced arrays
        sensFilename = fullfile(tmpDir,sprintf('sens_%d.mat',iSet));
        if WSettings.restart && exist(sensFilename,'file')
                load(sensFilename,'sens');
        else        
        
            % generate new "individuals" which differ by one parameter
            [parPaths,parValues,SensPointer] = generateSensitivityParameterSet(WSettings,sensParameterList,{},[]);
            nInd = size(parValues,1);
            
            % get result
            SimResult = generateSimResult(WSettings,MeanModelSet(iSet),nan,parPaths,parValues,{OutputList.pathID},nInd,zeros(1,nInd)); %#ok<NASGU>
            
            save(fullfile(tmpDir,'sens_simResult_1.mat'),'SimResult');
            
            % start the PKParameter calculation
            parPathsDM = {'Organism|Weight','Organism|Height'};
            parValuesDM(1,1) = getParameter('*Organism|Weight',1,'parametertype','readonly');
            parValuesDM(1,2) = getParameter('*Organism|Height',1,'parametertype','readonly');
            PKPListSens = calculatesPKParameterList(WSettings,MeanModelSet.name,MeanModelSet.calculatePKParameterFh,parPathsDM,parValuesDM,'sens_');
            
            save(fullfile(tmpDir,'sens_pKPList.mat'),'PKPListSens','SensPointer');
            load(fullfile(tmpDir,'pKPList.mat'),'PKPList');
            
            % calculate sesnitivity
            sens = calculateSensitivity(WSettings,PKPList,PKPListSens,SensPointer,1);
            save(sensFilename,'sens');        
        end
        
        % export results to csv
        load(fullfile('tmp',MeanModelSet.name,'outputList.mat'),'OutputList');
        exportSensitivity(WSettings,MeanModelSet.name,sens,OutputList,sensParameterList);
        
        
        header = sprintf('sensitivity analysis of %s ',MeanModelSet.reportName);
        FP = FP.iniCaptiontextFigtextArray(header,MeanModelSet.name);
        
        % create sensitivity Plots
        FP = plotSensitivity(WSettings,MeanModelSet.name,sens,OutputList,sensParameterList,FP,PKParameterTemplate);
        
        FP.saveCaptiontextArray;
    end
    
    writeToReportLog('INFO','finalize sensitivity analysis',false);

catch exception
        
    save(sprintf('exception_%s.mat',datestr(now,'ddmmyy_hhMM')),'exception');
    writeToReportLog('ERROR',exception.message,false);
    writeToReportLog('INFO',sprintf('Absorption plots finished with error \n'),false);
        
end       
    
return

function FP = plotSensitivity(WSettings,analysisName,sens,OutputList,sensParameterList,FP,PKParameterTemplate)


% fill structure
for iO = 1:length(OutputList)
    
    for iPK =   1:size(OutputList(iO).pKParameterList,2)
        [~,sortSumSensIx,iCut] = getListOfBestSensitivities(WSettings,sens(iO,iPK));
        
        FP = plotSensListMostSensitive(WSettings,FP,sens(iO,iPK),sortSumSensIx,iCut,sensParameterList,'o',[0 0 1],{},{});
        
        % get name and figure description
        jj = strcmp(OutputList(iO).pKParameterList{1,iPK},{PKParameterTemplate.name});
        PKReportName = PKParameterTemplate(jj).reportName;

        figureName = sprintf('%s_%s_%s_List',analysisName,PKReportName,OutputList(iO).reportName);
        figtxt = sprintf('Most sensitiv parameter for %s for %s', PKReportName,OutputList(iO).reportName);
        % save figure
        FP = FP.printFigure(figureName,figtxt);

        
    end
end

return


function exportSensitivity(WSettings,analysisName,sens,OutputList,sensParameterList) %#ok<INUSL>
% EXPORTSENSITIVITY export list of all caclualted sensitivities
%
% exportSensitivity(WSettings,analysisName,sens,OutputList,sensParameterList)
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

% Open Systems Pharmacology Suite;  http://open-systems-pharmacology.org



% initialize export structure
data(1,:) = {'Output','PK-Parameter','Parameter','sensitivity value','lower CI 95','upper CI 95','R-square','pValue'};
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
for iO = 1:length(OutputList)
    
    for iPK =   1:size(OutputList(iO).pKParameterList,2)
        
        offset = length(data{3,1});
        nPar = size(sensParameterList,1);
        
        data{3,1}(offset+[1:nPar],1) = {OutputList(iO).reportName};
        data{3,2}(offset+[1:nPar],1) = OutputList(iO).pKParameterList(1,iPK);
        data{3,3}(offset+[1:nPar],1) = sensParameterList(:,1);
        data{3,4}(offset+[1:nPar],1) = [sens{iO,iPK}.slope];
        data{3,5}(offset+[1:nPar],1) = [sens{iO,iPK}.slopeCILower];
        data{3,6}(offset+[1:nPar],1) = [sens{iO,iPK}.slopeCIUpper];
        data{3,7}(offset+[1:nPar],1) = [sens{iO,iPK}.rSquare];
        data{3,8}(offset+[1:nPar],1) = [sens{iO,iPK}.pValue];            
        
    end
    
end

% write file
if ~exist(fullfile('figures','sensitivity'),'dir')
    mkdir(fullfile('figures','sensitivity'));
end
fname = fullfile('figures','sensitivity',[analysisName,'.csv']);

writetab(fname, data, ';', 0, 0, 1,0);
return