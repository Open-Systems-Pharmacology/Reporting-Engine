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
    FP = ReportFigurePrint(fullfile(WSettings.figures,'sensitivity'),WSettings.printFormatList);
    
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
            PKPListSens = calculatesPKParameterList(WSettings,MeanModelSet(iSet).name,MeanModelSet(iSet).calculatePKParameterFh,parPathsDM,parValuesDM,'sens_');
            
            save(fullfile(tmpDir,'sens_pKPList.mat'),'PKPListSens','SensPointer');
            load(fullfile(tmpDir,'pKPList.mat'),'PKPList');
            
            % calculate sesnitivity
            sens = calculateSensitivity(WSettings,PKPList,PKPListSens,SensPointer,1);
            save(sensFilename,'sens');        
        end
        
        % export results to csv
        load(fullfile('tmp',MeanModelSet(iSet).name,'outputList.mat'),'OutputList');
        exportSensitivity(WSettings,MeanModelSet(iSet).name,sens,OutputList,sensParameterList);
        
        
        header = sprintf('Sensitivity analysis of %s ',MeanModelSet(iSet).reportName);
        FP = FP.iniCaptiontextFigtextArray(header,MeanModelSet(iSet).name);
        
        % create sensitivity Plots
        FP = plotSensitivity(WSettings,MeanModelSet(iSet).name,sens,OutputList,sensParameterList,FP,PKParameterTemplate);
        
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

        figureName = removeForbiddenLetters(sprintf('%s_%s_%s_List',analysisName,PKParameterTemplate(jj).name,OutputList(iO).reportName));
        figtxt = sprintf('Most sensitive parameter for %s for %s', PKReportName,OutputList(iO).reportName);
        % save figure
        FP = FP.printFigure(figureName,figtxt);

        
    end
end

return


function exportSensitivity(WSettings,analysisName,sens,OutputList,sensParameterList) 
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

% fill structure
for iO = 1:length(OutputList)
    
    for iPK =   1:size(OutputList(iO).pKParameterList,2)
        
        offset = size(data,1);
        nPar = size(sensParameterList,1);
        
        data(1,offset + (1:nPar)) = {OutputList(iO).reportName};
        data(2,offset + (1:nPar)) = {OutputList(iO).pKParameterList(1,iPK)};
        data(3,offset + (1:nPar)) = num2cell(sensParameterList(:,4));
        data(4,offset + (1:nPar)) = num2cell([sens{iO,iPK}.slope]);
        data(5,offset + (1:nPar)) = num2cell([sens{iO,iPK}.slopeCILower]);
        data(6,offset + (1:nPar)) = num2cell([sens{iO,iPK}.slopeCIUpper]);
        data(7,offset + (1:nPar)) = num2cell([sens{iO,iPK}.rSquare]);
        data(8,offset + (1:nPar)) = num2cell([sens{iO,iPK}.pValue]);
        
    end
    
end

% write file
if ~exist(fullfile(cd,WSettings.figures,'sensitivity'),'dir')
    mkdir(fullfile(WSettings.figures,'sensitivity'));
end
fname = fullfile(WSettings.figures,'sensitivity',[analysisName,'.csv']);

writeTabCellArray(data,fname);

return
