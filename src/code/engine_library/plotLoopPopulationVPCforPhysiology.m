function     FP = plotLoopPopulationVPCforPhysiology(WSettings,textFunctionHandle,Def,PopRunSet,FP)
% PLOTLOOPPOPULATIONVPCFORPHYSIOLOGY  does physiologocal plots for population VPC
%
%   FP = plotLoopPopulationVPCforPhysiology(WSettings,Def,PopRunSet,FP)
% 
% Inputs:
%       WSettings (structure)    definition of properties used in all
%                   workflow functions see GETDEFAULTWORKFLOWSETTINGS
%       textFunctionHandle (function handle) function handle to set figure text
%       Def (structure) contains information which plots should be
%               generated  see GETDEFAULTVPCPOPULATIONSETTINGS
%               (timeprofile)
%       PopRunSet (structure)   list of population simulations see GENERATEWORKFLOWINPUTFORPOPULATIONSIMULATION
%       FP (ReportFigurPrint)  objects with manages print of figures and
%                           tables
%   Outputs:
%       FP (ReportFigurPrint)  objects with manages print of figures and
%       tables


% Open Systems Pharmacology Suite;  http://open-systems-pharmacology.org


% get selectionList of population parameter
if isempty(Def.xList)
    parPathSelection = [Def.yList(:,[1 3])];
else
    parPathSelection = [Def.yList(:,[1 3]); Def.xList(:,[1 3])];
end
    
% get Reference population
if ~isempty(Def.ixPopRunSetRef)
    parValuesRef = loadSelectedPathsOfPopulation(WSettings,{PopRunSet(Def.ixPopRunSetRef).name},parPathSelection);
    refPopulationReportName = PopRunSet(Def.ixPopRunSetRef).popReportName;
else
    parValuesRef=[];
    refPopulationReportName = '';
end

% get loopIndices depending on popRunSetMerger
switch Def.popRunSetMerger
    case 'merged'
        nOuterLoop = 1;
        nInnerLoop = length(Def.ixOfPopRunSets);
        sheetNameList = {Def.name};
        populationReportName = {Def.reportName};
    case 'serial'
        nOuterLoop = length(Def.ixOfPopRunSets);
        nInnerLoop = 1;
        sheetNameList = {PopRunSet(Def.ixOfPopRunSets).name};
        populationReportName = {PopRunSet(Def.ixOfPopRunSets).popReportName};
    otherwise
        error('flag');

end

% start Loop
for iSet = 1:nOuterLoop
        
    %initialize new sheet
    header = sprintf('Physiology of %s',populationReportName{iSet});
    FP = FP.iniCaptiontextFigtextArray(header,sheetNameList{iSet});
    
    % merge population if necessary
    ixSet = Def.ixOfPopRunSets(iSet:max(nInnerLoop,iSet));
    [parValues,sourceIndex] = loadSelectedPathsOfPopulation(WSettings,{PopRunSet(ixSet).name},parPathSelection);
    
    if size(parValues,1) < WSettings.rangePlotsMin
        warningflagRangePlots = true;
    else
        warningflagRangePlots = false;
    end

    % load data if available
    [Data,~,Dict] = loadMergedData(WSettings,{PopRunSet(ixSet).name});
    
    for iY = 1:size(Def.yList,1)

        jY = strcmp(Def.yList{iY,1},parPathSelection(:,1));
        y = parValues(:,jY);
        
        % reference
        if ~isempty(Def.ixPopRunSetRef)
            yRef = parValuesRef(:,iY);
        else
            yRef = [];
        end
        
        % data
        if ~isempty(Data)
            jjDict = strcmp({Dict.pathID},parPathSelection{iY});
            yData = [Data.(Dict(jjDict).matlabID)]';
        else
            yData = [];
        end
        dataReportName = PopRunSet(ixSet(1)).dataReportName;
        
        % if no property for x is given or if it a
        % categoric property like gender do a histogramm
        if isempty(Def.xList) || ~isempty(Def.yList{iY,4})
            
            
            % add ref as additional population
            y = [y;yRef]; %#ok<AGROW>
            sourceIndexWithRef = [sourceIndex; repmat(max(sourceIndex+1),length(yRef),1)];
            popLabels = [{PopRunSet(ixSet).boxwhiskerLabel},{PopRunSet(Def.ixPopRunSetRef).boxwhiskerLabel}];
            nPop = arrayfun(@(x) sum(sourceIndexWithRef==x),[1:max(sourceIndexWithRef)]);
            
            % get name and figure description
            figureName = sprintf('h%d_%s',iY,removeForbiddenLetters(Def.yList{iY,2}));
            [figtxt,figtxtTable,legendEntries] = feval(textFunctionHandle,WSettings,'physHist',...
                {Def.yList{iY,2},populationReportName{iSet},dataReportName,popLabels,nPop,length(yData)});
            
            % do figure
            csv = plotReportHistogram(WSettings,FP.figureHandle,y,sourceIndexWithRef,yData,Def.yList{iY,2},...
                Def.yList{iY,3},Def.yList{iY,4},legendEntries);
            
            
            % save figure
            FP = FP.printFigure(figureName,figtxt,csv,figtxtTable);
            
                
        else
            
            for iX = 1:size(Def.xList,1)
                
                jX = strcmp(Def.xList{iX,1},parPathSelection);
                x = parValues(:,jX);
        
            
                if warningflagRangePlots
                   writeToLog(sprintf(['WARNING: you are creating Rangeplots with %d individuals, recommended are at least %d individual. ',...
                       'Statistic might be not sufficient'],size(parValues,1),WSettings.rangePlotsMin),WSettings.logfile,true,false); 
                    warningflagRangePlots = false;
                end
                
                
                % get name and figure description
                figureName = sprintf('shA%d_%d_%s_%s',iY,iX,removeForbiddenLetters(Def.yList{iY,2}),...
                    removeForbiddenLetters(Def.xList{iX,2}));

                [figtxt,figtxtTable,legendEntries] = feval(textFunctionHandle,WSettings,'physShA',...
                    {Def.xList{iX,2},Def.yList{iY,2},populationReportName{iSet},refPopulationReportName,dataReportName,...
                    [length(y),length(yRef),length(yData)]});
                
                
                csv = plotReportShadedArea(WSettings,FP.figureHandle,x,y,yRef,Def.xList{iX,2},Def.xList{iX,3},...
                    Def.yList{iY,2},Def.yList{iY,3},'lin',legendEntries);
                
                % save figure
                FP = FP.printFigure(figureName,figtxt,csv,figtxtTable);

                
            end
            
        end
    end
    
    
    FP.saveCaptiontextArray;
    
end
    
return
