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
%       PopRunSet (structure)   list of population simulations see GETDEFAULTPOPRUNSET
%       FP (ReportFigurPrint)  objects with manages print of figures and
%                           tables
%   Outputs:
%       FP (ReportFigurPrint)  objects with manages print of figures and
%       tables


% Open Systems Pharmacology Suite;  http://forum.open-systems-pharmacology.org
% Date: 27-July-2017

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
end

% start Loop
for iSet = 1:nOuterLoop
        
    %initialize new sheet
    header = sprintf('Physiology of %s',populationReportName{iSet});
    FP = FP.iniCaptiontextFigtextArray(header,sheetNameList{iSet});
    
    % merge population if necessayr
    ix = Def.ixOfPopRunSets(iSet:max(nInnerLoop,iSet));
    parValues = loadSelectedPathsOfPopulation(WSettings,{PopRunSet(ix).name},parPathSelection);

    for iY = 1:size(Def.yList,2)

        y = parValues(:,iY);
        if ~isempty(Def.ixPopRunSetRef)
            yRef = parValuesRef(:,iY);
        else
            yRef = [];
        end
        
        % if no property for x is give do a histogramm
        if isempty(Def.xList)
            
            % get name and figure description
            figureName = sprintf('h%d_%s',iY,removeForbiddenLetters(Def.yList{iY,2}));
            [figtxt,figtxtTable,legendEntries] = feval(textFunctionHandle,WSettings,'physHist',...
                {Def.yList{iY,2},populationReportName{iSet},refPopulationReportName},...
                [length(y),length(yRef)]);
            
            
            % do figure
            csv = plotReportHistogram(WSettings,FP.figureHandle,y,yRef,Def.yList{iY,2},Def.yList{iY,3},...
                legendEntries);
            
            
            % save figure
            FP = FP.printFigure(figureName,figtxt,csv,figtxtTable);
            
                
        else
            
            for iX = 1:size(Def.xList,2)
                plotLayoutShadedArea(WSettings,x,y,[],yRef,Def.x_list{iY,2},Def.x_list{iY,3},...
                    Def.yList{iY,2},Def.yList{iY,3});
            end
            
        end
    end
    
    
    FP.saveCaptiontextArray;
    
end
    
return
