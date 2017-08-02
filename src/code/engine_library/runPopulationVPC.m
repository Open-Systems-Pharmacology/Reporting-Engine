function runPopulationVPC(Settings,VPC,PopRunSet)
%RUNPOPULATIONVPC runs the visual predictive check for VPC

% Open Systems Pharmacology Suite;  http://forum.open-systems-pharmacology.org
% Date: 26-July-2017


writeToLog(sprintf('Start Population VPC'),Settings.logfile,true,false);


% demographic
% loop on demographic definitions
% Initialize figureDir
FP = ReportFigurePrint(fullfile('figures','demographic'),Settings.printFormatList);

for iD = 1:length(VPC.demographic)
    FP = runDemographicVPC(Settings,VPC.demographic(iD),PopRunSet,FP);
end

writeToLog(sprintf('Finalized Population VPC \n'),Settings.logfile,true,false);

return


function     FP = runDemographicVPC(Settings,demographic,PopRunSet,FP)

% get selectionList of population parameter
if isempty(demographic.xList)
    parPathSelection = [demographic.yList(:,[1 3])];
else
    parPathSelection = [demographic.yList(:,[1 3]); demographic.xList(:,[1 3])];
end
    
% get Reference population
if ~isempty(demographic.ixPopRunSetRef)
    parValuesRef = loadPop(Settings,{PopRunSet(demographic.ixPopRunSetRef).name},parPathSelection);
    refPopulationReportName = PopRunSet(demographic.ixPopRunSetRef).popReportName;
else
    parValuesRef=[];
    refPopulationReportName = '';
end

% get loopIndices depending on popRunSetMerger
switch demographic.popRunSetMerger
    case 'merged'
        nOuterLoop = 1;
        nInnerLoop = length(demographic.ixOfPopRunSets);
        sheetNameList = {demographic.name};
        populationReportName = {demographic.reportName};
    case 'serial'
        nOuterLoop = length(demographic.ixOfPopRunSets);
        nInnerLoop = 1;
        sheetNameList = {PopRunSet(demographic.ixOfPopRunSets).name};
        populationReportName = {PopRunSet(demographic.ixOfPopRunSets).popReportName};
end

% start Loop

for iSet = 1:nOuterLoop
        
    %initilaize new sheet
    header = sprintf('Physiology of %s',populationReportName{iSet});
    FP = FP.iniCaptiontextFigtextArray(sheetNameList{iSet},header);
    
    % merge population if necessayr
    ix = demographic.ixOfPopRunSets(iSet:max(nInnerLoop,iSet));
    parValues = loadPop(Settings,{PopRunSet(ix).name},parPathSelection);

    for iY = 1:size(demographic.yList,2)

        y = parValues(:,iY);
        if ~isempty(demographic.ixPopRunSetRef)
            yRef = parValuesRef(:,iY);
        else
            yRef = [];
        end
        
        % if no property for x is give do a histogramm
        if isempty(demographic.xList)
            
            % do figure
            csv = plotLayoutHistogram(Settings,FP.figure_handle,y,yRef,demographic.yList{iY,2},demographic.yList{iY,3},...
                populationReportName{iSet},refPopulationReportName);
            
            % get name and figure description
            figureName = sprintf('h%d_%s',iY,removeForbiddenLetters(demographic.yList{iY,2}));
            figtxt = sprintf('Distribution of %s for %s.',demographic.yList{iY,2},populationReportName{iSet});
            if ~isempty(yRef)
                figtxt = sprintf('%s in comparison to %s.',figtxt(1:end-1),populationReportName{iSet});
            end
            
            % save figure
            FP = FP.printFigure(figureName,figtxt,csv,figtxt);
            
                
        else
            
            for iX = 1:size(demographic.xList,2)
                plotLayoutShadedArea(Settings,x,y,[],yRef,demographic.x_list{iY,2},demographic.x_list{iY,3},...
                    demographic.yList{iY,2},demographic.yList{iY,3});
            end
            
        end
    end
    
    
    FP.saveCaptiontextArray;
    
end
    
return
    

function parValuesFinal = loadPop(Settings,listOfname,parPathSelection)

for iName = 1:length(listOfname)

    load(fullfile('tmp',listOfname{iName},'pop.mat'));
    [jj,ix] = ismember(parPathSelection(:,1),parPaths);
    
    if sum(jj) < size(parPathSelection,2)
       error('Error path selection for population parameter, was not correct');
    end

    if iName ==1
        % initialize parValuesFinal
        parValuesFinal = parValues(:,ix(jj)); %#ok<NODEF>
        
    else
        parValuesFinal = [parValuesFinal;parValues(:,ix(jj))]; %#ok<NODEF,AGROW>
    end
end
        
% get Unit factor
unit = unit(ix(jj)); %#ok<NODEF>

for iU = 1:length(unit)
    
    switch unit{iU}
        case 'none'
            unitFactor = 1;
        otherwise
            unitFactor = getUnitFactorForUnknownDimension(Settings,unit{iU},parPathSelection{iU,2});
    end
    
    parValuesFinal(:,iU) = parValuesFinal(:,iU).*unitFactor;

    
end
 

return

