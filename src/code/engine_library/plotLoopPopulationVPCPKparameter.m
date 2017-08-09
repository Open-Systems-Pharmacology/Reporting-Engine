function FP = plotLoopPopulationVPCPKparameter(WSettings,textFunctionHandle,Def,PopRunSet,FP) 
% PLOTLOOPPOPULATIONVPCPKPARAMETER generates pK Parameter plots for VPC
%
% FP = plotLoopPopulationVPCPKparameter(WSettings,Def,PopRunSet,FP) 
%
% Inputs: 
%       Settings (structure)    definition of properties used in all
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

if isempty(Def.xList)
    
    FP = doBWhiskers(WSettings,Def,textFunctionHandle,PopRunSet,FP,'absolute');

    if ~isempty(Def.ixPopRunSetRef)
        FP = doBWhiskers(WSettings,Def,textFunctionHandle,PopRunSet,FP,'ratio');
    end
else
    
    FP = doShadedArea(WSettings,Def,textFunctionHandle,PopRunSet,FP);
end

return

function FP = doBWhiskers(WSettings,Def,textFunctionHandle,PopRunSet,FP,flag)

% collect all simulated PK Parameter of interst
o = collectPK(PopRunSet(Def.ixOfPopRunSets),'matrix');

% get Reference simulation
if ~isempty(Def.ixPopRunSetRef)
    oRef = collectPK(PopRunSet(Def.ixPopRunSetRef),'matrix');
    reportNameRef = PopRunSet(Def.ixPopRunSetRef).reportName;
else
    oRef = [];
    reportNameRef = '';
end

reportNames = {PopRunSet(Def.ixOfPopRunSets).reportName};
popReportNames = {PopRunSet(Def.ixOfPopRunSets).popReportName};
xLabels = {PopRunSet(Def.ixOfPopRunSets).boxwhiskerLabel};

yscale = {'lin','log'};

for iO = 1:length(o)
    
    %initialize new sheet
    header = sprintf('PK parameter of %s',o(iO).reportName);
    sheet = removeForbiddenLetters(header);
    FP = FP.iniCaptiontextFigtextArray(header,sheet);
    
    for iPK = 1:size(o(iO).pKParameterList,2)
         
        switch flag
            case 'absolute'
                y = squeeze(o(iO).values(:,:,iPK));
                yLabel = o(iO).pKParameterList{4,iPK};
                yUnit = o(iO).pKParameterList{2,iPK};
            case 'ratio'
                
                jPKref = strcmp(o(iO).pKParameterList{1,iPK},oRef(iO).pKParameterList{1,iPK});
                
                if any(jPKref)
                    y = squeeze(o(iO).values(:,:,iPK))./squeeze(oRef(iO).values(:,:,jPKref));
                    yLabel = o(iO).pKParameterList{4,iPK};
                    yUnit = o(iO).pKParameterList{2,iPK};
                else
                    y=[];
                end
        end
        
        if ~isempty(y)
            % loop on scale
            for iScale = 1:length(yscale)
                
                switch flag
                    case 'absolute'
                        % get name and figure description
                        figureName = sprintf('bw%d_%s_%s',iPK,removeForbiddenLetters(o(iO).pKParameterList{1,iPK}),yscale{iScale});
                        [figtxt,figtxtTable] = feval(textFunctionHandle,WSettings,'pkBW',...
                            {yLabel,o(iO).reportName,reportNames,popReportNames,yscale{iScale}});
                        
                    case 'ratio'
                        % get name and figure description
                        figureName = sprintf('bwRatio%d_%s_%s',iPK,removeForbiddenLetters(o(iO).pKParameterList{1,iPK}),yscale{iScale});
                        [figtxt,figtxtTable] = feval(textFunctionHandle,WSettings,'pkBWRatio',...
                            {yLabel,o(iO).reportName,reportNames,popReportNames,yscale{iScale}});
                end
                
                % do figure
                [csv] = plotReportBoxwhisker(WSettings,FP.figureHandle,y,yLabel,yUnit,yscale{iScale},xLabels);
                
                % save figure
                if iScale ==1
                    FP = FP.printFigure(figureName,figtxt,csv,figtxtTable);
                else
                    FP = FP.printFigure(figureName,figtxt);
                end
            end
        end
    end

    % save caption array for this sheet
    FP.saveCaptiontextArray;

end

return

    
function FP = doShadedArea(WSettings,Def,textFunctionHandle,PopRunSet,FP)

o = collectPK(PopRunSet(Def.ixOfPopRunSets),'vector');

    
% get Reference simulation
if ~isempty(Def.ixPopRunSetRef)
    oRef = collectPK(PopRunSet(Def.ixPopRunSetRef),'vector');
    reportNameRef = PopRunSet(Def.ixPopRunSetRef).reportName;
else
    oRef = [];
    reportNameRef = '';
end

% get selectionList of population parameter
parPathSelection = Def.xList(:,[1 3]);

% merge population
[parValues] = loadSelectedPathsOfPopulation(WSettings,{PopRunSet(Def.ixOfPopRunSets).name},parPathSelection);
if size(parValues,1) < WSettings.rangePlotsMin
    warningflagRangePlots = true;
else
    warningflagRangePlots = false;
end


yscale = {'lin','log'};

for iO = 1:length(o)
    
    for iX = 1:size(Def.xList,1)

        %initialize new sheet
        header = sprintf('%s-dependence of %s',Def.xList{iX,2},o(iO).reportName);
        sheet = sprintf('%s_%s',Def.xList{iX,2},o(iO).reportName);
        FP = FP.iniCaptiontextFigtextArray(header,sheet);
        
        jX = strcmp(Def.xList{iX,1},parPathSelection);
        x = parValues(:,jX);
        
        for iPK = 1:size(o(iO).pKParameterList,2)
            
            y = o(iO).values(:,iPK);                
            yLabel = o(iO).pKParameterList{4,iPK};
            yUnit = o(iO).pKParameterList{2,iPK};
            
            if isempty(oRef)
                yRef = [];
            else
                yRef = oRef(iO).values(:,iPK);         
            end
       
            if warningflagRangePlots
                writeToLog(sprintf(['WARNING: you are creating Rangeplots with less than %d individual. ',...
                    'Statistic might be not sufficient'],WSettings.rangePlotsMin),WSettings.logfile,true,false);
                warningflagRangePlots = false;
            end

            % loop on scale
            for iScale = 1:length(yscale)

                % get name and figure description
                figureName = sprintf('shA%d_%d_%s_%s_%s',iX,iPK,removeForbiddenLetters(Def.xList{iX,2}),...
                    removeForbiddenLetters(o(iO).pKParameterList{1,iPK}),yscale{iScale});
                [figtxt,figtxtTable,legendEntries] = feval(textFunctionHandle,WSettings,'pkShA',...
                    {Def.xList{iX,2},yLabel,o(iO).reportName,Def.reportName,reportNameRef,yscale{iScale}});

                csv = plotReportShadedArea(WSettings,FP.figureHandle,x,y,yRef,Def.xList{iX,2},Def.xList{iX,3},...
                    yLabel,yUnit,yscale{iScale},legendEntries);
                
                % save figure
                FP = FP.printFigure(figureName,figtxt,csv,figtxtTable);
            end

            
        end
        
        
    end
end


return


function o = collectPK(PopRunSet,flag)

o = [];
for iSet = 1:length(PopRunSet)
    
    load(fullfile('tmp',PopRunSet(iSet).name,'outputList.mat'),'OutputList');
    PKPList = loadPKPList(PopRunSet(iSet).name);
    
    load(fullfile('tmp',PopRunSet(iSet).name,'applicationProtocol.mat'),'PKParameterTemplate')

    if isempty(o)
        
        for iO = 1:length(PKPList); 
            
            o(iO).reportName = OutputList(iO).reportName; %#ok<AGROW>
            
            o(iO).pKParameterList = OutputList(iO).pKParameterList; %#ok<AGROW>
            
            o(iO).values = []; %#ok<AGROW>
            
            % add displayname
            if ~isempty( OutputList(iO).pKParameterList)
                [~,ix] = ismember(o(iO).pKParameterList(1,:),{PKParameterTemplate.name});
                o(iO).pKParameterList(4,:) = {PKParameterTemplate(ix).reportName}; %#ok<AGROW>

                for iPK = 1:length(PKPList{iO})
                    switch flag
                        case 'matrix'
                            o(iO).values(:,iSet,iPK) = PKPList{iO}(iPK).value.*OutputList(iO).pKParameterList{3,iPK};%#ok<AGROW>
                        case 'vector'
                            o(iO).values(1:length(PKPList{iO}(iPK).value),iPK) = ...
                                PKPList{iO}(iPK).value.*OutputList(iO).pKParameterList{3,iPK}; %#ok<AGROW>
                    end
                end
            end
                
                
        end
    else
        for iO = 1:length(PKPList); 
            
            if ~isempty(PKPList{iO})
                % get indices
                [jj,ixTarget] = ismember(OutputList(iO).pKParameterList(1,:),o(iO).pKParameterList(1,:));
                ixSource = [find(jj) find(~jj)];
                ixTarget = ixTarget(jj);
                
                % add new values to structure
                for iNew = find(~jj)
                    ixTarget(end+1) =  size(o(iO).pKParameterList,2)+1;%#ok<AGROW>
                    o(iO).pKParameterList(1:3,ixTarget(end)) = OutputList(iO).pKParameterList(:,iNew);%#ok<AGROW>
                    
                    % add displayname
                    [~,ix] = ismember(o(iO).pKParameterList(1,iNew),{PKParameterTemplate.name});
                    o(iO).pKParameterList{4,ixTarget(end)} = {PKParameterTemplate(ix).reportName}; %#ok<AGROW>
                    
                end
                
                % add values
                nInd = length(PKPList{iO}(1).value);
                nIndOld = size(o(iO).values,1);
                switch flag
                    case 'matrix'
                        if nInd > nIndOld
                            o(iO).values(nIndOld+1:nInd,:,:) = nan; %#ok<AGROW>
                        end
                        
                        tmp = nan(max(nInd,nIndOld),1);
                end

                for k = 1:length(ixTarget)
                    tmp(1:nInd) = PKPList{iO}(ixSource(k)).value...
                        .*OutputList(iO).pKParameterList{3,ixTarget(k)};
                    switch flag
                        case 'matrix'                            
                            o(iO).values(:,iSet,ixTarget(k)) = tmp; %#ok<AGROW>
                        case 'vector'
                            o(iO).values(nIndOld+[1:nInd],ixTarget(k)) = tmp; %#ok<AGROW>
                    end
                end
            end
        end
    end
        
end

return

function PKPList = loadPKPList(simulationName)

fname = fullfile('tmp',simulationName,'pKPList.mat');

if exist(fname,'file')
    load(fname);
else
    PKPList = readPKParameter(simulationName);
end

return
