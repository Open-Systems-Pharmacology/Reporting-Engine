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
    
    FP = doShadedArea(WSettings,Def,textFunctionHandle,PopRunSet,FP,'absolute');
    
    if ~isempty(Def.ixPopRunSetRef)
        FP = doShadedArea(WSettings,Def,textFunctionHandle,PopRunSet,FP,'ratio');
    end

end

return

function FP = doBWhiskers(WSettings,Def,textFunctionHandle,PopRunSet,FP,flag)

% collect all simulated PK Parameter of interst
o = collectPK(PopRunSet,'matrix');

% get Reference simulation
if ~isempty(Def.ixPopRunSetRef)
    oRef = collectPK(PopRunSet(Def.ixPopRunSetRef),'matrix');
    reportNameRef = PopRunSet(Def.ixPopRunSetRef).reportName;
else
    oRef = [];
    reportNameRef = '';
end

switch flag
        case 'absolute'
            reportNames = {PopRunSet.reportName};
            popReportNames = {PopRunSet.popReportName};
            xLabels = {PopRunSet.boxwhiskerLabel};
    case 'ratio'
            reportNames = {PopRunSet(Def.ixOfPopRunSets).reportName};
            popReportNames = {PopRunSet(Def.ixOfPopRunSets).popReportName};
            xLabels = {PopRunSet(Def.ixOfPopRunSets).boxwhiskerLabel};        
    otherwise
        error('unknown flag');
end

yscale = {'lin','log'};

for iO = 1:length(o)
    
    %initialize new sheet
    switch flag
        case 'absolute'
            header = sprintf('PK parameter of %s',o(iO).reportName);
            sheet = removeForbiddenLetters(header);
        case 'ratio'
            header = sprintf('PK parameter of %s as fraction of %s',o(iO).reportName,lower(reportNameRef));
            sheet = removeForbiddenLetters(sprintf('PK_%s_as_fraction',o(iO).reportName));
        otherwise
            error('unknown flag');
    end
    FP = FP.iniCaptiontextFigtextArray(header,sheet);
    
    for iPK = 1:size(o(iO).pKParameterList,2)
         
         FP = FP.addSubSection(o(iO).pKParameterList{4,iPK},2);
        
        switch flag
            case 'absolute'
                % get Vector of PK Parameter
                y = squeeze(o(iO).values(:,:,iPK));
                yLabel = o(iO).pKParameterList{4,iPK};
                yUnit = o(iO).pKParameterList{2,iPK};
                
                % get popPKValues
                popPKValues = nan;
                popPKReference = '';
                if isfield(Def,'popPK')
                    jjO = [Def.popPK.values(:).iO] == iO;
                    jj = strcmp(o(iO).pKParameterList{1,iPK},{Def.popPK.values(jjO).par});
                    if any(jj)
                        popPKValues = Def.popPK.values(jj).v;
                        popPKReference = Def.popPK.reference;
                    end
                end
                
                
                
            case 'ratio'
                
                % get popPKValues
                popPKValues = nan;
                popPKReference = '';

                
                jPKref = strcmp(o(iO).pKParameterList{1,iPK},oRef(iO).pKParameterList(1,:));
                
                if any(jPKref)
                    y = squeeze(o(iO).values(:,Def.ixOfPopRunSets,iPK))./repmat(squeeze(oRef(iO).values(:,:,iPK)),1,length(Def.ixOfPopRunSets));
                    yLabel = o(iO).pKParameterList{4,iPK};
                    yUnit = sprintf('fraction of %s',lower(reportNameRef));
                else
                    y=[];
                end
            otherwise
                error('unknown flag');

        end
        
        if ~isempty(y)
            % loop on scale
            for iScale = 1:length(yscale)
                
                switch flag
                    case 'absolute'
                        % get name and figure description
                        figureName = sprintf('bw%d_%s_%s',iPK,removeForbiddenLetters(o(iO).pKParameterList{1,iPK}),yscale{iScale});
                        [figtxt,figtxtTable] = feval(textFunctionHandle,WSettings,'pkBW',...
                            {yLabel,o(iO).reportName,reportNames,popReportNames,yscale{iScale},popPKReference});
                        
                    case 'ratio'
                        % get name and figure description
                        figureName = sprintf('bwRatio%d_%s_%s',iPK,removeForbiddenLetters(o(iO).pKParameterList{1,iPK}),yscale{iScale});
                        [figtxt,figtxtTable] = feval(textFunctionHandle,WSettings,'pkBWRatio',...
                            {yLabel,o(iO).reportName,reportNames,popReportNames,yscale{iScale},reportNameRef});
                        
                    otherwise
                        error('unknown flag');


                end
                
                % do figure
                [csv] = plotReportBoxwhisker(WSettings,FP.figureHandle,y,yLabel,yUnit,yscale{iScale},xLabels,popPKValues,popPKReference);
                
                % save figure
                if iScale == length(yscale)
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

    
function FP = doShadedArea(WSettings,Def,textFunctionHandle,PopRunSet,FP,flag)

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
        switch flag
            case 'absolute'
                header = sprintf('%s%s-dependence of %s',upper(Def.xList{iX,2}(1)),Def.xList{iX,2}(2:end),o(iO).reportName);
                sheet = sprintf('%s_%s',Def.xList{iX,2},o(iO).reportName);
            case 'ratio'
                header = sprintf('%s%s-dependence of %s as fraction of %s',...
                    upper(Def.xList{iX,2}(1)),Def.xList{iX,2}(2:end),o(iO).reportName,reportNameRef);                
                sheet = removeForbiddenLetters(sprintf('%s-dependence of %s_fraction',Def.xList{iX,2},o(iO).reportName));
            otherwise
                error('unknown flag')
        end
        FP = FP.iniCaptiontextFigtextArray(header,sheet);

        
        jX = strcmp(Def.xList{iX,1},parPathSelection);
        x = parValues(:,jX);
        
        for iPK = 1:size(o(iO).pKParameterList,2)
            
            % add new subheader
             FP = FP.addSubSection(o(iO).pKParameterList{4,iPK},2);
            
            if isempty(oRef)
                yRef = [];
            else
                jPKref = strcmp(o(iO).pKParameterList{1,iPK},oRef(iO).pKParameterList(1,:));
                if any(jPKref)
                    yRef = oRef(iO).values(:,jPKref);
                else
                    yRef = [];
                end
            end
       
            
            y = o(iO).values(:,iPK);
            yLabel = o(iO).pKParameterList{4,iPK};
            switch flag
                case 'absolute'
                    yUnit = o(iO).pKParameterList{2,iPK};
                case 'ratio'
                    
                    if ~isempty(yRef)
                        [~,yMeanRef,~,legendTextMean] = getRangePlotPercentiles(WSettings,yRef);
                        y = y./yMeanRef*100;
                        yUnit = '%';                        
                    else
                        y=[];
                    end
                otherwise
                    error('unknown flag');
                    
            end
            
            
            
            
            if warningflagRangePlots
                writeToReportLog('WARNING',sprintf(['you are creating Rangeplots with %d individuals, recommended are at least %d individuals. ',...
                    'Statistic might be not sufficient'],length(y),WSettings.rangePlotsMin),false);
                warningflagRangePlots = false;
            end
            
            % get popPKValues
            popPKValues = nan;
            popPKReference = '';
            if isfield(Def,'popPK')
                popPKReference = Def.popPK.reference;
                jjO = [Def.popPK.values(:).iO] == iO;
                jj = strcmp(o(iO).pKParameterList{1,iPK},{Def.popPK.values(jjO).par});
                if any(jj)
                    popPKValues = Def.popPK.values(jj).v;
                    popPKReference = Def.popPK.reference;
                end
            end

            % loop on scale
            for iScale = 1:length(yscale)

                
                switch flag
                    case 'absolute'
                        % get name and figure description
                        figureName = sprintf('shA%d_%d_%s_%s_%s',iX,iPK,removeForbiddenLetters(Def.xList{iX,2}),...
                            removeForbiddenLetters(o(iO).pKParameterList{1,iPK}),yscale{iScale});
                        
                        [figtxt,figtxtTable,legendEntries] = feval(textFunctionHandle,WSettings,'pkShA',...
                            {Def.xList{iX,2},yLabel,o(iO).reportName,Def.reportName,reportNameRef,...
                            [length(y),length(yRef)],yscale{iScale},popPKReference});
                        
                        csv = plotReportShadedArea(WSettings,FP.figureHandle,x,y,yRef,Def.xList{iX,2},Def.xList{iX,3},...
                            yLabel,yUnit,yscale{iScale},legendEntries,popPKValues);

                        
                    case 'ratio'
                        % get name and figure description
                        figureName = sprintf('shARatio%d_%d_%s_%s_%s',iX,iPK,removeForbiddenLetters(Def.xList{iX,2}),...
                            removeForbiddenLetters(o(iO).pKParameterList{1,iPK}),yscale{iScale});

                        [figtxt,figtxtTable,legendEntries] = feval(textFunctionHandle,WSettings,'pkShARatio',...
                            {Def.xList{iX,2},yLabel,o(iO).reportName,Def.reportName,reportNameRef, ...
                            [length(y),length(yRef)],yscale{iScale},'',legendTextMean});
                                                
                        csv = plotReportShadedArea(WSettings,FP.figureHandle,x,y,[],Def.xList{iX,2},Def.xList{iX,3},...
                            yLabel,yUnit,yscale{iScale},legendEntries,nan);

                    otherwise
                        error('unknown flag');
                end

                

                
                % save figure
                if iScale == length(yscale)                    
                    FP = FP.printFigure(figureName,figtxt,csv,figtxtTable);
                else
                    FP = FP.printFigure(figureName,figtxt);
                end
            end

            
        end
        
        % save caption array for this sheet
        FP.saveCaptiontextArray;

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
        
        for iO = 1:length(PKPList)
            
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
                            o(iO).values(:,iSet,iPK) = [PKPList{iO}(:,iPK).value].*OutputList(iO).pKParameterList{3,iPK};
                        case 'vector'
                            o(iO).values(1:length([PKPList{iO}(:,iPK).value]),iPK) = ...
                                [PKPList{iO}(:,iPK).value].*OutputList(iO).pKParameterList{3,iPK}; 
                        otherwise
                            error('unknown flag');
                    end
                end
            end
                
                
        end
    else
        for iO = 1:length(PKPList)
            
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
                nInd = length([PKPList{iO}(:,1).value]);
                nIndOld = size(o(iO).values,1);
                if strcmp(flag,'matrix')
                        if nInd > nIndOld
                            o(iO).values(nIndOld+1:nInd,:,:) = nan; %#ok<AGROW>
                        end
                        
                        tmp = nan(max(nInd,nIndOld),1);
                end

                for k = 1:length(ixTarget)
                    tmp(1:nInd) = [PKPList{iO}(:,ixSource(k)).value]...
                        .*OutputList(iO).pKParameterList{3,ixTarget(k)};
                    switch flag
                        case 'matrix'                            
                            o(iO).values(:,iSet,ixTarget(k)) = tmp; %#ok<AGROW>
                        case 'vector'
                            o(iO).values(nIndOld + [1:nInd],ixTarget(k)) = tmp; %#ok<NBRAK,AGROW>
                            
                        otherwise
                            error('unknown flag');

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
    load(fname); %#ok<LOAD>
else
    PKPList = readPKParameter(simulationName);
end

return
