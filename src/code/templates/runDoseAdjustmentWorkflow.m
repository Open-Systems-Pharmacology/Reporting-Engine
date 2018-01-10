function runDoseAdjustmentWorkflow(WSettings,PopRunSet,Dosefit)


%% initialize workflow
[WSettings] = initializeWorkflow(WSettings);


for iSet = 1:length(PopRunSet)

    fname = fullfile('tmp',PopRunSet.name,'dosing.mat');

    
    doseFitForSet(WSettings,PopRunSet(iSet),Dosefit,fname);
    plotsForSet(WSettings,PopRunSet(iSet),Dosefit,fname);
end


function plotsForSet(WSettings,PopRunSet,Dosefit,fname)

% read population
load(fullfile('tmp',PopRunSet.name,'pop.mat'),'parPaths','parValues');
jj = strcmp('Organism|Weight',parPaths); 
weight =  parValues(:,jj); %#ok<NODEF>
jj = strcmp('Organism|BSA',parPaths);
bsa =  parValues(:,jj);

load(fname,'xMin','xMax');

lgtxt = {'all','C_{max} <= max(adult)','C_{tend} >= min(adult)','min(adult) <= AUC <= max(adult)'};


 FP = ReportFigurePrint(fullfile(WSettings.figures,'dosing'),WSettings.printFormatList);

 ax = getReportFigure(WSettings,2,1,FP.figureHandle);
 axes(ax(1))
hist(weight);
xlabel('body weight [kg]')
 axes(ax(2))
hist(bsa);
xlabel('body surface area [m^2]')
 
for iCond = 1:length(Dosefit.Condition)
    ax = getReportFigure(WSettings,2,1,FP.figureHandle);

    for iAX = 1:2
        
        switch iAX
            case 1
                y = weight;
            case 2
                y= bsa;
        end
        
        axes(ax(iAX))
        lgh = [];
        if ~isnan(Dosefit.Condition(iCond).values(1))
            jj = xMin(:,iCond)>0 & xMin(:,iCond) < Dosefit.adultDose;
            lgh(end+1) = scatter(y(jj),xMin(jj,iCond),'^g','markerfacecolor','g','displayname','minimal dose');
            jj = xMin(:,iCond) >= Dosefit.adultDose;
            scatter(y(jj),repmat(Dosefit.adultDose,1,sum(jj)),'^g','displayname','minimal dose');
            jj = xMin(:,iCond) <0;
            scatter(y(jj),repmat(-1,1,sum(jj)),'^g','displayname','minimal dose');
        end
        if ~isnan(Dosefit.Condition(iCond).values(2))
            jj = xMax(:,iCond)>0 & xMax(:,iCond) < Dosefit.adultDose;
            lgh(end+1) = scatter(y(jj),xMax(jj,iCond),'vr','markerfacecolor','r','displayname','maximal dose');
            jj = xMax(:,iCond) >= Dosefit.adultDose;
            scatter(y(jj),repmat(Dosefit.adultDose,1,sum(jj)),'vr','displayname','maximal dose');
            jj = xMin(:,iCond) <0;
            scatter(y(jj),repmat(-1,1,sum(jj)),'^g','displayname','maximal dose');
        end
        ylabel('dose [mg]')
        switch iAX
            case 1
                xlabel('Body weight [kg]');
                title(lgtxt{iCond+1})
                legend(lgh,get(lgh,'displayname'),4);
            case 2
                xlabel('Body surface area [m^2]');
        end
    end
         
end



jjInd = any(~isnan(xMin),2);
nInd = sum(jjInd);
weight = weight(jjInd);
bsa = bsa(jjInd);

jj = cellfun(@(x) ~isnan(x(1)),{Dosefit.Condition.values});
x = reshape(xMin(jjInd,jj),1,sum(jj)*nInd);
jj = cellfun(@(x) ~isnan(x(2)),{Dosefit.Condition.values});
x = [x,reshape(xMax(jjInd,jj),1,sum(jj)*nInd)];

jj = x>0 & x<= Dosefit.adultDose;
x=x(jj);
x = sort(x);


for iFlag = 2%1:2

    switch iFlag
        case 1
            % weighthBins = [0 15 inf];
            % weight_txt = {'BW<=15kg','BW>15kg'};
            weighthBins = [0 10:2:20 inf];
            for iBW = 1:length(weighthBins)-1
                weight_txt{iBW} = sprintf('BW %d-%d kg',weighthBins(iBW),weighthBins(iBW+1));
            end
            weight_txt{1} = sprintf('BW <= %d kg',weighthBins(2));
            weight_txt{end} = sprintf('BW > %d kg',weighthBins(end-1));
            
            prop = weight;
        case 2
            weighthBins = [0 0.5:0.1:0.8 inf];
            weight_txt = {};
            for iBW = 1:length(weighthBins)-1
                weight_txt{iBW} = sprintf('BSA %d-%d m^2',weighthBins(iBW),weighthBins(iBW+1));
            end
            weight_txt{1} = sprintf('BSA <= %d m^2',weighthBins(2));
            weight_txt{end} = sprintf('BSA > %d m^2',weighthBins(end-1));

            prop = bsa;
    end
    y=[];
    for iBW = 1:length(weight_txt)
        
        axBW = getReportFigure(WSettings,1,1,FP.figureHandle,'figureformat','landscape');
        lgh = [];
        
        for iCond = 0:length(Dosefit.Condition)
            
            for iX = 1:length(x)
                jjX = prop>weighthBins(iBW) & prop<= weighthBins(iBW+1);
                nBin = sum(jjX);
                
                iCList  = find(cellfun(@(x) ~isnan(x(1)),{Dosefit.Condition.values}));
                if iCond >0
                    iCList = intersect(iCList,iCond);
                end
                for iC = iCList
                    jjX = jjX & xMin(jjInd,iC) <= x(iX);
                end
                
                iCList = find(cellfun(@(x) ~isnan(x(2)),{Dosefit.Condition.values}));
                if iCond >0
                    iCList = intersect(iCList,iCond);
                end
                for iC = iCList
                    jjX = jjX & xMax(jjInd,iC) >= x(iX);
                end
                
                y(iX,iCond+1,iBW) = sum(jjX)./nBin*100;
            end
        end
    end
    
    colfield = getcolmarkForMap(jet,length(lgtxt));
    colBW = getcolmarkForMap(jet,length(weight_txt));
    
    for iBW = 100:length(weight_txt)
        
        ax = getReportFigure(WSettings,1,1,FP.figureHandle,'figureformat','landscape');
        
        for iCond = 0:length(Dosefit.Condition)
            plot(ax,x,y(:,1+iCond,iBW),'-','linewidth',2,'color',colfield(iCond+1,:),...
                'displayname',sprintf('%s',lgtxt{iCond+1}));
        end
        
        title(sprintf('schoolchildren %s',weight_txt{iBW}));
        xlabel(ax,'Dose [mg]')
        ylabel(ax,'Percent of children which fullfiles condition')
        set(ax,'xlim',[0 10])
        l = legend(ax,'show');
        
    end
    
    
    cut = 70;
    ax = getReportFigure(WSettings,1,1,FP.figureHandle,'figureformat','landscape');
    for iBW = 1:length(weight_txt)
        lgh(iBW) = plot(ax,x,y(:,1,iBW),'-','linewidth',2,'color',colBW(iBW,:),...
            'displayname',sprintf('%s',weight_txt{iBW}));
        ijcut = find(y(:,1,iBW) > cut);
        plot(ax,x(ijcut([1 1])),[0 cut],'--','linewidth',1,'color',colBW(iBW,:));
        plot(ax,x(ijcut([end end])),[0 cut],'--','linewidth',1,'color',colBW(iBW,:));
    end
    
    
    xlabel(ax,'Dose [mg]')
    ylabel(ax,'Percent of children which fullfiles condition')
    set(ax,'xlim',[0 10])
    legend(lgh,get(lgh,'displayname'));
end


return
   
function  doseFitForSet(WSettings,PopRunSet,Dosefit,fname)
    

writeToReportLog('INFO',sprintf('Start evaluating of %s',PopRunSet.name),false);

% read population
load(fullfile('tmp',PopRunSet.name,'pop.mat'),'parPaths','parValues');

% load OutputList
load(fullfile('tmp',PopRunSet.name,'outputList.mat'),'OutputList');

% load OutpulList
load(fullfile('tmp',PopRunSet.name,'applicationProtocol.mat'),'ApplicationProtocol','isValid');



jj = strcmp('Organism|Weight',parPaths);  %#ok<USENS>
weight =  parValues(:,jj); %#ok<NODEF>
jj = strcmp('Organism|Height',parPaths);
height =  parValues(:,jj);
jj = strcmp('Organism|BSA',parPaths);
bsa =  parValues(:,jj);


% initialize model
initStruct=[];
for iPar=1:length(parPaths) 
    initStruct=initParameter(initStruct,[ '*|' parPaths{iPar}],'always','throwWarningIfNotExisting',false);
end
initStruct=initParameter(initStruct,[ '*|' Dosefit.Parameter],'always','throwWarningIfNotExisting',false);
% initSimulation(PopRunSet.xml,initStruct);

% Only one simulation is intialized so the simulation index is always 1
simulationIndex = 1;

% set simulaion time if given
if isfield(PopRunSet,'simulationTime') && ~isempty(PopRunSet.simulationTime)
    setSimulationTime(PopRunSet.simulationTime,simulationIndex)
end

% getRowindex
rowIndex=nan(size(parPaths));
for iPar=1:length(parPaths)
    [ise,~,tmp] = existsParameter(['*|' parPaths{iPar}],simulationIndex);
    if ise
        rowIndex(iPar) = tmp;
    end
end
[ise,~,rowIndexFit] = existsParameter(['*|' Dosefit.Parameter],simulationIndex);
if ~ise
    error('Fit parameter %s does not exist');
end
rowIndexS = nan(size(OutputList));

% get ReferenceValue
refValue = getParameter(['*|' Dosefit.Parameter],simulationIndex,'rowIndex',rowIndexFit,'parametertype','reference');
refRange = refValue.*Dosefit.refRange;

% get individual ids
nInd = size(parValues,1);

if exist(fname,'file')
    load(fname,'xMin','xMax');
    iIndStart = find(all(isnan(xMin),2),1);
else
    xMin = nan(nInd,length(Dosefit.Condition));
    xMax = nan(nInd,length(Dosefit.Condition));
    iIndStart = 1;
end

for iInd = iIndStart:nInd
    % with row index use fast methode
    for iPar=find(~isnan(rowIndex))
        setParameter(parValues(iInd,iPar),'',simulationIndex,'speedy','variable',rowIndex(iPar));
    end

    fitValue = refValue./weight(iInd);
    fitValueToEvaluate = max(fitValue/2,Dosefit.refRange(1).*fitValue);
    goOn = true;
    fPar = nan(length(Dosefit.Condition),5);
    
    while goOn
        
        iRun = length(fitValue);
        
        setParameter(fitValue(iRun),['*|' Dosefit.Parameter],simulationIndex,'rowIndex',rowIndexFit);
        
        success=processSimulation;
        
        if success
            % get result
            for iO=1:length(rowIndexS)
                [simTime,simValues(:,1),rowIndexS(iO)] = getSimulationResult(['*|' OutputList(iO).pathID],simulationIndex,'rowindex',rowIndexS(iO));
                
                tmp = calculatePKParameterForApplicationProtocol(WSettings,ApplicationProtocol,simTime,simValues,weight(iInd),height(iInd));
                
                jj = [Dosefit.Condition.OutputID] ==iO;
                [jj,ix] = ismember({Dosefit.Condition(jj).field},{tmp.name});
                PKPList(jj) = tmp(ix(jj));
            end
            
            for iCond = 1:length(Dosefit.Condition)
                pK(iCond,iRun) = PKPList(iCond).value;
                
            end
        end
        
        % sort fitValues
        [fitValue,ix] = sort(fitValue);
        pK = pK(:,ix);
        
        
        % set predifined new value
        if ~isempty(fitValueToEvaluate)
            fitValue(end+1) = fitValueToEvaluate(1); %#ok<AGROW>
            fitValueToEvaluate = fitValueToEvaluate(2:end);
        else
            % evaluate conditions check for lower boundary
            jjCond = isnan(xMin(iInd,:)) & cellfun(@(x) ~isnan(x(1)),{Dosefit.Condition.values});
            if any(jjCond)
                for iCond = find(jjCond)
                    [newValuesToEvaluate,xMinTmp,fPar(iCond,:)] = ...
                        checkAcceptance(fitValue,pK(iCond,:),Dosefit.Condition(iCond).values(1),refRange,fPar(iCond,:));
                    if ~isempty(xMinTmp)
                        xMin(iInd,iCond) = xMinTmp; %#ok<AGROW>
                    end
                    fitValueToEvaluate = [fitValueToEvaluate, newValuesToEvaluate]; %#ok<AGROW>
                    
                end
            end
            % evaluate conditions check for upper boundary
            jjCond = isnan(xMax(iInd,:)) & cellfun(@(x) ~isnan(x(2)),{Dosefit.Condition.values});
            if any(jjCond)
                for iCond = find(jjCond)
                    % multiply with -1 to make upper boreder to lower
                    % border
                    [newValuesToEvaluate,xMaxTmp,fPar(iCond,:)] = ...
                        checkAcceptance(fitValue,pK(iCond,:),Dosefit.Condition(iCond).values(2),refRange,fPar(iCond,:));
                    if ~isempty(xMaxTmp)
                        xMax(iInd,iCond) = xMaxTmp;
                    end
                    fitValueToEvaluate = [fitValueToEvaluate, newValuesToEvaluate]; %#ok<AGROW>
                end
            end
            if isempty(fitValueToEvaluate)
                goOn = false;
            else
                fitValueToEvaluate = unique(fitValueToEvaluate);
            end
        end
        
        if iRun ==100
            goOn = false;
            save(fullfile('tmp',PopRunSet.name,sprintf('dosing_Ind%d.mat',iInd)),'fitValue','pK');
        end
        
    end
    
    ResInd(iInd).fPar = fPar;
    ResInd(iInd).fitValue = fitValue;
    ResInd(iInd).pK = pK;
    
    disp(sprintf('finished Ind %d with %d runs',iInd,iRun));
    save(fname,'xMin','xMax','ResInd');
    
end


return

function [newValuesToEvaluate,xMin,p] = checkAcceptance(fitValue,y,cutValue,refRange,p)
                 
newValuesToEvaluate = [];
xMin = [];



% CutValue part of y already ?
[minDist,iMin] = min(abs(cutValue-y(:)));
if abs(minDist./min(cutValue,y(iMin))) <= 0.01
    xMin =  fitValue(iMin);
    return
end

% estimate for x for cutValue
switch length(y)
    case 1
        error('Start with at least two values')
    case 2   
        p(1:2) = polyfit(fitValue,y,1);
        newValuesToEvaluate = (cutValue-p(2))/p(1);
    otherwise

        % get the slope direction of pk vs dose
        p(5) = sign(y(end) - y(1));
        y = y.^p(5);


        % dimension of p
        nP = min(length(y),4);
        % start value
        if any(isnan(p))
            p(1) = 0;
            p(2) = max(y);
            p(3) = mean(fitValue);
        end
        if isnan(p(4))
            p(4) = 1;
        end
        

        % get new fitparameter
        p(1:nP) = fminsearch(@(p) emaxHill(p,fitValue,y),p(1:nP));
        
        % check if cutValue can be reached
        if cutValue >= p(2)
            newValuesToEvaluate = refRange(2);
        elseif cutValue <= p(1)
            newValuesToEvaluate = refRange(1);
        else
            %get new Value from fit
            newValuesToEvaluate = inverseEmaxHill(p,cutValue);
        end
end
% limit Value to ranges
newValuesToEvaluate = min(max(newValuesToEvaluate,refRange(1)),refRange(2));

% already tested ?
[minDist,iMin] = min(abs(newValuesToEvaluate-fitValue(:)));
if abs(minDist./newValuesToEvaluate) <= 0.01
    newValuesToEvaluate = [];
    xMin =  -fitValue(iMin);
    return
end

return


function [newValuesToEvaluate,xMin] = checkAcceptance2(fitValue,y,cutValue,refRange)
                 
newValuesToEvaluate = [];
xMin = [];

% get the slope direction of pk vs dose
isIncreasing = y(1) < y(end);

% check for boundary
accept =  y <= cutValue;

% check if right side has to be increased
if (all( accept) && isIncreasing)  || (~any( accept) && ~isIncreasing)
    sens = range(y(end-1:end))./y(end-1) ./ (range(fitValue(end-1:end))./fitValue(end-1));
    if abs(fitValue(end) - refRange(end))/refRange(2)<0.01;
        xMin = refRange(2);
    elseif abs(sens)>0.05
        newValuesToEvaluate(end+1) = min(2*fitValue(end),refRange(2));
    else
        xMin = -1*fitValue(end);
    end
    % or if left side has to be increased
elseif (all( accept) && ~isIncreasing)  || (~any( accept) && isIncreasing)
    sens = range(y(1:2))./y(2) ./ (range(fitValue(1:2))./fitValue(2));
    if abs(fitValue(end) - refRange(1))/refRange(1)<0.01;
        xMin = refRange(1);
    elseif abs(sens)>0.05
        newValuesToEvaluate(end+1) = max(fitValue(1)/2,refRange(1));
    else
        xMin = -1*fitValue(1);
    end
    % or if if the value is within
else
    % already there?
    [minDist,iMin] = min(abs(cutValue-y(:)));
    if abs(minDist./cutValue) > 0.01
        newValuesToEvaluate(end+1) = interp1(y,fitValue,cutValue);
    else
        xMin =  fitValue(iMin); 
    end
end

return

function err = emaxHill(p,x,y)

yMin = p(1);
yMax = p(2);
ec50 = p(3);
if length(p)==4
    hill = p(4);
else
    hill = 1;
end

yfit = yMin + (yMax-yMin)./(1+(x./ec50).^hill);

err = norm(yfit-y);
return

function x = inverseEmaxHill(p,y)

yMin = p(1);
yMax = p(2);
ec50 = p(3);
if length(p)==4
    hill = p(4);
else
    hill = 1;
end

if y<yMin || y>yMax
    x=nan;
else
    x =  (((yMax-yMin)/(y-yMin)-1)^(1/hill)).*ec50;
end
