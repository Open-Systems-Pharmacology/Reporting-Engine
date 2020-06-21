function runPopulationSensitivity(WSettings,PopRunSet)
% RUNPOPULATIONSENSITIVITY calculates and plots senstivity for all poulations
%
% runPopulationSensitivity(WSettings,PopRunSet)
%
% Inputs
%       WSettings (structure)    definition of properties used in all
%                   workflow functions see GETDEFAULTWORKFLOWSETTINGS
%       PopRunSet (structure)   list of population simulations see GENERATEWORKFLOWINPUTFORPOPULATIONSIMULATION

% Open Systems Pharmacology Suite;  http://open-systems-pharmacology.org

try
    writeToReportLog('INFO','start sensitivity analysis',false);
    
    % percentiles for individual selection
    prctileSelection = WSettings.sensitivitPrctileSelection;
    
    %if strcmp(WSettings.workflowMode,'pediatric')
    if strcmp(WSettings.workflowMode,'ratioComparison')
        jjRef = [PopRunSet.isReference];
        setList = find(~jjRef);
    else
        setList = 1:length(PopRunSet);
    end
    
    % create resultdirector
    % Initialize figureDir
    FP = ReportFigurePrint(fullfile(WSettings.figures,'sensitivity'),WSettings.printFormatList);
    
    % get outputs and PKPList
    load(fullfile('tmp',PopRunSet(setList(1)).name,'outputList.mat'),'OutputList');
    % get corresponding application protocoll
    load(fullfile('tmp',PopRunSet(setList(1)).name,'applicationProtocol.mat'),'ApplicationProtocol','isValid','PKParameterTemplate');
    
    
    for iO = 1:length(OutputList)
        
        for iPK =  1:size(OutputList(iO).pKParameterList,2)
            
            % in restart mode use previously produced arrays
            sensFilename = fullfile('tmp',sprintf('sens_O%d_PK%d.mat',iO,iPK));
            sens = cell(length(prctileSelection),length(PopRunSet)); 
            if WSettings.restart
                if exist(sensFilename,'file')
                    load(sensFilename,'sens');
                    writeToReportLog('INFO',sprintf('restart: reload %s',sensFilename),false);
                end
            end

            
            
            for iSet = setList
                
                tmpDir = fullfile('tmp',PopRunSet(iSet).name);
                PKPList = readPKParameter(PopRunSet(iSet).name);
                PKParameterList = PKPList{iO};
                
                % get corresponding population
                load(fullfile(tmpDir,'pop.mat'),'parPaths','parValues');
                jj = strcmp(parPaths,'IndividualId');
                individualVector = parValues(:,jj); %#ok<IDISVAR,NODEF>
                jj = strcmp('Organism|Weight',parPaths);
                weight =  parValues(:,jj);
                jj = strcmp('Organism|Height',parPaths);
                height =  parValues(:,jj);
                
                
                % load sensitivity parameter
                load(fullfile(tmpDir,'sensitivity.mat'),'sensParameterList');
                
                % get percentile individual
                PKrefValues = [PKParameterList(:,iPK).value];
                tmp = prctile(PKrefValues,prctileSelection);
                ixInd = nan(size(prctileSelection));
                for iPct = 1:length(tmp)
                    [~,ixInd(iPct)] = min(abs(tmp(iPct)-PKrefValues));
                end
                
                writeToReportLog('INFO',sprintf('start sensitivity simulation for %s, %s %s',...
                    PopRunSet(iSet).name,OutputList(iO).reportName,OutputList(iO).pKParameterList{1,iPK}),false);

                
                % loop on Ind
                for iInd = 1:length(ixInd)
                    
                    if isempty(sens{iInd,iSet})
                    
                        % generate new "individuals" which differ by one parameter
                        [parPathsSens,parValuesSens,SensPointer] = generateSensitivityParameterSet(WSettings,...
                            sensParameterList,parPaths,parValues(ixInd(iInd),:));
                        nInd = size(parValuesSens,1);
                        
                        writeToReportLog('INFO',sprintf('start simulations for percentile %g individual ID %d with %d simulations',...
                            prctileSelection(iInd),ixInd(iInd),nInd),false);

                        % get result
                        SimResult = generateSimResult(WSettings,PopRunSet(iSet),nan,parPathsSens,parValuesSens,{OutputList(iO).pathID},...
                            nInd,ones(1,nInd).*individualVector(ixInd(iInd)));
                        
                        % start the PKParameter calculation
                        if ~isempty(PopRunSet(iSet).calculatePKParameterFh)
                            tmp = feval(PopRunSet(iSet).calculatePKParameterFh,WSettings,ApplicationProtocol,SimResult.time,SimResult.values{1},parPaths,parValues(ixInd(iInd),:));
                        else
%                             tmp = calculatePKParameterForApplicationProtocol(WSettings,ApplicationProtocol,...
%                                 SimResult.time,SimResult.values{1},weight(ixInd(iInd)),height(ixInd(iInd)),ixInd(iInd));
                            tmp = calculatePKParameterForApplicationProtocol(WSettings,ApplicationProtocol,...
                                  SimResult.time,SimResult.values{1},weight(ixInd(iInd)),height(ixInd(iInd)));
                        end
                        [jj,ix] = ismember(OutputList(iO).pKParameterList(1,iPK),{PKParameterTemplate.name});
                        PKPListSens = {tmp(ix(jj))};
                        
                        % construct referecne PK Parameter
                        PKPRef = PKPList{iO}(1,iPK);
                        PKPRef.value = PKrefValues(ixInd(iInd));
                        
                        tmp = calculateSensitivity(WSettings,{PKPRef},PKPListSens,SensPointer,1);
                        sens{iInd,iSet} =  tmp{1}; 
                    end
                    
                    save(sensFilename,'sens'); 
                end
            end
            
            % export results to csv
            analysisName = sprintf('O%d_%s',iO,PKParameterList(1,iPK).name);
            exportSensitivityForPopulation(WSettings,[analysisName '_detail'],sens(:,setList),sensParameterList,individualVector(ixInd),prctileSelection,{PopRunSet(setList).name});
            
            
            jj = strcmp(OutputList(iO).pKParameterList{1,iPK},{PKParameterTemplate.name});
            PKReportName = PKParameterTemplate(jj).reportName;
            
            header = sprintf('Sensitivity analysis of %s for %s',PKReportName,OutputList(iO).reportName);
            FP = FP.iniCaptiontextFigtextArray(header,analysisName);
            
            FP = plotSensitivityForPopulation(WSettings,FP,sens(:,setList),sensParameterList,prctileSelection,PopRunSet(setList),...
                OutputList(iO).reportName,PKReportName,analysisName);
            
            FP.saveCaptiontextArray;
        end
        
        
    end
    
    
    writeToReportLog('INFO','finalize sensitivity analysis',false);
  
catch exception
        
    save(sprintf('exception_%s.mat',datestr(now,'ddmmyy_hhMM')),'exception')
    writeToReportLog('ERROR',exception.message,false,exception);
    writeToReportLog('INFO',sprintf('sensitivity analysis finished with error \n'),false);
        
end  
    
    

return

function FP = plotSensitivityForPopulation(WSettings,FP,sensPop,sensParameterList,prctileSelection,PopRunSet,...
            OutputReportName,PKParameterReportName,analysisName)
        
% initialize color per percentile and population
nPrc=length(prctileSelection);
nPop=length(PopRunSet);
colorVectorPop= getcolmarkForMap(jet,max(nPop,nPrc));
% initialize marker per percentile
marker = 'ovs';
    
% construct percentile legend
lgtxtPrctl=cell(nPrc,1);
for iPrc=1:nPrc
    lgtxtPrctl{iPrc}=sprintf('%s percentile individual',getPercentilePotenzText(prctileSelection(iPrc)));
end

% get sorted sensitivities
[allSensValues,sortSumSensIx,iCut,CsortSumSensPop] = getListOfBestSensitivities(WSettings,sensPop);

% export overview
csv = createSensOverviewTable(sensParameterList,allSensValues,iCut,sortSumSensIx);
writeTabCellArray(csv,fullfile(FP.figureDir,[analysisName '_overview.csv']))

%Plot of the individual vs. mean sensitivities ----------------------
plotSensIndividualVsMean(WSettings,allSensValues,FP,marker,colorVectorPop,lgtxtPrctl,{PopRunSet.boxwhiskerLabel});

% get name and figure description
figureName = sprintf('%s_IndVsMean',analysisName);
figtxt = sprintf('Individual sensitivity vs mean sensitivity for %s for %s', PKParameterReportName,OutputReportName);
% save figure
FP = FP.printFigure(figureName,figtxt);

 %Plot of the sorted absolute sensitivity values ----------------------
FP = plotSensSortedValues(WSettings,allSensValues,sortSumSensIx,FP,marker,colorVectorPop,lgtxtPrctl,{PopRunSet.boxwhiskerLabel});

% get name and figure description
figureName = sprintf('%s_SortedAbsolute',analysisName);
figtxt = sprintf('Sorted absolute sensitivity for %s for %s', PKParameterReportName,OutputReportName);
% save figure
FP = FP.printFigure(figureName,figtxt);

%Plot of the sorted cumulated sensitivity values ---------------------------
FP = plotSensCumSortedValues(WSettings,CsortSumSensPop,FP,marker,colorVectorPop,lgtxtPrctl,{PopRunSet.boxwhiskerLabel},iCut);

% get name and figure description
figureName = sprintf('%s_CumSortedAbsolute',analysisName);
figtxt = sprintf('Cumulated sorted absolute sensitivity for %s for %s', PKParameterReportName,OutputReportName);
% save figure
FP = FP.printFigure(figureName,figtxt);
  
%Listing of the most sensitive parameters --------------------------- -----
plotSensListMostSensitive(WSettings,FP,sensPop,sortSumSensIx,iCut,...
    sensParameterList,marker,colorVectorPop,lgtxtPrctl,{PopRunSet.boxwhiskerLabel})

% get name and figure description
figureName = sprintf('%s_ListMostSensitive',analysisName);
figtxt = sprintf('Most sensitive parameter for %s for %s', PKParameterReportName,OutputReportName);
% save figure
FP = FP.printFigure(figureName,figtxt);

return
   
   

%-- plot individual vs mean
function  plotSensIndividualVsMean(WSettings,allSensValues,FP,marker,colorVectorPop,lgtxtPrctl,popLabels)


% create figure
ax = getReportFigure(WSettings,1,1,FP.figureHandle,'figureformat','portrait');

nPop = length(popLabels);
nPrc = length(lgtxtPrctl);

meanSens=mean(allSensValues(:,:,1),2);

% plot sensitivites
lgh = [];
for iPop=1:nPop
    for iPrc=1:nPrc
        
        if nPop==1
            col = colorVectorPop(iPrc,:);
        else
            col = colorVectorPop(iPop,:);
        end

        
        k=(iPop-1)*nPrc + iPrc;
        % plot sensitivites
        tmp = plot(meanSens,allSensValues(:,k,1),marker(iPrc),'color',col,'markerfacecolor',col);
        
        % set legend
        lgh(end+1) = tmp; %#ok<AGROW>
        set(tmp,'displayname',sprintf('%s %s',popLabels{iPop},lgtxtPrctl{iPrc}));
    end
end

% format figure
xlabel('mean sensitivity')
ylabel('individual sensitivity')
% set Axes scaling
xl = get(ax,'Xlim');
yl = get(ax,'Ylim');
xl = [-1 1]*max(abs([xl yl]));
setAxesScaling(ax,'xlim',xl,'xscale','lin',...
    'yscale','lin','ylim',xl);
grid on
set(ax,'dataaspectratio',[1 1 1])


% plot identity
plot(xl,xl,'k:');
% for all plot calculate one R^2
Y=reshape(allSensValues(:,:,1),nPop*nPrc*size(allSensValues(:,:,1),1),1);
n=length(Y);
X=[repmat(meanSens,nPrc*nPop,1),ones(n,1)];
[b,~,~,~,stats] = regress(Y,X);
lgh(end+1) = plot(xl,polyval(b,xl),'k-','displayname',sprintf('line of regression R^2=%.3g',stats(1)));

legend(lgh,get(lgh,'displayname'),'location','northoutside','fontsize',8)


return


%-- plots the sorted sensitivity
function FP = plotSensSortedValues(WSettings,allSensValues,sortSumSensIx,FP,marker,colorVectorPop,lgtxtPrctl,popLabels)
  
% create figure
ax = getReportFigure(WSettings,1,1,FP.figureHandle,'figureformat','portrait');

nPop = length(popLabels);
nPrc = length(lgtxtPrctl);

% plot sensitivity
yMin = min(min(abs(allSensValues(:,:,1))));
yMax = max(max(abs(allSensValues(:,:,1))));
lgh = [];
for iPop=1:nPop
    for iPrc=1:nPrc
        
        if nPop==1
            col = colorVectorPop(iPrc,:);
        else
            col = colorVectorPop(iPop,:);
        end
        
        k=(iPop-1)*nPrc + iPrc;
        % plot sensitivites
        tmp = plot(abs(allSensValues(sortSumSensIx,k,1)),marker(iPrc),'color',col,'markerfacecolor',col);
        
        % set legend
        lgh(end+1) = tmp; %#ok<AGROW>
        set(tmp,'displayname',sprintf('%s %s',popLabels{iPop},lgtxtPrctl{iPrc}));

                 
    end
end

% figureFormat
xlabel('parameters sorted after maximum absolute individual sensitivity')
%         xlabel('parameters sorted after absolute median sensitivity')
%         xlabel('parameters sorted after absolute mean ranking')
ylabel('absolute sensitivity')

% set axes scaling
yl = [yMin*0.9,yMax*1.1];
xl = [0 length(sortSumSensIx)+1];
setAxesScaling(ax,'xlim',xl,'xscale','lin',...
    'yscale','log','ylim',yl);
grid on
       
legend(lgh,get(lgh,'displayname'),'location','northoutside','fontsize',8)

return


% -- plot cumulated sorted sensitivity
function FP = plotSensCumSortedValues(WSettings,CsortSumSensPop,FP,marker,colorVectorPop,lgtxtPrctl,popLabels,iCut)


% create figure
ax = getReportFigure(WSettings,1,1,FP.figureHandle,'figureformat','portrait');

% get number of sensitivity vectors
nPop=size(CsortSumSensPop,2);
nPrc=size(CsortSumSensPop,1);

% plot sensitivity
yMin=inf;
yMax=0;
lgh = [];
for iPop=1:nPop
    for iPrc=1:nPrc

        if nPop==1
            col = colorVectorPop(iPrc,:);
        else
            col = colorVectorPop(iPop,:);
        end
        
        tmp = plot(CsortSumSensPop(iPrc,iPop).slope,marker(iPrc),'color',col,'markerfacecolor',col);
        hold on;
        
        % set legend
        lgh(end+1) = tmp; %#ok<AGROW>
        set(tmp,'displayname',sprintf('%s %s',popLabels{iPop},lgtxtPrctl{iPrc}));

        
        yMin=min(yMin,CsortSumSensPop(iPrc,iPop).slope(end));
        yMax=max(yMax,CsortSumSensPop(iPrc,iPop).slope(1));
    end
end

% figureFormat
xlabel('parameters sorted after maximum absolute individual sensitivity')
% xlabel('parameters sorted after absolute median sensitivity')
% xlabel('parameters sorted after absolute mean ranking')

ylabel('cumulated sensitivities')

% set axes scaling
yl = [yMin*0.9,yMax*1.1];
xl = [0 length(CsortSumSensPop(1,1).slope)+1];
setAxesScaling(ax,'xlim',xl,'xscale','lin',...
    'yscale','log','ylim',yl);
grid on

% plot Cut 
plot(xl,[1 1]*(1-WSettings.sensitivityCutoff),'k','linewidth',2);
yl = get(ax,'Ylim');
if ~isempty(iCut)
    lgh(end+1) = plot(iCut+0.5*[1 1],yl,'k','linewidth',2,'displayname','Cut off line');
end

legend(lgh,get(lgh,'displayname'),'location','northoutside','fontsize',8)

return


% generates the overview tables
function csv = createSensOverviewTable(parameterList,allSensValues,iCut,sortSumSensIx)

% calculates statistics:
mean_value=mean(allSensValues(sortSumSensIx,:,1),2);
std_value=std(allSensValues(sortSumSensIx,:,1),0,2);
max_value=max(abs(allSensValues(sortSumSensIx,:,1)),[],2);
     
% prepare xls export
csv(1,:) = {'Parameter name','sensitivity mean ± std (|max|)'};

for iPar=1:iCut
    
    csv(iPar+1,:) = {parameterList{sortSumSensIx(iPar),1},sprintf('%.3g±%.3g (%.3g)',mean_value(iPar),std_value(iPar),max_value(iPar))};
end

return
