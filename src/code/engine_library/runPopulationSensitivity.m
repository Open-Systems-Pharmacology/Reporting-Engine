function runPopulationSensitivity(WSettings,PopRunSet,workflowType)
% RUNPOPULATIONSENSITIVITY calculates and plots senstivity for all poulations
%
% runPopulationSensitivity(WSettings,PopRunSet,workflowType)
%
% Inputs
%       WSettings (structure)    definition of properties used in all
%                   workflow functions see GETDEFAULTWORKFLOWSETTINGS
%       PopRunSet (structure)   list of population simulations see GENERATEWORKFLOWINPUTFORPOPULATIONSIMULATION
%       workflowType (string) type of workflow

% Open Systems Pharmacology Suite;  http://open-systems-pharmacology.org

writeToLog('start sensitivity analysis',WSettings.logfile,true,false);

% percentiles for individual selection
prctileSelection = WSettings.sensitivitPrctileSelection;

if strcmp(workflowType,'pediatric')
    jjRef = [PopRunSet.isReference];
    setList = find(~jjRef);
else
    setList = 1:length(PopRunSet);
end

% create resultdirector
% Initialize figureDir
FP = ReportFigurePrint(fullfile('figures','sensitivity'),WSettings.printFormatList);

% get outputs and PKPList
load(fullfile('tmp',PopRunSet(setList(1)).name,'outputList.mat'),'OutputList');
% get corresponding application protocoll
load(fullfile('tmp',PopRunSet(setList(1)).name,'applicationProtocol.mat'),'ApplicationProtocol','isValid','PKParameterTemplate');


for iO = 1:length(OutputList)
        
    for iPK =  1:size(OutputList(iO).pKParameterList,2)

        for iSet = setList
            
            tmpDir = fullfile('tmp',PopRunSet(iSet).name);            
            PKPList = readPKParameter(PopRunSet(iSet).name);
            PKParameterList = PKPList{iO};

            % get corresponding population
            load(fullfile(tmpDir,'pop.mat'),'parPaths','parValues');
            jj = strcmp(parPaths,'IndividualId');
            individualVector = parValues(:,jj); %#ok<NODEF>
            jj = strcmp('Organism|Weight',parPaths);
            weight =  parValues(:,jj);
            jj = strcmp('Organism|Height',parPaths);
            height =  parValues(:,jj);


            % load sensitivity parameter
            load(fullfile(tmpDir,'senssitivity.mat'),'sensParameterList');
                   
            % get percentile individual
            tmp = prctile(PKParameterList(iPK).value,prctileSelection);  
            ixInd = nan(size(prctileSelection));
            for iPct = 1:length(tmp)
                [~,ixInd(iPct)] = min(abs(tmp(iPct)-PKParameterList(iPK).value));
            end
            
            % loop on Ind
            for iInd = 1:length(ixInd)
                % generate new "individuals" which differ by one parameter
                [parPathsSens,parValuesSens,SensPointer] = generateSensitivityParameterSet(WSettings,...
                    sensParameterList,parPaths,parValues(ixInd(iInd),:));
                nInd = size(parValuesSens,1);
                            
                % get result
                SimResult = generateSimResult(WSettings,PopRunSet(iSet),nan,parPathsSens,parValuesSens,{OutputList(iO).pathID},...
                    nInd,ones(1,nInd).*individualVector(ixInd(iInd))); 
                                
                % start the PKParameter calculation
                if ~isempty(PopRunSet(iSet).calculatePKParameterFh)
                    tmp = feval(PopRunSet(iSet).calculatePKParameterFh,WSettings,ApplicationProtocol,SimResult.time,SimResult.values{1},parPaths,parValues);            
                else    
                    tmp = calculatePKParameterForApplicationProtocol(WSettings,ApplicationProtocol,...
                        SimResult.time,SimResult.values{1},weight(ixInd(iInd)),height(ixInd(iInd)));
                end
                [jj,ix] = ismember(OutputList(iO).pKParameterList(1,:),{PKParameterTemplate.name});
                PKPListSens = {tmp(ix(jj))};

                tmp = calculateSensitivity(WSettings,PKPList(iO),PKPListSens,SensPointer,ixInd(iInd));
                sens{iInd,iSet} =  tmp{1,iPK}; %#ok<AGROW>
            end
        end
        
        % export results to csv
        analysisName = sprintf('%s_%s',removeForbiddenLetters(OutputList(iO).reportName),PKParameterList(iPK).name);
        exportSensitivityForPopulation(WSettings,[analysisName '_detail'],sens(:,setList),sensParameterList,individualVector(ixInd),prctileSelection,{PopRunSet(setList).name});
        
        
        jj = strcmp(OutputList(iO).pKParameterList{1,iPK},{PKParameterTemplate.name});
        PKReportName = PKParameterTemplate(jj).reportName;

        header = sprintf('sensitivity analysis of %s of %s',PKReportName,OutputList(iO).reportName);
        FP = FP.iniCaptiontextFigtextArray(header,analysisName);

        FP = plotSensitivityForPopulation(WSettings,FP,sens(:,setList),sensParameterList,prctileSelection,PopRunSet(setList),...
            OutputList(iO).reportName,PKReportName,analysisName);

        FP.saveCaptiontextArray;
    end
            
            
end


writeToLog('finalize sensitivity analysis',WSettings.logfile,true,false);

return

function FP = plotSensitivityForPopulation(WSettings,FP,sensPop,sensParameterList,prctileSelection,PopRunSet,...
            OutputReportName,PKParameterReportName,analysisName)
        
% initialize color per percentile and population
nPrc=length(prctileSelection);
nPop=length(PopRunSet);
colorVectorPop= getcolmarkForMap(jet,nPop);
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
figtxt = sprintf('sorted absloute sensitivity for %s for %s', PKParameterReportName,OutputReportName);
% save figure
FP = FP.printFigure(figureName,figtxt);

%Plot of the sorted cumulated sensitivity values ---------------------------
FP = plotSensCumSortedValues(WSettings,CsortSumSensPop,FP,marker,colorVectorPop,lgtxtPrctl,{PopRunSet.boxwhiskerLabel},iCut);

% get name and figure description
figureName = sprintf('%s_CumSortedAbsolute',analysisName);
figtxt = sprintf('cumulated sorted absloute sensitivity for %s for %s', PKParameterReportName,OutputReportName);
% save figure
FP = FP.printFigure(figureName,figtxt);
  
%Listing of the most sensitive parameters --------------------------- -----
plotSensListMostSensitive(WSettings,FP,sensPop,sortSumSensIx,iCut,...
    sensParameterList,marker,colorVectorPop,lgtxtPrctl,{PopRunSet.boxwhiskerLabel})

% get name and figure description
figureName = sprintf('%s_ListMostSensitive',analysisName);
figtxt = sprintf('Most sensitiv parameter for %s for %s', PKParameterReportName,OutputReportName);
% save figure
FP = FP.printFigure(figureName,figtxt);

return
   
   

%-- plot individual vs mean
function  plotSensIndividualVsMean(WSettings,allSensValues,FP,marker,colorVectorPop,lgtxtPrctl,popLabels)


% create figure
ax = getReportFigure(WSettings,1,1,FP.figureHandle,'figureformat','portrait');

nPop = size(colorVectorPop,1);
nPrc = length(lgtxtPrctl);

meanSens=mean(allSensValues(:,:,1),2);

% plot sensitivites
lgh = [];
for iPop=1:nPop
    lgtxt=cell(nPrc,1);
    for iPrc=1:nPrc
        
        k=(iPop-1)*nPrc + iPrc;
        % plot sensitivites
        tmp = plot(meanSens,allSensValues(:,k,1),marker(iPrc),'color',colorVectorPop(iPop,:),'markerfacecolor',colorVectorPop(iPop,:));
        
        % set legend
         lgh = addToLegendPopulationSensitivity(lgh,tmp,popLabels,lgtxtPrctl,iPop,iPrc);
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

nPop = size(colorVectorPop,1);
nPrc = length(lgtxtPrctl);

% plot sensitivity
yMin = min(min(abs(allSensValues(:,:,1))));
yMax = max(max(abs(allSensValues(:,:,1))));
lgh = [];
for iPop=1:nPop
    for iPrc=1:nPrc
        
        k=(iPop-1)*nPrc + iPrc;
        % plot sensitivites
        tmp = plot(abs(allSensValues(sortSumSensIx,k,1)),marker(iPrc),'color',colorVectorPop(iPop,:),'markerfacecolor',colorVectorPop(iPop,:));
        
        % set legend
         lgh = addToLegendPopulationSensitivity(lgh,tmp,popLabels,lgtxtPrctl,iPop,iPrc);
                 
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
nPop=size(CsortSumSensPop,1);
nPrc=size(CsortSumSensPop,2);

% plot sensitivity
yMin=inf;
yMax=0;
lgh = [];
for iPop=1:nPop
    for iPrc=1:nPrc

        tmp = plot(CsortSumSensPop(iPop,iPrc).slope,marker(iPrc),'color',colorVectorPop(iPop,:),'markerfacecolor',colorVectorPop(iPop,:));
        hold on;
        
        % set legend
         lgh = addToLegendPopulationSensitivity(lgh,tmp,popLabels,lgtxtPrctl,iPop,iPrc);
        
        yMin=min(yMin,CsortSumSensPop(iPop,iPrc).slope(end));
        yMax=max(yMax,CsortSumSensPop(iPop,iPrc).slope(1));
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
