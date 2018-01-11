function FP = plotLoopVPCTimeprofiles(WSettings,textFunctionHandle,Def,RunSet,FP) 
% PLOTLOOPVPCTIMEPROFILES generates timeprofile plots for a population VPC
%
% FP = plotLoopVPCTimeprofiles(WSettings,Def,RunSet,FP) 
%
% Inputs: 
%       Settings (structure)    definition of properties used in all
%                   workflow functions see GETDEFAULTWORKFLOWSETTINGS
%       textFunctionHandle (function handle) function handle to set figure text
%       Def (structure) contains information which plots should be
%               generated  see GETDEFAULTVPCPOPULATIONSETTINGS
%               (timeprofile)
%       RunSet (structure)   list of population simulations or MeanModel simulations 
%                see GENERATEWORKFLOWINPUTFORPOPULATIONSIMULATION
%                see GENERATEWORKFLOWINPUTFORMEANMODELSIMULATION
%       FP (ReportFigurPrint)  objects with manages print of figures and
%                           tables
%   Outputs:
%       FP (ReportFigurPrint)  objects with manages print of figures and

% Open Systems Pharmacology Suite;  http://open-systems-pharmacology.org

% get Reference simulation
if ~isempty(Def.ixRunSetRef)
    reportNameRef = RunSet(Def.ixRunSetRef).reportName;
    reportLabelRef = RunSet(Def.ixRunSetRef).popReportName;
else
    reportNameRef = '';
    reportLabelRef = '';
end


% get factor to translate time in  diplay units
% it may ba a cell array the for each Timerange a special unit is defined
if iscell(Def.timeDisplayUnit)
    timeUnitFactor = cellstr(@(x) getUnitFactor('',x,'Time'),Def.timeDisplayUnit);    
else    
    timeUnitFactor = getUnitFactor('',Def.timeDisplayUnit,'Time');
end

% set scaling for y axes
yscale = Def.yScale;

% start Loop on popsets
for iSet = [Def.ixRunSetRef Def.ixOfRunSets]
    
    
    %initialize new section per Pop runset
    header = sprintf('Time profiles for %s',RunSet(iSet).reportName);
    sheet = RunSet(iSet).name;
    FP = FP.iniCaptiontextFigtextArray(header,sheet);
    
    % get population specific names
    if isfield (RunSet,'popReportName')
        popReportName = RunSet(iSet).popReportName;
    else
        popReportName = '';
    end

    % get population shortname
    if isfield (RunSet,'boxwhiskerLabel')
        popLabelName = RunSet(iSet).boxwhiskerLabel;
    else
        popLabelName = popReportName;
    end

    
    % load OutputList of population
    OutputList = loadOutputList(RunSet(iSet).name);
    
    % load data if available
    [~,TP,~,dataReportName] = loadMergedData(WSettings,RunSet(iSet));
    if iscell(dataReportName)
        dataReportName = dataReportName{1};
    end
    
    
    % initialize residual vector for all outputs
    [residuals,residualsUnScaled,Jacob] = getResidualsForAllOutputs(WSettings,OutputList,RunSet(iSet),TP,1);

    % initialize goodness table for qc check
    goodness = {'output','plotType','messsage','success'};
    
    % Plot loop on Outputs
    for iO = 1:length(OutputList)
        
        %add new subcaption
        header = sprintf('Time profiles of %s',OutputList(iO).reportName);
        FP = FP.addSubSection(header,2);


        % load Simulation of main population, reference population and mean model
        Sim = loadAllSimulations(RunSet,iSet,iO,Def);
            
        
        if ~isempty(Jacob)
                % todo check calculation, and put outside loop of outputs
                VPCresult = getVPCRange(Sim.values,Sim.sensitivities,residualsUnScaled,Jacob,OutputList(iO).residualScale);
        else 
            VPCresult = [];
        end            

        
        % check if it is a multi application to generate timeLimits
        [timelimit,timeRangetxt,startTimes] = getTimeLimits(RunSet(iSet).name,Sim.time,Def);
        
        % loop on time ranges
        for iT = 1:size(timelimit,1)

            %add new subcaption
            if ~isempty(timeRangetxt{iT})
                FP = FP.addSubSection([upper(timeRangetxt{iT}(1)),timeRangetxt{iT}(2:end)],3);
            end

            % chekc if reference population shall be plotted
            doPlotReference = ~ismember(iSet,Def.ixRunSetRef) & ~isempty(Sim.valuesRef);
            reportNameRefTmp = reportNameRef;
            if ~doPlotReference
                reportNameRefTmp = '';
            end

            % get simulated time and y vectors within defined timerange and
            % convert to Units
            timeUnitFactorTL = timeUnitFactor(min(iT,length(timeUnitFactor)));
            [SimTL,VPCresultTL,timeLabel,startTimeTL] = getTimeVectorForTimeLimit(timelimit(iT,:),startTimes,...
                Sim,VPCresult,timeUnitFactorTL,OutputList(iO).unitFactor,doPlotReference);
            
            
            % get data for this timelimit           
            [DataTP,lloq] = prepareDataForTimeRange(WSettings,TP,iO,[startTimeTL timelimit(iT,2)],...
                SimTL.time./timeUnitFactorTL,SimTL.y, timeUnitFactorTL);
            
            % loop on scale
            for iScale = 1:length(yscale)
                
                %% timeprofile
                if ismember('yVsTime',Def.plotTypes) 

                    % take range plots for pop
                    if strcmp(WSettings.workflowType,'popModel') || isempty(VPCresult)
                    
                        % get name and figure description
                        figureName = sprintf('TP_%s_O%d_%s',RunSet(iSet).name,iO,yscale{iScale});
                        
                        [figtxt,figtxtTable,legendEntries] = feval(textFunctionHandle,WSettings,'tpShadedArea',...
                            {OutputList(iO).reportName,RunSet(iSet).reportName,reportNameRefTmp,dataReportName,timeRangetxt{iT},yscale{iScale},...
                            lloq,[size(SimTL.y,2),size(SimTL.yRef,2),length(DataTP)],...
                            popLabelName, reportLabelRef});
                       
                        
                        % do the plot
                        [csv,goodnessTmp] =  plotReportTimeProfile(WSettings,FP.figureHandle,SimTL,DataTP,timeLabel,Def.timeDisplayUnit,...
                            OutputList(iO).reportName,OutputList(iO).displayUnit,legendEntries,yscale{iScale},lloq);

                        if ~isnan(goodnessTmp{2}) && iScale==1
                            goodness(end+1,:) = [{OutputList(iO).reportName,'timeprofile with population range'}, goodnessTmp]; %#ok<AGROW>
                        end
                        
                        % save figure
                        if iScale==1
                            FP = FP.printFigure(figureName,figtxt,csv,figtxtTable);
                        else
                            FP = FP.printFigure(figureName,figtxt);
                        end

                    else
                        
                        % VPC parameter / data /  data and parameter uncertainty ------------------------------------------
                        VPCTypeTxt = { 'tpVPCpar','tpVPCdata','tpVPCdataPar'};
                        for iVPCFlag = 1:length(VPCTypeTxt)
                            % get name and figure description
                            figureName = sprintf('%s_%s_O%d_%s',VPCTypeTxt{iVPCFlag},RunSet(iSet).name,iO,yscale{iScale});
                            [figtxt,figtxtTable,legendEntries] = feval(textFunctionHandle,WSettings,VPCTypeTxt{iVPCFlag},...
                                {OutputList(iO).reportName,RunSet(iSet).reportName,dataReportName,yscale{iScale},...
                                lloq,length(DataTP),...
                                popReportName});
                        
                            
                            % do the plot
                            csv =  plotReportTimeProfile(WSettings,FP.figureHandle,SimTL,DataTP,timeLabel,Def.timeDisplayUnit,...
                                OutputList(iO).reportName,OutputList(iO).displayUnit,legendEntries,yscale{iScale},lloq,VPCresultTL,iVPCFlag);
                            
                            
                            
                            % save figure
                            if iScale==1
                                FP = FP.printFigure(figureName,figtxt,csv,figtxtTable);
                            else
                                FP = FP.printFigure(figureName,figtxt);
                            end
                        end
                        
                    end
                
                end
            end % iScale

            % comparison plots predicted vs observed
            if ~isempty(DataTP) && any(cellfun(@(x) any(~isnan(x)),{DataTP.predicted}))
                %% predicted vs observed
                if ismember('predVsObs',Def.plotTypes)
                    
                    % get name and figure description
                    figureName = sprintf('PvO_%s_O%d_%s',RunSet(iSet).name,iO,yscale{iScale});
                    [figtxt,~,legendEntries] = feval(textFunctionHandle,WSettings,'tpPredVsObs',...
                        {OutputList(iO).reportName,RunSet(iSet).reportName,dataReportName,yscale{iScale},...
                        lloq,popReportName});
                    
                    
                    % do the plot
                    plotReportPredictedVsObserved(WSettings,FP.figureHandle,DataTP,...
                        OutputList(iO).reportName,OutputList(iO).displayUnit,legendEntries,OutputList(iO).residualScale,lloq);
                    
                    % save figure
                    FP = FP.printFigure(figureName,figtxt);
                    
                end
                %% residuals vs time
                if ismember('resVsTime',Def.plotTypes)
                    
                    % get name and figure description
                    figureName = sprintf('RvT_%s_O%d_%s',RunSet(iSet).name,iO,yscale{iScale});
                    [figtxt,~,legendEntries] = feval(textFunctionHandle,WSettings,'tpResVsTime',...
                        {OutputList(iO).reportName,RunSet(iSet).reportName,dataReportName,yscale{iScale},...
                        lloq,popReportName});
                    
                    % do the plot
                    goodnessTmp = plotReportResiduals(WSettings,FP.figureHandle,DataTP,timeLabel,Def.timeDisplayUnit,...
                        OutputList(iO).reportName,OutputList(iO).displayUnit,OutputList(iO).residualScale,'vsTime',legendEntries);
                    
                    if ~isnan(goodnessTmp{2})
                        goodness(end+1,:) = [{OutputList(iO).reportName,'residuals vs time'}, goodnessTmp]; %#ok<AGROW>
                    end
                    
                    % save figure
                    FP = FP.printFigure(figureName,figtxt);
                end
                %% residuals vs y
                if ismember('resVsY',Def.plotTypes)
                    
                    % get name and figure description
                    figureName = sprintf('RvY_%s_O%d_%s',RunSet(iSet).name,iO,yscale{iScale});
                    [figtxt,~,legendEntries]  = feval(textFunctionHandle,WSettings,'tpResVsY',...
                        {OutputList(iO).reportName,RunSet(iSet).reportName,dataReportName,OutputList(iO).residualScale,...
                        lloq,popReportName});
                    
                    % do the plot
                    goodnessTmp = plotReportResiduals(WSettings,FP.figureHandle,DataTP,timeLabel,Def.timeDisplayUnit,...
                        OutputList(iO).reportName,OutputList(iO).displayUnit,OutputList(iO).residualScale,'vsY',legendEntries);
                    
                    if ~isnan(goodnessTmp{2})
                        goodness(end+1,:) = [{OutputList(iO).reportName,'residuals vs Y'}, goodnessTmp]; %#ok<AGROW>
                    end
                end

                
                % save figure
                FP = FP.printFigure(figureName,figtxt);
            end
            
            
            
        end % iT

        
            
    end % iO
    
    
    % histogramm of all Residuals
    if ~isempty(residuals) && any(ismember({'histRes','qqPlotRes'},Def.plotTypes)) 

        FP = FP.addSubSection('Overview on all residuals',2);

    
    
        % histogramm of all Residuals
        if ismember({'histRes'},Def.plotTypes)   
            
            % get name and figure description
            figureName = 'ResidualHist';
            [figtxt,~,legendEntries] = feval(textFunctionHandle,WSettings,'histRes',...
                {{OutputList.reportName},RunSet(iSet).reportName,popReportName});
            
            % do the plot
            plotReportHistogram(WSettings,FP.figureHandle,residuals(:,1) ,residuals(:,2),[],[],...
                'Residuals','',{},legendEntries,'isResiduals');
            
            % save figure
            FP = FP.printFigure(figureName,figtxt);
            
        end
        
        % histogramm of all Residuals
        if ismember({'qqPlotRes'},Def.plotTypes)  
            
            % get name and figure description
            figureName = 'ResidualQQplot';
            figtxt = feval(textFunctionHandle,WSettings,'qqRes',...
                {{OutputList.reportName},RunSet(iSet).reportName,popReportName});
            
            % do the plot
            plotReportQQPlot(WSettings,FP.figureHandle,residuals(:,1));
            
            % save figure
            FP = FP.printFigure(figureName,figtxt);
            
        end
    end
    
    FP.saveCaptiontextArray;
            
            

    
    % write goodness
    if size(goodness,1)>1
        fname = fullfile(FP.figureDir,sprintf('goodnessOfFit_%s.csv',RunSet(iSet).name));
        writeTabCellArray(goodness,fname);
    end

end % iSet export


return

function  [DataTP,lloq] = prepareDataForTimeRange(WSettings,TP,iO,timelimit,time,y,timeUnitFactor)
% all input time variables, timlimit, time and TP.time are given in min
% DataTP is returned in display units


if size(y,2)>1
    [~,y] = getRangePlotPercentiles(WSettings,y');
end


% initialize return values
DataTP = [];
lloq = nan;

% if no data are avilableeeee nothing to do
if isempty(TP) || isempty(TP{iO})
    return
end


%  get timefield
if isfield(TP{iO},'tad') && timelimit(1,1)>0
    timefield = 'tad';
    timeoffset =0;    
else
    timefield = 'time';
    timeoffset = timelimit(1,1);        
end
    
% get stucture to plot
for iInd = 1:length(TP{iO})
        
    jjT = TP{iO}(iInd).time >= timelimit(1,1) & ...
        TP{iO}(iInd).time <= timelimit(1,2);

    if any(jjT)
        DataTP(iInd).time = (TP{iO}(iInd).(timefield)(jjT)-timeoffset).*timeUnitFactor; %#ok<AGROW>
        DataTP(iInd).y(:,1) = TP{iO}(iInd).dv(jjT); %#ok<AGROW>
        
        % set lloq
        isLloq = TP{iO}(iInd).isLloq(jjT);
        DataTP(iInd).lloq(:,1) = nan(length(isLloq),1); %#ok<AGROW>
        if any(isLloq)
            DataTP(iInd).y(isLloq,1) = nan; %#ok<AGROW>
            lloq = [TP{iO}(iInd).lloq];
            DataTP(iInd).lloq(isLloq,1) = lloq /2; %#ok<AGROW>
        end
        
        % get predicted
        DataTP(iInd).predicted(:,1) = interp1(time.*timeUnitFactor,y,DataTP(iInd).time); %#ok<AGROW>
    else
        DataTP(iInd).time = []; %#ok<AGROW>
        DataTP(iInd).y = [];  %#ok<AGROW>
        DataTP(iInd).predicted = [];  %#ok<AGROW>
        DataTP(iInd).lloq = [];  %#ok<AGROW>
    end
        
end

return



function [residuals,residualsUnScaled,Jacob] = getResidualsForAllOutputs(WSettings,OutputList,RunSet,TP,timeUnitFactor)
    
residuals = [];
residualsUnScaled = [];
Jacob = [];
    
% Get residuals and Jacobian
for iO = 1:length(OutputList)
    
    % load Simulation
    [simTime,simValues,~,~,~,sensitivities] = loadSimResult(RunSet.name,iO);
    
    % prepare residualOveralPlot
    if  ~isempty(TP)
        if ~isempty(TP{iO})
            
            % get Y values for time range
            y = simValues.*OutputList(iO).unitFactor;
            % get data for this timelimit
            [DataTP] = prepareDataForTimeRange(WSettings,TP,iO,[0 simTime(end)],simTime,y,timeUnitFactor);
            
            ixRes = [];
            
            
            for iInd = 1:length(DataTP)
                
                offset = size(residuals,1);
                
                switch OutputList(iO).residualScale
                    case 'lin'
                        jj = ~isnan(DataTP(iInd).y) & ~isnan(DataTP(iInd).predicted);
                        residualsUnScaled(offset+[1:sum(jj)],1) =  (DataTP(iInd).y(jj)-DataTP(iInd).predicted(jj))./OutputList(iO).unitFactor; %#ok<NBRAK,AGROW>
                        residuals(offset+[1:sum(jj)],1) =  (DataTP(iInd).y(jj)-DataTP(iInd).predicted(jj)); %#ok<NBRAK,AGROW>
                    case 'log'
                        jj = DataTP(iInd).y > 0 & DataTP(iInd).predicted > 0;
                        residualsUnScaled(offset+[1:sum(jj)],1) =  log(DataTP(iInd).y(jj))-log(DataTP(iInd).predicted(jj)); %#ok<NBRAK,AGROW>
                        residuals(offset+[1:sum(jj)],1) =  log(DataTP(iInd).y(jj))-log(DataTP(iInd).predicted(jj)); %#ok<NBRAK,AGROW>
                    otherwise
                        error('scale');
                        
                end
                residuals(offset+[1:sum(jj)],2) = iO; %#ok<NBRAK,AGROW>
                ixRes(offset+[1:sum(jj)]) = round(interp1(simTime.*timeUnitFactor,1:length(simTime),DataTP(iInd).time(jj)));%#ok<NBRAK,AGROW>
            end
            
            
            if ~isempty(sensitivities)
                Jacob = [Jacob,sensitivities(ixRes,:)]; %#ok<AGROW>
            end
        end
        
    end
end

return



 % load Simulation of main population, reference population and mean model
function Sim = loadAllSimulations(RunSet,iSet,iO,Def)
   
[Sim.time,Sim.values,pathID,~,~,Sim.sensitivities] = loadSimResult(RunSet(iSet).name,iO);
        
% load SImulation of reference simulation
if isfield(Def,'ixRunSetRef') &&  ~isempty(Def.ixRunSetRef)
    [Sim.timeRef,Sim.valuesRef] = loadSimResult(RunSet(Def.ixRunSetRef).name,nan,pathID);
else
    Sim.timeRef = [];
    Sim.valuesRef = [];
end

% load mean model        
if isfield(Def,'plotMeanModel') &&  Def.plotMeanModel
    [~,Sim.valuesMeanModel] = loadSimResult(RunSet(iSet).name,iO,'','meanModel_');
else
    Sim.valuesMeanModel = [];
end

return


function [timelimit,timeRangetxt,startTimes] = getTimeLimits(RunSetName,simTime,Def)
  
% get start times of application
[ApplicationProtocol,isValid] = loadApplicationProtocoll(RunSetName);
if isValid
    startTimes = unique([ApplicationProtocol.startTime]);
else
    startTimes = simTime(1);
end

% user defined timelimits
if ~isempty(Def.timelimit)
    timelimit = Def.timelimit;
            
    % take only ranges where endpoint of simulationrange is larger than beginning of timelimit range
    jj = timelimit(:,1) < simTime(end);
    timelimit = timelimit(jj,:);
    for ij = find(~jj)'
        writeToReportLog('WARNING',sprintf('Time Range %d - %d outside simulation range %d %d. Will be skipped.',...
            timelimit(ij,1),timelimit(ij,2),simTime(1),simTime(end)),false);
    end
    % shorten ranges where endpoint of simulationrange is less than end of timelimit range
    jj = timelimit(:,2) > simTime(end);
    for ij = find(jj)'
        writeToReportLog('WARNING',sprintf('Time Range %d - %d partly outside simulation range %d %d. Will be shortened.',...
            timelimit(ij,1),timelimit(ij,2),simTime(1),simTime(end)),false);
        timelimit(ij,2) = simTime(end);
    end
    
    % set text
    if ~isfield(Def,'timeRangetxt')
        for iT = 1:size(timelimit,1)
            timeRangetxt{iT} = sprintf('Simulation time range %g - %g %s',timelimit(iT,1).*timeUnitFactor,...
                timelimit(iT,2).*timeUnitFactor,Def.timeDisplayUnit); %#ok<AGROW>
        end
    else
        timeRangetxt = Def.timeRangetxt;
    end
% for invalid applications or single applications take the total simulation range    
elseif ~isValid || length(startTimes) ==1
    timelimit = [simTime(1) simTime(end)];
    timeRangetxt = {''};
% for multi applciations takt the total, the first and the last applciation range    
else
    timelimit = [simTime([1 end]);...
        startTimes(1:2);...
        startTimes(end) simTime(end)];
    timeRangetxt = {'for total simulation time range','for first application range','for last application range'};
end

return


function [SimTL,VPCresultTL,timeLabel,startTimeTL] = getTimeVectorForTimeLimit(timelimit,startTimes,Sim,VPCresult,...
    timeUnitFactor,OutputListUnitFactor,doPlotReference)
            
% find first dose in timelimit
idxNextDose = find(startTimes >= timelimit(1,1),1);
startTimeTL = startTimes(idxNextDose);

% get indices for time range and adjust time 
jjT = Sim.time >= timelimit(1,1) & Sim.time <= timelimit(1,2);
SimTL.time = (Sim.time(jjT) - startTimeTL).*timeUnitFactor;

% get Y values for time range
SimTL.y = Sim.values(jjT,:).*OutputListUnitFactor;


% get indices for time range and adjust time and y  for reference simulation
if ~isempty(Sim.timeRef) && doPlotReference
    jjTRef = Sim.timeRef >= timelimit(1,1) & Sim.timeRef <= timelimit(1,2);
    SimTL.timeRef = (Sim.timeRef(jjTRef) - startTimeTL).*timeUnitFactor;
    SimTL.yRef = Sim.valuesRef(jjTRef,:).*OutputListUnitFactor;
else
    SimTL.timeRef = [];
    SimTL.yRef = [];
end

% check if Mean Model Result is loaded
SimTL.yMeanModel = [];
if ~isempty(Sim.valuesMeanModel)
    SimTL.yMeanModel = Sim.valuesMeanModel(jjT,:).*OutputListUnitFactor;
end

% check if VPC output exists
if ~isempty(VPCresult)
    VPCresultTL = VPCresult;
    VPCresultTL.yMin = VPCresultTL.yMin(jjT,:).*OutputListUnitFactor;
    VPCresultTL.yMax = VPCresultTL.yMax(jjT,:).*OutputListUnitFactor;
else
    VPCresultTL = [];
end

% switch tad or Time
if startTimes(idxNextDose)>0
    timeLabel = 'time after dose';
else
    timeLabel = 'time';
end

