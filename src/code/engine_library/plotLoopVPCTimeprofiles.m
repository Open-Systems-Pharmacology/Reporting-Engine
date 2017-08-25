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
if isfield(Def,'ixRunSetRef') &&  ~isempty(Def.ixRunSetRef)
    popReportNameRef = RunSet(Def.ixRunSetRef).popReportName;
    reportNameRef = RunSet(Def.ixRunSetRef).reportName;
else
    popReportNameRef = '';
    reportNameRef = '';
end


% get factor to translate time in  diplay units
timeUnitFactor = getUnitFactor('',Def.timeDisplayUnit,'Time');

% set scaling for y axes
yscale = Def.yScale;

% start Loop on popsets
for iSet = Def.ixOfRunSets

    
    % get population specific names
    if isfield (RunSet,'popReportName')
        popReportName = RunSet(iSet).popReportName;
    else
        popReportName = '';
    end
    
    % load OutputList of population
    load(fullfile('tmp',RunSet(iSet).name,'outputList.mat'));
    
    % load data if available
    [~,TP,Dict] = loadMergedData(WSettings,{RunSet(iSet).name});
    dataReportName = RunSet(iSet).dataReportName;
    
    % get Unitfactor factor for data
    if ~isempty(Dict)
        jj = strcmp({Dict.matlabID},'time');
        dataTimeUnitFactor = getUnitFactor(Dict(jj).nonmemUnit,Def.timeDisplayUnit,'Time');
    else
        dataTimeUnitFactor = nan;
    end
    
    
    % initialize residual vector for all outputs
    residuals = [];
    residualsUnScaled = [];
    Jacob = [];
    goodness = {'output','plotType','messsage','success'};
    
    % Get residuals and VPC ranges
    for iO = 1:length(OutputList)

        % load Simulation
        [simTime,simValues,~,~,~,sensitivities] = loadSimResult(RunSet(iSet).name,iO);
        
        % prepare residualOveralPlot
        if  ~isempty(TP) 
            if ~isempty(TP{iO})
        
                % get Y values for time range
                y = simValues.*OutputList(iO).unitFactor;
                % get data for this timelimit
                [DataTP] = prepareDataForTimeRange(WSettings,TP,iO,[0 simTime(end)],simTime.*timeUnitFactor,y,dataTimeUnitFactor);
                
                ixRes = [];
                
                
                for iInd = 1:length(DataTP)
                    
                    offset = size(residuals,1);
                    
                    switch OutputList(iO).residualScale
                        case 'lin'
                            jj = ~isnan(DataTP(iInd).y) & ~isnan(DataTP(iInd).predicted);
                            residualsUnScaled(offset+[1:sum(jj)],1) =  (DataTP(iInd).y(jj)-DataTP(iInd).predicted(jj))./OutputList(iO).unitFactor; %#ok<AGROW>
                            residuals(offset+[1:sum(jj)],1) =  (DataTP(iInd).y(jj)-DataTP(iInd).predicted(jj)); %#ok<AGROW>
                        case 'log'
                            jj = DataTP(iInd).y > 0 & DataTP(iInd).predicted > 0;
                            residualsUnScaled(offset+[1:sum(jj)],1) =  log(DataTP(iInd).y(jj))-log(DataTP(iInd).predicted(jj)); %#ok<AGROW>
                            residuals(offset+[1:sum(jj)],1) =  log(DataTP(iInd).y(jj))-log(DataTP(iInd).predicted(jj)); %#ok<AGROW>
                        otherwise
                            error('scale');
                            
                    end
                    residuals(offset+[1:sum(jj)],2) = iO; %#ok<AGROW>
                    ixRes(offset+[1:sum(jj)]) = round(interp1(simTime.*timeUnitFactor,1:length(simTime),DataTP(iInd).time(jj)));%#ok<AGROW>
                end
                
                
                if ~isempty(sensitivities)
                    Jacob = [Jacob,sensitivities(ixRes,:)]; %#ok<AGROW>
                end
            end
            
        end
    end

    

    
    % Plot loop on Outputs
    for iO = 1:length(OutputList)

        % load Simulation
        [simTime,simValues,pathID,~,~,sensitivities] = loadSimResult(RunSet(iSet).name,iO);
        if isfield(Def,'ixRunSetRef') &&  ~isempty(Def.ixRunSetRef)
            [simTimeRef,simValuesRef] = loadSimResult(RunSet(Def.ixRunSetRef).name,nan,pathID);   
        else
            simTimeRef = [];
            simValuesRef = [];
        end
        
        if ~isempty(Jacob)
                % todo check calculation, and put outside loop of outputs
                VPCresult = getVPCRange(simValues,sensitivities,residualsUnScaled,Jacob,OutputList(iO).residualScale);
        else 
            VPCresult = [];
        end            

        
        % check if it is a multi application to generate timeLimits
        load(fullfile('tmp',RunSet(iSet).name,'applicationProtocol.mat'));
        if isValid
            startTimes = unique([ApplicationProtocol.startTime]);
        end
        
        if ~isempty(Def.timelimit)
            timelimit = Def.timelimit;
            timeRangetxt = sprintf('simulation time range %g - %g %s',timelimit(1),timelimit(2),Def.timeDisplayUnit);
        elseif ~isValid || length(startTimes) ==1
            timelimit = [simTime(1) simTime(end)];
            timeRangetxt = {};
        else
            timelimit = [simTime([1 end])';...
                startTimes(1:2);...
                startTimes(end) simTime(end)];
            timeRangetxt = {'total simulation time range','first application range','last application range'};
        end
        
    
        
        for iT = 1:size(timelimit,1)

            % get indices for time range and adjust time
            jjT = simTime >= timelimit(iT,1) & simTime <= timelimit(iT,2);
            time = (simTime(jjT) - timelimit(iT,1)).*timeUnitFactor;

            % get indices for time range and adjust time for reference simulation
            if ~isempty(simTimeRef)
                jjTRef = simTimeRef >= timelimit(iT,1) & simTimeRef <= timelimit(iT,2);
                timeRef = (simTimeRef(jjT) - timelimit(iT,1)).*timeUnitFactor;
            else
                timeRef = [];
            end
            
            % switch tad or Time
            if timelimit(iT,1)>0;
                timeLabel = 'Time after dose';
            else
                timeLabel = 'Time';
            end
            
            
            %initialize new sheet
            header = sprintf('Time profiles of %s for a %s',OutputList(iO).reportName,RunSet(iSet).reportName);
            sheet = RunSet(iSet).name;
            if ~isempty(timeRangetxt)
                header = sprintf('%s of the %s',header,timeRangetxt{iT});
                sheet = sprintf('%s_%s_%s',sheet,removeForbiddenLetters(OutputList(iO).reportName),removeForbiddenLetters(timeRangetxt{iT}));
            end
            FP = FP.iniCaptiontextFigtextArray(header,sheet);

            
        
            % get Y values for time range
            y = simValues(jjT,:).*OutputList(iO).unitFactor;
            
            % check if  reference output exists
            yRef = [];
            if ~isempty(simValuesRef)
                yRef = simValuesRef(jjTRef,:).*OutputList(iO).unitFactor;
            end

            % check if VPC output exists
            if ~isempty(VPCresult)
                VPCresult_time = VPCresult;
                VPCresult_time.yMin = VPCresult_time.yMin(jjT,:).*OutputList(iO).unitFactor;
                VPCresult_time.yMax = VPCresult_time.yMax(jjT,:).*OutputList(iO).unitFactor;
            end

            
            % get data for this timelimit           
            [DataTP,lloq] = prepareDataForTimeRange(WSettings,TP,iO,timelimit(1,:),time,y,dataTimeUnitFactor);
            
            % loop on scale
            for iScale = 1:length(yscale)
                
                %% timeprofile
                if ismember('yVsTime',Def.plotTypes) 

                    % take range plots for pop
                    if strcmp(WSettings.workflowType,'popModel') || isempty(VPCresult)
                    
                        % get name and figure description
                        figureName = sprintf('O%d_TP_%s_%s',iO,removeForbiddenLetters(OutputList(iO).reportName),yscale{iScale});
                        [figtxt,figtxtTable,legendEntries] = feval(textFunctionHandle,WSettings,'tpShadedArea',...
                            {OutputList(iO).reportName,RunSet(iSet).reportName,dataReportName,yscale{iScale},...
                            lloq,length(DataTP),...
                            popReportName,reportNameRef, popReportNameRef});
                        
                        
                        % do the plot
                        [csv,goodnessTmp] =  plotReportTimeProfile(WSettings,FP.figureHandle,time,y,timeRef,yRef,DataTP,timeLabel,Def.timeDisplayUnit,...
                            OutputList(iO).reportName,OutputList(iO).displayUnit,legendEntries,yscale{iScale},lloq);

                        if ~isnan(goodnessTmp{2}) && iScale==1
                            goodness(end+1,:) = [{OutputList(iO).reportName,'timeprofile with poulation range'}, goodnessTmp]; %#ok<AGROW>
                        end
                        
                        % save figure
                        if iScale==1
                            FP = FP.printFigure(figureName,figtxt,csv,figtxtTable);
                        else
                            FP = FP.printFigure(figureName,figtxt);
                        end

                    else
                        
                        % VPC parameter uncertainty ------------------------------------------
                        % get name and figure description
                        figureName = sprintf('O%d_TPVPCpar_%s_%s',iO,removeForbiddenLetters(OutputList(iO).reportName),yscale{iScale});
                        [figtxt,figtxtTable,legendEntries] = feval(textFunctionHandle,WSettings,'tpVPCpar',...
                            {OutputList(iO).reportName,RunSet(iSet).reportName,dataReportName,yscale{iScale},...
                            lloq,length(DataTP),...
                            popReportName});
                        
                        
                        % do the plot
                        csv =  plotReportTimeProfile(WSettings,FP.figureHandle,time,y,timeRef,yRef,DataTP,timeLabel,Def.timeDisplayUnit,...
                            OutputList(iO).reportName,OutputList(iO).displayUnit,legendEntries,yscale{iScale},lloq,VPCresult_time,1);


                        
                        % save figure
                        if iScale==1
                            FP = FP.printFigure(figureName,figtxt,csv,figtxtTable);
                        else
                            FP = FP.printFigure(figureName,figtxt);
                        end
                        
                        % VPC data uncertainty ------------------------------------------
                        % get name and figure description
                        figureName = sprintf('O%d_TPVPCdata_%s_%s',iO,removeForbiddenLetters(OutputList(iO).reportName),yscale{iScale});
                        [figtxt,figtxtTable,legendEntries] = feval(textFunctionHandle,WSettings,'tpVPCdata',...
                            {OutputList(iO).reportName,RunSet(iSet).reportName,dataReportName,yscale{iScale},...
                            lloq,length(DataTP),...
                            popReportName});
                        
                        
                        % do the plot
                        csv =  plotReportTimeProfile(WSettings,FP.figureHandle,time,y,timeRef,yRef,DataTP,timeLabel,Def.timeDisplayUnit,...
                            OutputList(iO).reportName,OutputList(iO).displayUnit,legendEntries,yscale{iScale},lloq,VPCresult_time,2);

                        % save figure
                        if iScale==1
                            FP = FP.printFigure(figureName,figtxt,csv,figtxtTable);
                        else
                            FP = FP.printFigure(figureName,figtxt);
                        end

                        % VPC data and parameter uncertainty ------------------------------------------
                        % get name and figure description
                        figureName = sprintf('O%d_TPVPCdataPar_%s_%s',iO,removeForbiddenLetters(OutputList(iO).reportName),yscale{iScale});
                        [figtxt,figtxtTable,legendEntries] = feval(textFunctionHandle,WSettings,'tpVPCdataPar',...
                            {OutputList(iO).reportName,RunSet(iSet).reportName,dataReportName,yscale{iScale},...
                            lloq,length(DataTP),...
                            popReportName});
                        
                        
                        % do the plot
                        [csv,goodnessTmp] =  plotReportTimeProfile(WSettings,FP.figureHandle,time,y,timeRef,yRef,DataTP,timeLabel,Def.timeDisplayUnit,...
                            OutputList(iO).reportName,OutputList(iO).displayUnit,legendEntries,yscale{iScale},lloq,VPCresult_time,3);
                        
                        if ~isnan(goodnessTmp{2}) && iScale==1
                            goodness(end+1,:) = [{OutputList(iO).reportName,'timeprofile with uncertainty'}, goodnessTmp]; %#ok<AGROW>
                        end

                        
                        % save figure
                        if iScale==1
                            FP = FP.printFigure(figureName,figtxt,csv,figtxtTable);
                        else
                            FP = FP.printFigure(figureName,figtxt);
                        end
                        
                    end
                
                end
            end % iScale

            %% predicted vs observed
            if ismember('predVsObs',Def.plotTypes) && ~isempty(DataTP)
                
                % get name and figure description
                figureName = sprintf('O%d_PvO_%s_%s',iO,removeForbiddenLetters(OutputList(iO).reportName),OutputList(iO).residualScale);
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
            if ismember('resVsTime',Def.plotTypes)  && ~isempty(DataTP)
                
                % get name and figure description
                figureName = sprintf('O%d_RvT_%s_%s',iO,removeForbiddenLetters(OutputList(iO).reportName),yscale{iScale});
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
            if ismember('resVsY',Def.plotTypes)  && ~isempty(DataTP)
                
                % get name and figure description
                figureName = sprintf('O%d_RvY_%s_%s',iO,removeForbiddenLetters(OutputList(iO).reportName),yscale{iScale});
                [figtxt,~,legendEntries]  = feval(textFunctionHandle,WSettings,'tpResVsY',...
                    {OutputList(iO).reportName,RunSet(iSet).reportName,dataReportName,OutputList(iO).residualScale,...
                    lloq,popReportName});
                
                % do the plot
                goodnessTmp = plotReportResiduals(WSettings,FP.figureHandle,DataTP,timeLabel,Def.timeDisplayUnit,...
                    OutputList(iO).reportName,OutputList(iO).displayUnit,OutputList(iO).residualScale,'vsY',legendEntries);
                
                if ~isnan(goodnessTmp{2})
                    goodness(end+1,:) = [{OutputList(iO).reportName,'residuals vs Y'}, goodnessTmp]; %#ok<AGROW>
                end

                
                % save figure
                FP = FP.printFigure(figureName,figtxt);
            end
            
            
            FP.saveCaptiontextArray;
            
            
            
        end % iT

        
            
    end % iO
    
    % histogramm of all Residuals
    if ismember({'histRes'},Def.plotTypes)   && ~isempty(residuals)
        
        % get name and figure description
        figureName = 'ResidualHist';
        [figtxt,~,legendEntries] = feval(textFunctionHandle,WSettings,'histRes',...
            {{OutputList.reportName},RunSet(iSet).reportName,popReportName});
        
        % do the plot
        plotReportHistogram(WSettings,FP.figureHandle,residuals(:,1) ,residuals(:,2),[],...
            'Residuals','',{},legendEntries,'isResiduals');
        
        % save figure
        FP = FP.printFigure(figureName,figtxt);
        
    end
    
    % histogramm of all Residuals
    if ismember({'qqPlotRes'},Def.plotTypes)   && ~isempty(residuals)
        
        % get name and figure description
        figureName = 'ResidualQQplot';
        figtxt = feval(textFunctionHandle,WSettings,'qqRes',...
            {{OutputList.reportName},RunSet(iSet).reportName,popReportName});
        
        % do the plot
        plotReportQQPlot(WSettings,FP.figureHandle,residuals(:,1));
        
        % save figure
        FP = FP.printFigure(figureName,figtxt);
        
    end
    
    % write goodness
    if size(goodness,1)>1
        writeGoodness(goodness,RunSet(iSet).name,FP.figureDir);
    end

end % iSet export


return
    

function writeGoodness(goodness,simulationName,figureDir)

data(1,:) = goodness(1,:);
data(2,:) = {'string','string','string','double'};
for iCol = 1:size(data,2)
    switch  data{2,iCol}
        case 'string'
            data{3,iCol} =goodness(2:end,iCol);
        case 'double'
            data{3,iCol} = double(cell2mat(goodness(2:end,iCol)));
        otherwise
            error('unknown datatype')
    end
end
fname = fullfile(figureDir,sprintf('goodnessOfFit_%s.csv',simulationName));
writetab(fname, data, ';', 0, 0, 1,0);

return

function  [DataTP,lloq] = prepareDataForTimeRange(WSettings,TP,iO,timelimit,time,y,dataTimeUnitFactor)


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
        
    jjT = TP{iO}(iInd).time.*dataTimeUnitFactor >= timelimit(1,1) & ...
        TP{iO}(iInd).time.*dataTimeUnitFactor <= timelimit(1,2);

    DataTP(iInd).time = TP{iO}(iInd).(timefield).*dataTimeUnitFactor-timeoffset; %#ok<AGROW>
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
    DataTP(iInd).predicted(:,1) = interp1(time,y,DataTP(iInd).time); %#ok<AGROW>
        
end

return
