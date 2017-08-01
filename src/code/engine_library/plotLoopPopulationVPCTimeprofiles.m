function FP = plotLoopPopulationVPCTimeprofiles(WSettings,textFunctionHandle,Def,PopRunSet,FP) 
% PLOTLOOPPOPULATIONVPCTIMEPROFILES generates timeprofile plots for a population VPC
%
% FP = plotLoopPopulationVPCTimeprofiles(WSettings,Def,PopRunSet,FP) 
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

% Open Systems Pharmacology Suite;  http://open-systems-pharmacology.org


% get Reference simulation
if ~isempty(Def.ixPopRunSetRef)
    SimResultRef = loadSim(PopRunSet(Def.ixPopRunSetRef).name);
    popReportNameRef = PopRunSet(Def.ixPopRunSetRef).popReportName;
    reportNameRef = PopRunSet(Def.ixPopRunSetRef).reportName;
else
    SimResultRef = [];
    popReportNameRef = '';
    reportNameRef = '';
end

% get factor to translate time in  diplay units
timeUnitFactor = getUnitFactor('',Def.timeDisplayUnit,'Time');

% set scaling for y axes
yscale = {'lin','log'};

% start Loop on popsets
for iSet = Def.ixOfPopRunSets


    % load Simulation
    SimResult = loadSim(PopRunSet(iSet).name);
    
    % check if it is a multi application to generate timeLimits
    load(fullfile('tmp',PopRunSet(iSet).name,'applicationProtocol.mat'));
    if isValid
        startTimes = unique([ApplicationProtocol.startTime]);
    end
        
    if ~isempty(Def.timelimit)
        timelimit = Def.timelimit;
        timeRangetxt = sprintf('simulation time range %g - %g %s',timelimit(1),timelimit(2),Def.timeDisplayUnit);
    elseif ~isValid || length(startTimes) ==1
        timelimit = SimResult.time([1 end])';
        timeRangetxt = {};
    else
        timelimit = [SimResult.time([1 end])';...
            startTimes(1:2);...
            startTimes(end) SimResult.time(end)];
        timeRangetxt = {'total simulation time range','first application range','last application range'};
    end
    
    % load OutputList of population
    load(fullfile('tmp',PopRunSet(iSet).name,'outputList.mat'));
    
    % load data if available
    [~,TP,~] = loadMergedData(WSettings,{PopRunSet(iSet).name});
    dataReportName = PopRunSet(iSet).dataReportName;
   
    
    for iT = 1:size(timelimit,1)
        
        %initialize new sheet
        header = sprintf('Time profiles of %s',PopRunSet(iSet).reportName);
        sheet = PopRunSet(iSet).name;
        if ~isempty(timeRangetxt)
            header = sprintf('%s of the %s',header,timeRangetxt{iT});
            sheet = sprintf('%s_%s',sheet,removeForbiddenLetters(timeRangetxt{iT}));
        end
        FP = FP.iniCaptiontextFigtextArray(header,sheet);
        
        % get indices for time range and adjust time
        jjT = SimResult.time >= timelimit(iT,1) & SimResult.time <= timelimit(iT,2);
        time = (SimResult.time(jjT) - timelimit(iT,1)).*timeUnitFactor;

        % get indices for time range and adjust time for reference simulation
        if ~isempty(SimResultRef)
            jjTRef = SimResultRef.time >= timelimit(iT,1) & SimResultRef.time <= timelimit(iT,2);
            timeRef = (SimResultRef.time(jjT) - timelimit(iT,1)).*timeUnitFactor;
        else
            timeRef = [];
        end

        % switch tad or Time
        if timelimit(iT,1)>0;
            timeLabel = 'Time after dose';
        else
            timeLabel = 'Time';
        end
            
        
        % loop on Outputs
        for iO = 1:length(SimResult.outputPathList)
        
            % get Y values for time range
            y = SimResult.values{iO}(jjT,:).*OutputList(iO).unitFactor;
            
            % check if  reference output exists
            yRef = [];
            if ~isempty(SimResultRef)
                jj = strcmp(OutputList(iO).pathID,SimResultRef.outputPathList);
                if any(jj)
                    yRef = SimResultRef.values{jj}(jjTRef,:).*OutputList(iO).unitFactor;
                end
            end
            
            % get data for this timlimit           
            [TP_tmp,lloq] = prepareDataForTimeRange(TP,iO,timelimit(1,:));
            
            % without data LLOQ = nan
    
            % loop on scale
            for iScale = 1:length(yscale)
                
                
                % get name and figure description
                figureName = sprintf('O%d_%s_%s',iO,removeForbiddenLetters(OutputList(iO).reportName),yscale{iScale});
                [figtxt,figtxtTable,legendEntries] = feval(textFunctionHandle,WSettings,'tpShadedArea',...
                    {OutputList(iO).reportName,PopRunSet(iSet).reportName,PopRunSet(iSet).popReportName,...
                        reportNameRef, popReportNameRef,dataReportName,yscale{iScale},...
                    lloq,length(TP_tmp)});
                
                
                % do the plot
                csv =  plotReportTimeProfile(WSettings,FP.figureHandle,time,y,timeRef,yRef,TP_tmp,timeLabel,Def.timeDisplayUnit,...
                    OutputList(iO).reportName,OutputList(iO).displayUnit,legendEntries,yscale{iScale},lloq);
            
                
                % save figure
                if iScale==1
                    FP = FP.printFigure(figureName,figtxt,csv,figtxtTable);
                else
                    FP = FP.printFigure(figureName,figtxt);
                end

                
                
            end % iScale    
            
        end % iO

    end % iT
    
end % iSet 


return
    

function SimResult = loadSim(simulationName)

fname = fullfile('tmp',simulationName,'simResult.mat');

if ~exist(fname,'file')
    SimResult = readPopulationResultfile(simulationName);
else
    load(fname);
end



function  [TP_tmp,lloq] = prepareDataForTimeRange(TP,iO,timelimit)

    % initialize return values
TP_tmp = [];
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
            
% get stucture to plo
for iInd = 1:length(TP{iO})
    
    jjT = TP{iO}(iInd).time >= timelimit(1,1) & TP{iO}(iInd).time <= timelimit(1,2);

    TP_tmp(iInd).time = TP{iO}(iInd).(timefield)-timeoffset; %#ok<AGROW>
    TP_tmp(iInd).y(:,1) = TP{iO}(iInd).dv(jjT); %#ok<AGROW>
    
    % set lloq
    isLloq = TP{iO}(iInd).isLloq(jjT);
    TP_tmp(iInd).lloq(:,1) = nan(length(isLloq),1); %#ok<AGROW>
    if any(isLloq)
        TP_tmp(iInd).y(isLloq,1) = nan; %#ok<AGROW>
        TP_tmp(iInd).lloq(~isLloq,1) = TP(iInd).lloq/2; %#ok<AGROW>
        lloq = TP(iInd).lloq;
    end
end

return
