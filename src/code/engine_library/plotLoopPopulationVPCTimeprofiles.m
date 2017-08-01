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
%       PopRunSet (structure)   list of population simulations see GETDEFAULTPOPRUNSET
%       FP (ReportFigurPrint)  objects with manages print of figures and
%                           tables
%   Outputs:
%       FP (ReportFigurPrint)  objects with manages print of figures and

% Open Systems Pharmacology Suite;  http://forum.open-systems-pharmacology.org
% Date: 31-July-2017

% get Reference simulation
if ~isempty(Def.ixPopRunSetRef)
    SimResultRef = loadSim({PopRunSet(Def.ixPopRunSetRef).name});
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
        
    if ~isValid || length(startTimes) ==1
        timelimit = SimResult.time([ 1 end])';
        timeRangetxt = {};
    else
        timelimit = [SimResult.time([ 1 end])';...
            startTimes(1:2);...
            startTimes(end) SimResult.time(end)];
        timeRangetxt = {'total simulation range','first application range','last application range'};
    end
    
    % load OutputList of population
    load(fullfile('tmp',PopRunSet(iSet).name,'outputList.mat'));
    
    for iT = 1:size(timelimit,1)
        
        %initialize new sheet
        header = sprintf('Time profiles of %s',PopRunSet(iSet).reportName);
        sheet = PopRunSet(iSet).name;
        if ~isempty(timeRangetxt)
            header = sprintf('%s of the %s',header,timeRangetxt{iT});
            sheet = sprintf('%s_%s',sheet,removeForbiddenLetters(timeRangetxt{iT}));
        end
        FP = FP.iniCaptiontextFigtextArray(header,sheet);
        
        % get index for time range and adjust time
        jj_t = SimResult.time >= timelimit(iT,1) & SimResult.time <= timelimit(iT,2);
        time = (SimResult.time(jj_t) - timelimit(iT,1)).*timeUnitFactor;

        % get index for time range and adjust time for reference simulation
        if ~isempty(SimResultRef)
            jj_t_ref = SimResultRef.time >= timelimit(iT,1) & SimResultRef.time <= timelimit(iT,2);
            timeRef = (SimResultRef.time(jj_t) - timelimit(iT,1)).*timeUnitFactor;
        else
            timeRef = [];
        end

        % loop on Outputs
        for iO = 1:length(SimResult.outputPathList)
        
            % get Y values for time range
            y = SimResult.values{iO}(jj_t,:).*OutputList(iO).unitFactor;
            
            % check if  reference output exists
            yRef = [];
            if ~isempty(SimResultRef)
                jj = strcmp(OutputList(iO).pathID,SimResultRef.outputPathList);
                if any(jj)
                    yRef = SimResult.values{jj}(jj_t_ref,:).*OutputList(iO).unitFactor;
                end
            end
            
            % without data LLOQ = nan
            lloq = nan;
    
            % loop on scale
            for iScale = 1:length(yscale)
                
                
                % get name and figure description
                figureName = sprintf('O%d_%s_%s',iO,removeForbiddenLetters(OutputList(iO).reportName),yscale{iScale});
                [figtxt,figtxtTable,legendEntries] = feval(textFunctionHandle,WSettings,'tpShadedArea',...
                    {OutputList(iO).reportName,PopRunSet(iSet).reportName,PopRunSet(iSet).popReportName,...
                        reportNameRef, popReportNameRef, yscale{iScale}},...
                    []);
                
                
                % do the plot
                csv =  plotReportTimeProfile(WSettings,FP.figureHandle,time,y,timeRef,yRef,Def.timeDisplayUnit,...
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