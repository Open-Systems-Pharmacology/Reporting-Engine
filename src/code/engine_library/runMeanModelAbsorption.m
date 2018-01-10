function runMeanModelAbsorption(WSettings,MeanModelSet)
% RUNMEANMODELABSORPTION gets absorption characteristic of the applications at time 0
%
% runMeanModelAbsorption(WSettings,MeanModelSet)
%
% Inputs
%       WSettings (structure)    definition of properties used in all
%                   workflow functions see GETDEFAULTWORKFLOWSETTINGS
%       MeanModelSet (structure)   list of population simulations see GENERATEWORKFLOWINPUTFORMEANMODELSIMULATION

% Open Systems Pharmacology Suite;  http://open-systems-pharmacology.org

try
    writeToReportLog('INFO',sprintf('Start generate absorption plots'),false);
    
    % Initialize figureDir
    FP = ReportFigurePrint(fullfile(WSettings.figures,'absorption'),WSettings.printFormatList);
    FP = FP.iniCaptiontextFigtextArray('Absorption characteristics','absorption');
    
    for iSet = 1:length(MeanModelSet)
        
        [time,R] = simulateAbsorption(WSettings,MeanModelSet(iSet));
        
        if isempty(R)
            writeToReportLog('WARNING',sprintf('No compound is absorbed in %s',MeanModelSet(iSet).name),false);
        else
            FP = plotAbsorption(WSettings,MeanModelSet(iSet).name,time,R,FP);
        end
    end
    
    FP.saveCaptiontextArray;
    
    
    writeToReportLog('INFO',sprintf('Finalized Absorption plots \n'),false);
    
 catch exception
        
    save(sprintf('exception_%s.mat',datestr(now,'ddmmyy_hhMM')),'exception');
    writeToReportLog('ERROR',exception.message,false);
    writeToReportLog('INFO',sprintf('Absorption plots finished with error \n'),false);
        
end   
    
return


function FP = plotAbsorption(WSettings,simulationName,time,R,FP)

% get factor to translate time in  diplay units
displayUnit = 'h';
timeUnitFactor = getUnitFactor('',displayUnit,'Time');
time = time.*timeUnitFactor;

for iComp = 1:length(R)

    % create figure
    ax = getReportFigure(WSettings,1,1,FP.figureHandle,'figureformat','landscape');
    
    
    % legend handle array for plots
    lgh=[];
    
    
    lgh(1) = plot(time(2:end),R(iComp).fDiss(2:end),'k-','linewidth',2,'displayname','dissolved');
    lgh(2) = plot(time(2:end),R(iComp).fAbs(2:end),'r--','linewidth',2,'displayname','absorbed to mucosa');
    lgh(3) = plot(time(2:end),R(iComp).fGI(2:end),'g:','linewidth',2,'displayname','absorbed to portal vein');
    lgh(4) = plot(time(2:end),R(iComp).fBio(2:end),'b-.','linewidth',2,'displayname','absorbed to venous blood');
    lgh(5) = plot(time(2:end),R(iComp).fExcreted(2:end),'c-','linewidth',2,'displayname','excreted to feces');
    
    setAxesScaling(ax,'timeUnit',displayUnit,'xlim',[time(1) time(end)]);
    xlabel(sprintf('Time [%s]',displayUnit));
    ylabel('Fraction of drugmass');
    
    % set legend
    jj = lgh>0;
    legend(lgh(jj),get(lgh(jj),'displayname'),'location','northoutside');

    figureName = sprintf('%s_%s',simulationName,removeForbiddenLetters(R(iComp).compound));
    figtxt = sprintf('Absorption of %s',R(iComp).compound);
    
    csv = [{sprintf('time [%s]',displayUnit),'fraction dissolved','fraction absorbed','fraction absorbed to portal vein',...
        'fraction absorbed to venous blood','fraction excreted to feces'};
    num2cell([time;R(iComp).fDiss;R(iComp).fAbs;R(iComp).fGI;R(iComp).fBio;R(iComp).fExcreted])'];
    
    
    FP = FP.printFigure(figureName,figtxt,csv,figtxt);

end

return

function [time,R] = simulateAbsorption(WSettings,MeanModelSet) %#ok<INUSL>


initStruct=[];
initStruct=initParameter(initStruct,'*|Organism|Lung|Blood flow rate','always','throwWarningIfNotExisting',false);
initStruct=initParameter(initStruct,'*|Organism|PortalVein|Blood flow rate','always','throwWarningIfNotExisting',false);
initSimulation(MeanModelSet.xml,initStruct,'report','none');

% Only one simulation is intialized so the simulation index is always 1
simulationIndex = 1;


% load applicationProtocoll
load(fullfile('tmp',MeanModelSet.name,'applicationProtocol.mat'),'ApplicationProtocol','isValid');

% find interval of first application
if isValid
    
    % get startin application
    jjZero = [ApplicationProtocol.startTime] ==0;         %#ok<NODEF>
    compounds = unique({ApplicationProtocol(jjZero).compound});
    
    % get list of applicated compounds
    for iComp = 1:length(compounds)
        jjComp = jjZero & strcmp({ApplicationProtocol(:).compound},compounds{iComp}); %#ok<IDISVAR>
        R(iComp).compound = compounds{iComp}; %#ok<AGROW>
        R(iComp).drugmass = sum(ApplicationProtocol(jjComp).drugMass); %#ok<AGROW>
        R(iComp).fDiss = []; %#ok<AGROW>
        R(iComp).fAbs = []; %#ok<AGROW>
        R(iComp).fGI = []; %#ok<AGROW>
        R(iComp).fBio = []; %#ok<AGROW>
        R(iComp).fExcreted = []; %#ok<AGROW>
    end
    
    % shortn simulation time
    if any(~jjZero)
        startTimes = sort([ApplicationProtocol.startTime]);
        ix = find(startTimes>0,1);
        time = getSimulationTime;
        jjT = time <= startTimes(ix);
        setSimulationTime(time(jjT));
    end
    
else
    error('drugmass can not be analysed')
end

setRelativeParameter(1,'*|Organism|Lung|Blood flow rate',simulationIndex);
setRelativeParameter(1,'*|Organism|PortalVein|Blood flow rate',simulationIndex);

processSimulation(simulationIndex);
time = getSimulationTime(simulationIndex);

% fraction dissolved
for iComp = 1:length(R)
    
    % fraction dissolved
    pts = sprintf('*|Organism|Lumen|%s|Fraction dissolved',R(iComp).compound);
    if existsObserver(pts,simulationIndex)
        [time,tmpValues] = getSimulationResult(pts,simulationIndex);
        R(iComp).fDiss = tmpValues; 
    else
        R(iComp).fDiss = nan(size(time));
    end
    
    
    % fraction absorbed
    pts = sprintf('*|Organism|Lumen|%s|Fraction of oral drug mass absorbed into mucosa',R(iComp).compound);
    if existsObserver(pts,simulationIndex)
        [time,tmpValues] = getSimulationResult(pts,simulationIndex);
        R(iComp).fAbs = tmpValues; 
    else
        R(iComp).fAbs = nan(size(time));
    end
    
    % fraction excretd
    pts = sprintf('*|Organism|Lumen|Feces|%s|Fraction excreted to feces',R(iComp).compound);
    if existsObserver(pts,simulationIndex)
        [time,tmpValues] = getSimulationResult(pts,simulationIndex);
        R(iComp).fExcreted = tmpValues; 
    else
        R(iComp).fExcreted = nan(size(time));
    end
    
end

% get rid of non absorbed
jj = ~cellfun(@(x) all(x==0),{R.fAbs});
R = R(jj);

if isempty(R)
    return
end
    
% check first pass
setParameter(0,'*|Organism|Lung|Blood flow rate',simulationIndex);
processSimulation(simulationIndex);
for iComp = 1:length(R)
    
    [time,tmpValues] = getSimulationResult(['*|Organism|VenousBlood|*|' R(iComp).compound],1);
    
    R(iComp).fBio = sum(tmpValues)./ R(iComp).drugmass; 
end


% check metabolisation rate
setParameter(0,'*|Organism|PortalVein|Blood flow rate',simulationIndex);
processSimulation(simulationIndex);
for iComp = 1:length(R)
    
    [time,tmpValues] = getSimulationResult(['*|Organism|PortalVein|*|' R(iComp).compound],1);
    
    R(iComp).fGI = sum(tmpValues)./R(iComp).drugmass; 
end


return
