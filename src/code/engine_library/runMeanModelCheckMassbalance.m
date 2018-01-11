function runMeanModelCheckMassbalance(WSettings,MeanModelSet,MBS)
% RUNMEANMODELCHECKMASSBALANCE gets absorption characteristic of the applications at time 0
%
% runMeanModelCheckMassbalance(WSettings,MeanModelSet,MassbalanceSettings)
%
% Inputs
%       WSettings (structure)    definition of properties used in all
%                   workflow functions see GETDEFAULTWORKFLOWSETTINGS
%       MeanModelSet (structure)   list of population simulations see GENERATEWORKFLOWINPUTFORMEANMODELSIMULATION
%       MassbalanceSettings (structure) settings of massbalance plots

% Open Systems Pharmacology Suite;  http://open-systems-pharmacology.org

writeToReportLog('INFO',sprintf('Start generate Massbalance plots'),false);

% Initialize figureDir
FP = ReportFigurePrint(fullfile(WSettings.figures,'massbalance'),WSettings.printFormatList);

for iSet = 1:length(MeanModelSet)

    [time,drugmass,R] =simulateModel(WSettings,MeanModelSet(iSet),MBS);

    % start plots
    
    col = getcolmarkForMap(jet,length(R));
    FP = FP.iniCaptiontextFigtextArray(sprintf('Massbalance for %s',MeanModelSet(iSet).reportName),MeanModelSet(iSet).name);
    
    % timprofile  amount lin
    csv = plotTP(WSettings,FP.figureHandle,time,1,R,'lin',MBS,'amount [µmol]',col);
    
    figureName = sprintf('TPamount_%s_lin',MeanModelSet(iSet).name);
    figtxt = sprintf('Amount of drug vs time within the different compartments. Time profiles are in a linear scale.');
        
    FP = FP.printFigure(figureName,figtxt,csv,figtxt);

    % timprofile  amount log
    plotTP(WSettings,FP.figureHandle,time,1,R,'log',MBS,'amount [µmol]',col);
    
    figureName = sprintf('TPamount_%s_log',MeanModelSet(iSet).name);
    figtxt = sprintf('Amount of drug vs time within the different compartments. Time profiles are in a logarithmic scale.');
        
    FP = FP.printFigure(figureName,figtxt);
    
    % timprofile  fraction lin
    plotTP(WSettings,FP.figureHandle,time,drugmass,R,'lin',MBS,'amount [fraction of drugmass]',col);
    
    figureName = sprintf('TPfraction_%s_lin',MeanModelSet(iSet).name);
    figtxt = sprintf('Amount of drug vs time within the different compartments normalized to applicated drugmass. Time profiles are in a linear scale.');
        
    FP = FP.printFigure(figureName,figtxt);

    % timprofile  fraction log
    plotTP(WSettings,FP.figureHandle,time,drugmass,R,'log',MBS,'amount [fraction of drugmass]',col);
    
    figureName = sprintf('TPfraction_%s_log',MeanModelSet(iSet).name);
    figtxt = sprintf('Amount of drug vs time within the different compartments normalized to applicated drugmass. Time profiles are in a logarithmic scale.');
        
    FP = FP.printFigure(figureName,figtxt);
    
    % cum timprofile  amount
    plotTPCum(WSettings,FP.figureHandle,time,1,R,'lin',MBS,'cumulated amount [µmol]',col);
    
    figureName = sprintf('cumTPamount_%s_lin',MeanModelSet(iSet).name);
    figtxt = sprintf('Cumulated amount of drug vs time within the different compartments. Time profiles ar in a linear scale.');

    FP = FP.printFigure(figureName,figtxt);
        
    % cum timprofile  fraction
    plotTPCum(WSettings,FP.figureHandle,time,drugmass,R,'lin',MBS,'cumulated amount [fraction of drugmass]',col);

    figureName = sprintf('cumTPfraction_%s_lin',MeanModelSet(iSet).name);
    figtxt = sprintf('Cumulated amount of drug vs time within the different compartments normalized to applicated drugmass. Time profiles ar in a linear scale.');
    
    FP = FP.printFigure(figureName,figtxt);

    % pie plots
    for iT = 1:length(MBS.piePlotTime)
        
        if MBS.piePlotTime==-1
            tPoint = time(end); 
        else
            tPoint = MBS.piePlotTime(iT);
        end
        
        plotPie(WSettings,FP.figureHandle,time,tPoint,R,col);
        
        figureName = sprintf('pie_%s_t%g',MeanModelSet(iSet).name,tPoint);
        figtxt = sprintf('Fraction of drug  within the different compartments a %g%s.',tPoint,MBS.displayUnit);
    
        FP = FP.printFigure(figureName,figtxt);

    end
    
    FP.saveCaptiontextArray;
    
    
end


writeToReportLog('INFO',sprintf('Finalized Massbalance plots \n'),false);

return

function plotPie(WSettings,fh,time,tPoint,R,col)

getReportFigure(WSettings,1,1,fh,'figureformat','portrait');

x = (reshape([R.v],length(time),length(R)));
xPoint = interp1(time,x,tPoint);
xPoint = xPoint./sum(xPoint);

theta0 = pi/2;
maxpts = 100;


for iR=1:length(xPoint)
  n = max(1,ceil(maxpts*xPoint(iR)));
  r = [0;ones(n+1,1);0];
  theta = theta0 + [0;xPoint(iR)*(0:n)'/n;0]*2*pi;
  [xx,yy] = pol2cart(theta,r);
  theta0 = max(theta);
  lgh(iR) = patch(xx,yy,col(iR,:),'displayname',sprintf('%s %s (%.2g%%)',R(iR).compound,R(iR).name,xPoint(iR)*100)); %#ok<AGROW>
  hold on;
end

legend(lgh,get(lgh,'Displayname'),'location','northoutside');

return



function  plotTPCum(WSettings,fh,time,drugmass,R,yscale,MBS,ylabelTxt,col)

ax = getReportFigure(WSettings,1,1,fh,'figureformat','portrait');

vUpper = R(1).v./drugmass;
if strcmp(yscale,'log')
    jj =vUpper>0;
    vLower=ones(size(time)).*min(vUpper(jj))/2;
    vUpper(~jj) = vLower(1);
else
    vLower=zeros(size(time));
end
yl(1) = vLower(1);

for iR =1:length(R)
    
    lgh(iR) = patch([time time(end:-1:1)],[vUpper vLower(end:-1:1)],col(iR,:),'displayname',sprintf('%s %s',R(iR).compound,R(iR).name)); %#ok<AGROW>
    
    if iR < length(R)
        vLower = vUpper;
        vUpper = vUpper + R(iR+1).v./drugmass;
    end
end
yl(2) = max(vUpper);

ylabel(ax,ylabelTxt);
xlabel(ax, sprintf('time [%s]',MBS.displayUnit));

setAxesScaling(ax,'timeUnit',MBS.displayUnit,'xlim',[0 time(end)]);
set(ax,'yscale',yscale,'ylim',yl);

legend(lgh,get(lgh,'displayname'),'location','northoutside');

return


function csv = plotTP(WSettings,fh,time,drugmass,R,yscale,MBS,ylabelTxt,col)

csv = {};

ax = getReportFigure(WSettings,1,1,fh,'figureformat','portrait');

for iR=1:length(R)
    
    lgh(iR)=plot(ax,time,R(iR).v./drugmass,'color',col(iR,:),'linewidth',2,'displayname',sprintf('%s %s',R(iR).compound,R(iR).name)); %#ok<AGROW>
end

ylabel(ax,ylabelTxt);
xlabel(ax, sprintf('time [%s]',MBS.displayUnit));

set(ax,'yscale',yscale);
setAxesScaling(ax,'timeUnit',MBS.displayUnit,'xlim',[0 time(end)]);

legend(lgh,get(lgh,'displayname'),'location','northoutside');


% construct csv
csv{1,1} = sprintf('time [%s]',MBS.displayUnit);
csv(2:(1+length(time)),1) = cellfun(@num2str,num2cell(time),'uniformoutput',false);

for iR=1:length(R)
    csv{1,1+iR} = sprintf('%s %s [µmol]',R(iR).compound,R(iR).name); 
    csv(2:(1+length(time)),1+iR) = cellfun(@num2str,num2cell(R(iR).v),'uniformoutput',false);  
end

return



function [time,drugmass,R] = simulateModel(WSettings,MeanModelSet,MBS) %#ok<INUSL>


initStruct=[];
initSimulation(MeanModelSet.xml,initStruct,'report','none');

% Only one simulation is intialized so the simulation index is always 1
simulationIndex = 1;


% load applicationProtocoll
load(fullfile('tmp',MeanModelSet.name,'applicationProtocol.mat'),'ApplicationProtocol','isValid');

% find interval of first application
if isValid
    
    % get application for relevant compounds
    if isempty(MBS.excludedCompounds)
        jjComp = true(1,length(ApplicationProtocol));
    else
        jjComp = ~ismember({ApplicationProtocol.compound},MBS.excludedCompounds);
    end
    
    drugmass = sum([ApplicationProtocol(jjComp).drugMass]); 
        
else
    error('drugmass can not be analysed')
end

% do the simulation and return simulation time
processSimulation(simulationIndex);

time = getSimulationTime.*getUnitFactor('',MBS.displayUnit,'time');

% get List of all Compounds 
[~,desc]=existsSpeciesInitialValue('*',1,'parametertype','readonly');
tmp = regexp(desc(2:end,strcmp(desc(1,:),'Path')),'\|','split');
compoundList = unique(cellfun(@(x) x{end},tmp,'uniformoutput',false));
compoundList = setdiff(compoundList,{'Liquid','Oral mass absorbed segment','Suspension fraction Cecum',...
    'Suspension fraction Colon Ascendens','Suspension fraction Colon Descendens',...
    'Suspension fraction Colon Sigmoid','Suspension fraction Colon Transversum',...
    'Suspension fraction Duodenum','Suspension fraction Lower Ileum',...
    'Suspension fraction Lower Jejunum','Suspension fraction Rectum',...
    'Suspension fraction Stomach','Suspension fraction Upper Ileum',...
    'Suspension fraction Upper Jejunum','Tablet_Location'});


% delete excluded species
compoundList = setdiff(compoundList,MBS.excludedCompounds);

R = struct('name','','compound','','v',[]);
R=R([]);

% collect compounds
for iC = 1:length(compoundList)
    
    % get vector of all IDs
    [isExisting,desc] = existsSpeciesInitialValue(['*|Organism|*' compoundList{iC}],1,'parametertype','readonly');
    if isExisting
        Sdesc = cell2struct(desc(2:end,:),desc(1,:),2);
        tmp = regexp({Sdesc.Path},'\|','split');
        
        % separate Lumen
        jjLumen = cellfun(@(x) strcmp(x{end-2},'Lumen'),tmp);
        jjFeces = cellfun(@(x) strcmp(x{end-1},'Feces'),tmp);
        
        % get compartments
        compartments = setdiff(unique(cellfun(@(x) x{end-1},tmp(~jjLumen),'uniformoutput',false)),{'Saliva','SalivaGland'});
        
        if any(jjLumen)
            R = addSplit([Sdesc(jjLumen & ~jjFeces).ID],'Lumen',compoundList{iC},R);
        end
        
        for iComp = 1:length(compartments)
            
            jj = cellfun(@(x) strcmp(x{end-1},compartments{iComp}),tmp);
            
            R = addSplit([Sdesc(jj).ID],compartments{iComp},compoundList{iC},R);
            
        end
        if any(jjFeces)
            R = addSplit([Sdesc(jjFeces).ID],'Feces',compoundList{iC},R);
        end
    end
end

return

function [R] = addSplit(IDList,name,compound,R)

[~,v] = getSimulationResult(IDList,1);
v = sum(v,1);

% construct structure
if any(v) > 0
    
    R(end+1).name = name;
    R(end).compound = compound;
    R(end).v = v;
end

return
