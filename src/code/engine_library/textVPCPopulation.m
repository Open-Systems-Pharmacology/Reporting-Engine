function [figtxt,figtxtTable,legendEntries] = textVPCPopulation(WSettings,figureType,inputArray)
% TEXTVPCPOPULATION creates text fo figures and tables
%
% Inputs: 
%       WSettings (structure)    definition of properties used in all
%                   workflow functions see GETDEFAULTWORKFLOWSETTINGS
%

% Open Systems Pharmacology Suite;  http://open-systems-pharmacology.org



% initialize outputs
figtxt = '';
figtxtTable = '';
legendEntries = {};

switch figureType
    
    % histogram for physiological properties
    case 'physHist'
        [figtxt,figtxtTable,legendEntries] = textVPCPopulationPhysHist(inputArray{1},inputArray{2},inputArray{3},...
            inputArray{4},inputArray{5},inputArray{6});

    % shaded Area for physiological properties
    case 'physShA'
        [figtxt,figtxtTable,legendEntries] = textVPCPopulationphysShA(inputArray{1},inputArray{2},inputArray{3},inputArray{4},...
            inputArray{5},inputArray{6});

    % Time profile
    case 'tpShadedArea'
        [figtxt,figtxtTable,legendEntries] = textVPCPopulationTpShadedArea(inputArray{1},inputArray{2},inputArray{3},...
            inputArray{4},inputArray{5},inputArray{6},inputArray{7},inputArray{8},inputArray{9});
        
    % Boxwhisker of PK Parmeter
    case 'pkBW'
        [figtxt,figtxtTable,legendEntries] = textVPCPopulationPkBW(WSettings,inputArray{1},inputArray{2},inputArray{3},inputArray{4},inputArray{5});
        
    % shaded Area for PK Parameter
    case 'pkShA'
        [figtxt,figtxtTable,legendEntries] = textVPCPopulationPkShA(inputArray{1},inputArray{2},inputArray{3},inputArray{4},...
            inputArray{5},inputArray{6});

end


return

function  [figtxt,figtxtTable,legendEntries] = textVPCPopulationPhysHist(physProb,population,dataSource,popLabels,nPop,nData)


figtxt = sprintf('Percentile of %s for %s.',physProb,population);
if nData>0
    figtxt = sprintf('%s Data %s.',figtxt(1:end-1),dataSource);
end

figtxtTable = figtxt;

for iPop = 1:length(popLabels)
    legendEntries{iPop} = sprintf('%s: n=%d',popLabels{iPop},nPop(iPop)); %#ok<AGROW>
end
if nData>0
    legendEntries{end+1} = sprintf('%s: n=%d',dataSource,nData);
end

return


function  [figtxt,figtxtTable,legendEntries] = textVPCPopulationphysShA(physProbX,physProbY,population,referencePopulation,dataSource,nInd)


figtxt = sprintf('%s-dependence of %s for %s.',physProbX,physProbY,population);
% Reference
if nInd(2)>0
    figtxt = sprintf('%s in comparison to %s.',figtxt(1:end-1),referencePopulation);
end
% Data
if nInd(3)>0
    figtxt = sprintf('%s Data %s.',figtxt(1:end-1),dataSource);
end

figtxtTable = figtxt;

legendEntries = {sprintf('%s: n=%d',population,nInd(1)),sprintf('%s: n=%d',referencePopulation,nInd(2)),...
    sprintf('%s: n=%d',dataSource,nInd(3))};

return

function [figtxt,figtxtTable,legendEntries] = textVPCPopulationTpShadedArea(output,simulation,pop,...
    simulationRef,popRef,popData,scale,lloq,nData)

            
% get name and figure description
figtxtTable = sprintf('Time profile of %s for %s for %s.',...
    output,simulation,pop);
                
if isempty(simulationRef)
    figtxt = figtxtTable;
else
    figtxt = sprintf('%s Reference: %s for %s.',figtxtTable,simulationRef,...
        popRef);
end
if nData>0
    figtxt = sprintf('%s Data source: %s.',figtxt,popData);
end    
if ~isnan(lloq)
    figtxt = sprintf('%s Data below lower limit of quantification (lloq) are plotted as open symbols as lloq/2.',figtxt,popData);
end

switch scale
    case 'lin'
        figtxt = sprintf('%s Time profiles are plotted in a linear scale.',figtxt);
    case 'log'
        figtxt = sprintf('%s Time profiles are plotted in a logarithmic scale.',figtxt);
end


% legend
if isempty(simulationRef)
    legendEntries = {''};
else
    legendEntries = {pop,popRef};
end
if nData>0
    legendEntries{3} = popData;
end
if ~isnan(lloq)
     legendEntries{end+1} = 'lower limit of quantification';
end

return

function  [figtxt,figtxtTable,legendEntries] = textVPCPopulationPkBW(WSettings,yLabel,output,reportNames,popReportNames,scale)


% get table header
figtxtTable = sprintf('Percentile of %s of %s.',yLabel,output);

% check if extrameas are plotted
extremaTxt = '';
if WSettings.boxwhiskerWithExtrema
 extremaTxt = ' and extremes as open circles';
end

% get figure description
figtxt = sprintf('%s of %s shown as box whisker plot, which indicate the percentiles %d, %d, %d, %d, and %d%s',...
    yLabel,output,WSettings.displayPercentiles(1),WSettings.displayPercentiles(2),WSettings.displayPercentiles(3),...
    WSettings.displayPercentiles(4),WSettings.displayPercentiles(5),extremaTxt);

% set scale text
switch scale
    case 'lin'
        figtxt = sprintf('%s in a linear scale.',figtxt);
    case 'log'
        figtxt = sprintf('%s in a logarithmic scale.',figtxt);
end

% get xLabels
if length(reportNames) ==1
    legendEntries = {''};
elseif length(unique(popReportNames)) == length(popReportNames);
    legendEntries = popReportNames;
elseif  length(unique(reportNames)) == length(reportNames);
     legendEntries = reportNames;
else
     legendEntries = strcat(popReportNames,'; ',reportNames);
end

return


function  [figtxt,figtxtTable,legendEntries] = textVPCPopulationPkShA(xPhys,yLabel,output,reportNames,refReportName,scale)


% get table header
figtxtTable = sprintf('%s dependency of %s for %s for a %s',xPhys,yLabel,output,reportNames);


% set scale text
switch scale
    case 'lin'
        figtxt = sprintf('%s in a linear scale.',figtxtTable);
    case 'log'
        figtxt = sprintf('%s in a logarithmic scale.',figtxtTable);
end

legendEntries={reportNames,refReportName};

return