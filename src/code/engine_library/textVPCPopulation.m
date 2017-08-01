function [figtxt,figtxtTable,legendEntries] = textVPCPopulation(WSettings,figureType,strArray,numVal)
% TEXTVPCPOPULATION creates text fo figures and tables
%
% Inputs: 
%       WSettings (structure)    definition of properties used in all
%                   workflow functions see GETDEFAULTWORKFLOWSETTINGS
%

% Open Systems Pharmacology Suite;  http://forum.open-systems-pharmacology.org
% Date: 1-Aug-2017


% initialize outputs
figtxt = '';
figtxtTable = '';
legendEntries = {};

switch figureType
    
    % histogram for physiological properties
    case 'physHist'
        [figtxt,figtxtTable,legendEntries] = textVPCPopulationPhysHist(strArray{1},strArray{2},numVal(1),strArray{3},numVal(2));

    % Time profile
    case 'tpShadedArea'
        [figtxt,figtxtTable,legendEntries] = textVPCPopulationTpShadedArea(strArray{1},strArray{2},strArray{3},...
            strArray{4},strArray{5},strArray{6});
        
    % Boxwhisker of PK Parmeter
    case 'pkBW'
        [figtxt,figtxtTable,legendEntries] = textVPCPopulationPkBW(WSettings,strArray{1},strArray{2},strArray{3},strArray{4},strArray{5});
end


return

function  [figtxt,figtxtTable,legendEntries] = textVPCPopulationPhysHist(physProb,population,nInd,referencePopulation,nIndRef)

figtxt = sprintf('Percentile of %s for %s.',physProb,population);
if ~isempty(referencePopulation)
    figtxt = sprintf('%s in comparison to %s.',figtxt(1:end-1),referencePopulation);
end

figtxtTable = figtxt;

legendEntries = {sprintf('%s: n=%d',population,nInd),sprintf('%: n=%d',referencePopulation,nIndRef)};

return


function [figtxt,figtxtTable,legendEntries] = textVPCPopulationTpShadedArea(output,simulation,pop,simulationRef,popRef,scale)

            
% get name and figure description
figtxtTable = sprintf('Time profile of %s for %s for %s.',...
    output,simulation,pop);
                
if isempty(simulationRef)
    figtxt = figtxtTable;
else
    figtxt = sprintf('Reference: %s for %s.',figtxtTable,simulationRef,...
        popRef);
end
                
switch scale
    case 'lin'
        figtxt = sprintf('%s Time profiles are plotted on a linear scale.',figtxt);
    case 'log'
        figtxt = sprintf('%s Time profiles are plotted on a logarithmic scale.',figtxt);
end

if isempty(simulationRef)
    legendEntries = {''};
else
    legendEntries = {pop,popRef};
end

return

function  [figtxt,figtxtTable,legendEntries] = textVPCPopulationPkBW(WSettings,yLabel,output,reportNames,popReportNames,scale)


% get table header
figtxtTable = sprintf('Percentile of %s of %s.',yLabel,output);

% check if extrameas are plotted
extremaTxt = '';
if WSettings.withExtrema
 extremaTxt = ' and extremes as open circles';
end

% get figure description
figtxt = sprintf('%s of %s shown as box whisker plot, which indicate the percentiles %d, %d, %d, %d, and %d%s',...
    yLabel,output,WSettings.displayPercentiles(1),WSettings.displayPercentiles(2),WSettings.displayPercentiles(3),...
    WSettings.displayPercentiles(4),WSettings.displayPercentiles(5),extremaTxt);

% set scale text
switch scale
    case 'lin'
        figtxt = sprintf('%s on a linear scale.',figtxt);
    case 'log'
        figtxt = sprintf('%s on a logarithmic scale.',figtxt);
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