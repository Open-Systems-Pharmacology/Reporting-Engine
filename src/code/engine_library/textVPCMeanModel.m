function [figtxt,figtxtTable,legendEntries] = textVPCMeanModel(WSettings,figureType,inputArray) %#ok<INUSL>
% TEXTVPCMEANMODEL creates text fo figures and tables
%
% Inputs: 
%       WSettings (structure)    definition of properties used in all
%                   workflow functions see GETDEFAULTWORKFLOWSETTINGS
%       figureType (string) type of figure for which the text is needed
%       inputArray (cellarray) inputs, differ for each histogram type
%

% Open Systems Pharmacology Suite;  http://open-systems-pharmacology.org



% initialize outputs
figtxt = '';
figtxtTable = '';
legendEntries = {};

switch figureType
    
    % Time profile
    case 'tpShadedArea'
        [figtxt,figtxtTable,legendEntries] = textVPCTimeprofileTpShadedArea(inputArray{1},inputArray{2},inputArray{3},...
            inputArray{4},inputArray{5},inputArray{6});
        
    % Predicted vs observed
    case 'tpPredVsObs'
        [figtxt,legendEntries] = textVPCPredictedVsObserved(inputArray{1},inputArray{2},inputArray{3},inputArray{4},inputArray{5});
         
     % residuals vs Time 
    case 'tpResVsTime'
        [figtxt,legendEntries] = textVPCResidualsVsTime(inputArray{1},inputArray{2},inputArray{3},inputArray{4},inputArray{5});
        
     % residuals vs Prediction
    case 'tpResVsY'
        [figtxt,legendEntries] = textVPCResidualsVsY(inputArray{1},inputArray{2},inputArray{3},inputArray{4},inputArray{5});
        
    % residuals as Histogram
    case 'histRes'
        [figtxt,legendEntries] = textVPCResidualsAsHistogram(inputArray{1},inputArray{2});


    % residuals as qqPlot
    case 'qqRes'
        [figtxt,legendEntries] = textVPCResidualsAsQQPlot(inputArray{1},inputArray{2});

end


return


function [figtxt,figtxtTable,legendEntries] = textVPCTimeprofileTpShadedArea(output,simulation,nameData,scale,lloq,nData)

            
% get name and figure description
figtxtTable = sprintf('Time profile of %s for %s.',...
    output,simulation);
                

figtxt = figtxtTable;
if nData>0
    figtxt = sprintf('%s Data source: %s.',figtxt,nameData);
end    
if ~isnan(lloq)
    figtxt = sprintf('%s Data below lower limit of quantification (lloq) are plotted as open symbols as lloq/2.',figtxt);
end

switch scale
    case 'lin'
        figtxt = sprintf('%s Time profiles are plotted in a linear scale.',figtxt);
    case 'log'
        figtxt = sprintf('%s Time profiles are plotted in a logarithmic scale.',figtxt);
end


% legend % second entry Ref not needed for simulation
if nData>0
    legendEntries{3} = nameData;
end
if ~isnan(lloq)
     legendEntries{end+1} = 'lower limit of quantification';
end

return

function  [figtxt,legendEntries] = textVPCPredictedVsObserved(output,simulation,nameData,scale,lloq)


% get name and figure description
figtxt = sprintf('Predicted vs observed of %s for %s. Data source: %s.',output,simulation,nameData);
                
if ~isnan(lloq)
    figtxt = sprintf('%s Data below lower limit of quantification (lloq) are plotted as open symbols as lloq/2.',figtxt);
end

switch scale
    case 'lin'
        figtxt = sprintf('%s Prediction and observations are plotted in a linear scale.',figtxt);
    case 'log'
        figtxt = sprintf('%s Prediction and observations are plotted in a logarithmic scale.',figtxt);
end


% legend % second entry Ref not needed for simulation
legendEntries{1} = nameData;
legendEntries{2} = 'line of identity';
if ~isnan(lloq)
     legendEntries{end+1} = 'lower limit of quantification';
end

return

function [figtxt,legendEntries] = textVPCResidualsVsY(output,simulation,nameData,scale,lloq)

% get name and figure description
switch scale
    case 'lin'
        figtxt = sprintf('Linear residuals of %s vs predicted values for %s. Data source: %s.',output,simulation,nameData);
    case 'log'
        figtxt = sprintf('Logarithmic resiudals vs predicted values of %s for %s. Data source: %s.',output,simulation,nameData);
end

if ~isnan(lloq)
    figtxt = sprintf('%s For data below lower limit of quantification (lloq) the residuals were calculated using lloq/2 and as open symbols.',figtxt);
end


legendEntries{1} = nameData;

return

function [figtxt,legendEntries] = textVPCResidualsVsTime(output,simulation,nameData,scale,lloq)

% get name and figure description
switch scale
    case 'lin'
        figtxt = sprintf('Linear resiudals vs time of %s for %s. Data source: %s.',output,simulation,nameData);
    case 'log'
        figtxt = sprintf('Logarithmic resiudals vs time of %s for %s. Data source: %s.',output,simulation,nameData);
end

if ~isnan(lloq)
    figtxt = sprintf('%s For data below lower limit of quantification (lloq) the residuals were caluclated using lloq/2 and as open symbols.',figtxt);
end


legendEntries{1} = nameData;

return

function [figtxt,legendEntries] = textVPCResidualsAsHistogram(outputNameList,simulation)

figtxt = sprintf('Distribution of residuals for %s',simulation);


legendEntries = outputNameList;

return

function [figtxt,legendEntries] = textVPCResidualsAsQQPlot(outputNameList,simulation)

figtxt = sprintf('Residuals for %s as quantile-quantile plot.',simulation);


legendEntries = outputNameList;

return

