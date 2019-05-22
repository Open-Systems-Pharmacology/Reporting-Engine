function plotQualificationTimeProfile(WSettings,figureHandle,TimeProfile,ObservedDataSets,SimulationMappings, Curves, AxesOptions, PlotSettings, REInputPath)
%PLOTQUALIFICATIONTIMEPROFILE Plots the time profile of a population in comparison to a reference population
%
% plotQualificationTimeProfile(WSettings,figureHandle,SimTL,DataTP, Curves, AxesOptions,PlotSettings)
%
% Inputs:
%       WSettings (structure)    definition of properties used in all
%                   workflow functions see GETDEFAULTWORKFLOWSETTINGS
%   figureHandle ( handle) handle of figure
%   TimeProfile (structure) with Project and Simulation to be mapped
%   ObservedDataSets (structure) of all loaded observed data
%   SimulationMappings (structure) to mapping linking simulations to there path
%   Curves (structure) to internal path of simulations or observations
%   AxesOptions (structure) to set plot options
%   PlotSettings (structure) to set plot options
%
% Open Systems Pharmacology Suite;  http://open-systems-pharmacology.org

%---------------------------------------------

% Create figure with first setting from WSettings using getReportFigure
% To be updated using the Configuration plan Settings as optional arguments

ax = getReportFigureQP(WSettings,1,1,figureHandle,PlotSettings);

[xAxesOptions, yAxesOptions, yyAxesOptions] = setFigureOptions(AxesOptions);

% Load the mapped Time Profile Simulation Results
[csvSimFile, xmlfile] = getSimFile(TimeProfile, SimulationMappings, REInputPath);
if isempty(csvSimFile)
    ME = MException('plotQualificationTimeProfile:notFoundInPath', ...
        'In Time Profile plot %d, Project "%s" or Simulation "%s" were not found in SimulationMappings', figureHandle, TimeProfile.Project, TimeProfile.Simulation);
    throw(ME);
end
SimResult = loadSimResultcsv(csvSimFile, TimeProfile);

% Initialize simulation, and get Molecular Weight in g/mol for correct use of getUnitFactor
initSimulation(xmlfile,'none');

% Get Molecular Weight for Conversion
for i=1:length(Curves)
    OutputType = getElementsfromPath(Curves(i).Y);
    OutputType = OutputType{2};
    if strcmp(OutputType, 'ObservedData')
        MW = getMolecularWeightForPath(Curves(i).Y);
    end
    if exist('MW')
        break
    end
end


% Perform the plot based on Curves indications
legendLabels={};
hold on;

% Plot the Curves indicated by the array of structures Curves
for i=1:length(Curves)
    
    % ObservedData type is indicated as 2nd element of Path
    OutputType = getElementsfromPath(Curves(i).Y);
    OutputType = OutputType{2};
    
    if strcmp(OutputType, 'ObservedData')
        % For observation: Get the right observation curve with right unit
        [p_handle, legLabel] = testandplotObservations(Curves(i), ObservedDataSets, MW, xAxesOptions, yAxesOptions, yyAxesOptions);
        % If the output was not found
        if isempty(p_handle)
            ME = MException('plotQualificationTimeProfile:notFoundInPath', ...
                'In plot %d, Curves %d : %s not found', figureHandle, i, Curves(i).Y);
            throw(ME);
        end
    else
        % For simulations: Get the right simulation curve
        [p_handle, legLabel] = testandplotSimResults(Curves(i), SimResult, MW, xAxesOptions, yAxesOptions, yyAxesOptions);
        % If the output was not found
        if isempty(p_handle)
            ME = MException('plotQualificationTimeProfile:notFoundInPath', ...
                'In plot %d, Curves %d : %s not found', figureHandle, i, Curves(i).Y);
            throw(ME);
        end
    end
    legendLabels=[legendLabels legLabel];
    
end
legend(legendLabels, 'Location', 'northoutside');
%legend('off')

% ------------------------- Auxiliary functions -------------------------
% For simulations: Get the right simulation curve with right unit
function [p_handle, legendLabels] = testandplotSimResults(Curves, SimResult, MW, xAxesOptions, yAxesOptions, yyAxesOptions)

p_handle=[];
legendLabels={};

for j = 1:length(SimResult.outputPathList)
    
    findPathOutput = contains(Curves.Y,SimResult.outputPathList{j});
    if findPathOutput
        legendLabels=Curves.Name;
        
        % For right Y axis
        if isfield(Curves.CurveOptions, 'yAxisType')
            if strcmp(Curves.CurveOptions.yAxisType, 'Y2')
                yyaxis right
                % Get the dimension from the Unit
                % A flag can be output if the unit was not found
                YDimension = findDimensionfromUnit(yyAxesOptions.Unit);
                Yfactor=getUnitFactor(SimResult.outputUnit{j},yyAxesOptions.Unit,YDimension, 'MW',MW);
                
                % Convert units to reference unit
                XDimension = findDimensionfromUnit(xAxesOptions.Unit);
                Xfactor=getUnitFactor(SimResult.timeUnit,xAxesOptions.Unit,XDimension);
                
                if isfield(Curves, 'Type')
                    if strcmp(Curves.Type, 'Population')
                        [p_handle, legendLabels] = plotPopulationStatistics(SimResult.time.*Xfactor, SimResult.y{j}.*Yfactor, Curves);
                    end
                else
                    p_handle = plot(SimResult.time.*Xfactor, SimResult.y{j}.*Yfactor);
                    yyaxis left
                    break
                end
            end
        else
            % For left Y axis
            
            YDimension = findDimensionfromUnit(yAxesOptions.Unit);
            Yfactor=getUnitFactor(SimResult.outputUnit{j},yAxesOptions.Unit,YDimension, 'MW',MW);
            
            % Convert units to reference unit
            XDimension = findDimensionfromUnit(xAxesOptions.Unit);
            Xfactor=getUnitFactor(SimResult.timeUnit,xAxesOptions.Unit,XDimension);
            if isfield(Curves, 'Type')
                if strcmp(Curves.Type, 'Population')
                    [p_handle, legendLabels] = plotPopulationStatistics(SimResult.time.*Xfactor, SimResult.y{j}.*Yfactor, Curves);
                end
            else
                p_handle = plot(SimResult.time.*Xfactor, SimResult.y{j}.*Yfactor);
                setCurveOptions(p_handle, Curves.CurveOptions);
                break
            end
        end
    end
end

% For observation: Get the right observation curve with right unit
function [p_handle, legendLabels] = testandplotObservations(Curves, ObservedDataSets, MW, xAxesOptions, yAxesOptions, yyAxesOptions)

% Split the Curve path into its elements
CurveElements = getElementsfromPath(Curves.Y);

for j = 1:length(ObservedDataSets)
    findPathOutput = strcmp(CurveElements{1},ObservedDataSets(j).Id);
    if findPathOutput
        % Get the legend
        legendLabels=Curves.Name;
        
        % For right Y axis
        if isfield(Curves.CurveOptions, 'yAxisType')
            if strcmp(Curves.CurveOptions.yAxisType, 'Y2')
                yyaxis right
                % Get the dimension from the Unit
                % A flag can be output if the unit was not found
                YDimension = findDimensionfromUnit(yyAxesOptions.Unit);
                Yfactor=getUnitFactor(ObservedDataSets(j).outputUnit{1},yyAxesOptions.Unit,YDimension, 'MW',MW);
                
                % Convert units to reference unit
                XDimension = findDimensionfromUnit(xAxesOptions.Unit);
                Xfactor=getUnitFactor(ObservedDataSets(j).timeUnit,xAxesOptions.Unit,XDimension);
                
                % Check if the observation names match (many warning message
                % appears and are left aside so far)
                %{
            if ~strcmp(ObservedDataSets(j).outputPathList{1}, CurveElements{end})
                writeToReportLog('WARNING', sprintf('Warning: Curve %s in TimeProfile plot \n Curve path does not match ObservedData Path. \n Curve path: %s \n ObservedData path: %s \n',...
                    Curves.Name, CurveElements{end}, ObservedDataSets(j).outputPathList{1}));
            end
                %}
                
                % Convert the output to the correct unit
                ObservedTime = ObservedDataSets(j).time.*Xfactor;
                ObservedOutput = ObservedDataSets(j).y{1}.*Yfactor;
                
                % Plot the output
                p_handle = plot(ObservedTime, ObservedOutput);
                
                % Check if error bars to be plotted
                if length(ObservedDataSets(j).outputPathList)>1
                    if strcmp('Fraction', ObservedDataSets(j).outputDimension{2})
                        % Geometric SD is assumed for no dimension unit
                        p_handle2=errorbar(ObservedTime, ObservedOutput, ...
                            ObservedOutput - ObservedOutput./ObservedDataSets(j).y{2}, ObservedOutput.*ObservedDataSets(j).y{2} - ObservedOutput);
                        setCurveOptions(p_handle2, Curves.CurveOptions);
                    else
                        errorfactor=getUnitFactor(ObservedDataSets(j).outputUnit{2},yAxesOptions.Unit,YDimension, 'MW',MW);
                        p_handle2=errorbar(ObservedTime, ObservedOutput, ObservedDataSets(j).y{2}.*errorfactor);
                        setCurveOptions(p_handle2, Curves.CurveOptions);
                    end
                    
                end
                
                yyaxis left
                break
            end
        else
            % For left Y axis
            
            % A flag can be output if the unit was not found
            YDimension = findDimensionfromUnit(yAxesOptions.Unit);
            Yfactor=getUnitFactor(ObservedDataSets(j).outputUnit{1},yAxesOptions.Unit,YDimension, 'MW',MW);
            
            % Convert units to reference unit
            XDimension = findDimensionfromUnit(xAxesOptions.Unit);
            Xfactor=getUnitFactor(ObservedDataSets(j).timeUnit,xAxesOptions.Unit,XDimension);
            
            % Check if the observation names match (many warning message
            % appears and are left aside so far)
            %{
            if ~strcmp(ObservedDataSets(j).outputPathList{1}, CurveElements{end})
                writeToReportLog('WARNING', sprintf('Warning: Curve %s in TimeProfile plot \n Curve path does not match ObservedData Path. \n Curve path: %s \n ObservedData path: %s \n',...
                    Curves.Name, CurveElements{end}, ObservedDataSets(j).outputPathList{1}));
            end
            %}
            % Convert the output to the correct unit
            ObservedTime = ObservedDataSets(j).time.*Xfactor;
            ObservedOutput = ObservedDataSets(j).y{1}.*Yfactor;
            
            % Plot the output
            p_handle = plot(ObservedTime, ObservedOutput);
            
            
            % Check if error bars to be plotted
            if length(ObservedDataSets(j).outputPathList)>1
                if strcmp('Fraction', ObservedDataSets(j).outputDimension{2})
                    % Geometric SD is assumed for fraction dimension
                    p_handle2=errorbar(ObservedTime, ObservedOutput, ...
                        ObservedOutput - ObservedOutput./ObservedDataSets(j).y{2}, ObservedOutput.*ObservedDataSets(j).y{2} - ObservedOutput);
                    setCurveOptions(p_handle2, Curves.CurveOptions);
                else
                    errorfactor=getUnitFactor(ObservedDataSets(j).outputUnit{2},yAxesOptions.Unit,YDimension, 'MW',MW);
                    p_handle2=errorbar(ObservedTime, ObservedOutput, ObservedDataSets(j).y{2}.*errorfactor);
                    setCurveOptions(p_handle2, Curves.CurveOptions);
                end
                
            end
            
        end
        break
        
    end
    
end
if exist('p_handle')
    setCurveOptions(p_handle, Curves.CurveOptions);
else
    p_handle=[];
    legendLabels=[];
end

function [pp, legendLabels] = plotPopulationStatistics(time, Y, Curves)
pp=[];
legendLabels={};

time = reshape(time, 1, []);
ll = size(Y,1);
if ll == length(time)
    Y=Y';
end

for i=1:length(Curves.Statistics)
    if strcmp(Curves.Statistics(i).Id, 'ArithmeticMean')
        pp = plot(time, mean(Y));
        legendLabels{length(legendLabels)+1}=sprintf('Arithmetic Mean %s', Curves.Y);
        setCurveOptions(pp, Curves.Statistics(i));
    end
    if strcmp(Curves.Statistics(i).Id, 'ArithmeticStandardDeviation')
        pp = plot([time NaN time], [mean(Y)-std(Y) NaN mean(Y)+std(Y)]);
        legendLabels{length(legendLabels)+1}=sprintf('Arithmetic Standard Deviation %s', Curves.Y);
        setCurveOptions(pp, Curves.Statistics(i));
    end
    if strcmp(Curves.Statistics(i).Id, 'GeometricMean')
        pp = plot(time, geomean(Y));
        legendLabels{length(legendLabels)+1}=sprintf('Geometric Mean %s', Curves.Y);
        setCurveOptions(pp, Curves.Statistics(i));
    end
    if strcmp(Curves.Statistics(i).Id, 'GeometricStandardDeviation')
        pp = plot([time NaN time], [exp(mean(log(Y))-std(log(Y))) NaN exp(mean(log(Y))+std(log(Y)))]);
        legendLabels{length(legendLabels)+1}=sprintf('Geometric Standard Deviation %s', Curves.Y);
        setCurveOptions(pp , Curves.Statistics(i));
    end
    if strcmp(Curves.Statistics(i).Id, 'Median')
        pp = plot(time, median(Y));
        legendLabels{length(legendLabels)+1}=sprintf('Median %s', Curves.Y);
        setCurveOptions(pp, Curves.Statistics(i));
    end
    if strcmp(Curves.Statistics(i).Id, 'Min')
        pp = plot(time, min(Y));
        legendLabels{length(legendLabels)+1}=sprintf('Min %s', Curves.Y);
        setCurveOptions(pp, Curves.Statistics(i));
    end
    if strcmp(Curves.Statistics(i).Id, 'Max')
        pp = plot(time, max(Y));
        legendLabels{length(legendLabels)+1}=sprintf('Max %s', Curves.Y);
        setCurveOptions(pp, Curves.Statistics(i));
    end
    if contains(Curves.Statistics(i).Id, 'Percentile')
        perc=sscanf(Curves.Statistics(i).Id, 'Percentile_%d')/100;
        pp = plot(time, quantile(Y, perc));
        legendLabels{length(legendLabels)+1}=sprintf('Percentile %d %s ', perc*100, Curves.Y);
        setCurveOptions(pp, Curves.Statistics(i));
    end
    if contains(Curves.Statistics(i).Id, 'Range')
        perc = sscanf(Curves.Statistics(i).Id, 'Range%d')/100;
        ran = quantile(Y, [(1-perc)/2 (1+perc)/2]);
        pp = patch([time time(end:-1:1)], [ran(1,:) ran(2,end:-1:1)], [perc perc perc]);
        legendLabels{length(legendLabels)+1}=sprintf('Range %d %s ', perc*100, Curves.Y);
        setCurveOptions(pp, Curves.Statistics(i));
    end
end
