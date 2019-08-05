function plotQualificationTimeProfile(WSettings,plotIndex,TimeProfile,ObservedDataSets,SimulationMappings, Curves, AxesOptions, PlotSettings, REInputPath)
%PLOTQUALIFICATIONTIMEPROFILE Plots time profile from Configuration Plan
%
% plotQualificationTimeProfile(WSettings,plotIndex,
%   TimeProfile,ObservedDataSets,SimulationMappings, Curves, AxesOptions, PlotSettings, REInputPath)
%
% Inputs:
%   WSettings (structure)    definition of properties used in all
%                   workflow functions see GETDEFAULTWORKFLOWSETTINGS
%   plotIndex (integer) index of plot
%   TimeProfile (structure) TimeProfile plot information
%   Curves (structure) Curves information
%   ObservedDataSets (structure) Observed data
%   SimulationMappings (structure) Map simulation results to project
%   AxesOptions (structure) to set axes options
%   PlotSettings (structure) to set plot options
%   REInputPath (string) path of RE input files and folders
% Output
%

% Open Systems Pharmacology Suite;  http://open-systems-pharmacology.org

%---------------------------------------------

% Create figure with first setting from WSettings using getReportFigure
% To be updated using the Configuration plan Settings as optional arguments

ax = getReportFigureQP(WSettings,1,1,[],PlotSettings);

[xAxesOptions, yAxesOptions, yyAxesOptions] = setFigureOptions(AxesOptions);

% Load the mapped Time Profile Simulation Results
[csvSimFile, xmlfile] = getSimFile(TimeProfile, SimulationMappings, REInputPath);
if isempty(csvSimFile)
    ME = MException('plotQualificationTimeProfile:notFoundInPath', ...
        'In Time Profile Plot %d, Project "%s" or Simulation "%s" were not found in SimulationMappings', plotIndex, TimeProfile.Project, TimeProfile.Simulation);
    throw(ME);
end
SimResult = loadSimResultcsv(csvSimFile, TimeProfile);

if isempty(SimResult.outputPathList)
    ME = MException('plotQualificationTimeProfile:emptyOutputPathInSimulation', ...
        'In Time Profile Plot %d Project "%s" Simulation "%s", OutputPath is empty', plotIndex, TimeProfile.Project, TimeProfile.Simulation);
    throw(ME);
end

% Initialize simulation, and get Molecular Weight in g/mol for correct use of getUnitFactor
initSimulation(xmlfile,'none');

% Perform the plot based on Curves indications
curvesLegend={};
curvesHandle=[];

% Plot the Curves indicated by the array of structures Curves
for i=1:length(Curves)
    
    % Get Molecular Weight for Conversion
    MW = getMolecularWeightForPath(Curves(i).Y);
    
    % ObservedData type is indicated as 2nd element of Path
    OutputType = getElementsfromPath(Curves(i).Y);
    OutputType = OutputType{2};
    
    if strcmp(OutputType, 'ObservedData')
        % For observation: Get the right observation curve with right unit
        [p_handle, legLabel] = testandplotObservations(Curves(i), ObservedDataSets, MW, xAxesOptions, yAxesOptions, yyAxesOptions);
        % If the output was not found
        if isempty(p_handle)
            ME = MException('plotQualificationTimeProfile:notFoundInPath', ...
                'In Time Profile Plot %d, Curves %d : %s not found', plotIndex, i, Curves(i).Y);
            throw(ME);
        end
    else
        % For simulations: Get the right simulation curve
        [p_handle, legLabel] = testandplotSimResults(Curves(i), SimResult, MW, xAxesOptions, yAxesOptions, yyAxesOptions);
        % If the output was not found
        if isempty(p_handle)
            ME = MException('plotQualificationTimeProfile:notFoundInPath', ...
                'In Time Profile Plot %d, Curves %d : %s not found', plotIndex, i, Curves(i).Y);
            throw(ME);
        end
    end
    curvesLegend=[curvesLegend legLabel];
    curvesHandle=[curvesHandle p_handle];
end
legend(curvesHandle, curvesLegend);

% Scale X-Axis according to data if xlim is not provided
if ~isfield(xAxesOptions, 'Min') && ~isfield(xAxesOptions, 'Max')
    PlottedData = get(gca, 'children');
    Xmin = min(cellfun(@min, arrayfun(@(x) getfield(x, 'XData'), PlottedData, 'UniformOutput', false)));
    Xmax = max(cellfun(@max, arrayfun(@(x) getfield(x, 'XData'), PlottedData, 'UniformOutput', false)));
    xlim([Xmin Xmax]);
end

% ------------------------- Auxiliary functions -------------------------
% For simulations: Get the right simulation curve with right unit
function [p_handle, legendLabels] = testandplotSimResults(Curves, SimResult, MW, xAxesOptions, yAxesOptions, yyAxesOptions)

p_handle=[];

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
                        [p_handle, legendLabels] = plotPopulationStatistics(SimResult.time.*Xfactor, SimResult.y{j}.*Yfactor, Curves, yyAxesOptions.Scaling);
                        yyaxis left
                        break
                    end
                else
                    p_handle = plot(SimResult.time.*Xfactor, SimResult.y{j}.*Yfactor);
                    setCurveOptions(p_handle, Curves.CurveOptions);
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
                    [p_handle, legendLabels] = plotPopulationStatistics(SimResult.time.*Xfactor, SimResult.y{j}.*Yfactor, Curves, yAxesOptions.Scaling);
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

p_handle=[];
legendLabels={};

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
                
                % Convert the output to the correct unit
                ObservedTime = ObservedDataSets(j).time.*Xfactor;
                ObservedOutput = ObservedDataSets(j).y{1}.*Yfactor;
                
                % Plot the output
                p_handle = plot(ObservedTime, ObservedOutput);
                setCurveOptions(p_handle, Curves.CurveOptions);
                
                % Check if error bars to be plotted
                if length(ObservedDataSets(j).outputPathList)>1
                    if strcmp('Fraction', ObservedDataSets(j).outputDimension{2})
                        % Geometric SD is assumed for no dimension unit if
                        % SD>=1 else it is arithmetic
                        if min(ObservedDataSets(j).y{2})>1
                            p_handle2=errorbar(ObservedTime, ObservedOutput, ...
                                ObservedOutput.*(1-1./ObservedDataSets(j).y{2}), ObservedOutput.*(ObservedDataSets(j).y{2}-1), 'HandleVisibility','off');
                            setCurveOptions(p_handle2, p_handle);
                        else
                            p_handle2=errorbar(ObservedTime, ObservedOutput, ObservedDataSets(j).y{2}, 'HandleVisibility','off');
                            setCurveOptions(p_handle2, p_handle);
                        end
                    else
                        errorfactor=getUnitFactor(ObservedDataSets(j).outputUnit{2},yyAxesOptions.Unit,YDimension, 'MW',MW);
                        p_handle2=errorbar(ObservedTime, ObservedOutput, ObservedDataSets(j).y{2}.*errorfactor, 'HandleVisibility','off');
                        setCurveOptions(p_handle2, p_handle);
                    end
                end
                % Check for LLOQ to be plotted
                if ~isempty(ObservedDataSets(j).LLOQ)
                    LLOQfactor = getUnitFactor(ObservedDataSets(j).LLOQUnit,yyAxesOptions.Unit,YDimension, 'MW',MW);
                    yLLOQ = LLOQfactor.*ObservedDataSets(j).LLOQ*[1 1];
                    xLLOQ = get(gca, 'xlim');
                    LLOQ_handle = plot(xLLOQ, yLLOQ, 'k-', 'HandleVisibility','off');
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
            
            % Convert the output to the correct unit
            ObservedTime = ObservedDataSets(j).time.*Xfactor;
            ObservedOutput = ObservedDataSets(j).y{1}.*Yfactor;
            
            % Plot the output
            p_handle = plot(ObservedTime, ObservedOutput);
            setCurveOptions(p_handle, Curves.CurveOptions);
            
            % Check if error bars to be plotted
            if length(ObservedDataSets(j).outputPathList)>1
                if strcmp('Fraction', ObservedDataSets(j).outputDimension{2})
                    % Geometric SD is assumed for no dimension unit if
                    % SD>=1 else it is arithmetic
                    if min(ObservedDataSets(j).y{2})>1
                        p_handle2=errorbar(ObservedTime, ObservedOutput, ...
                            ObservedOutput.*(1-1./ObservedDataSets(j).y{2}), ObservedOutput.*(ObservedDataSets(j).y{2}-1), 'HandleVisibility','off');
                        setCurveOptions(p_handle2, p_handle);
                    else
                        p_handle2=errorbar(ObservedTime, ObservedOutput, ObservedDataSets(j).y{2}, 'HandleVisibility','off');
                        setCurveOptions(p_handle2, p_handle);
                    end
                else
                    errorfactor=getUnitFactor(ObservedDataSets(j).outputUnit{2},yAxesOptions.Unit,YDimension, 'MW',MW);
                    p_handle2=errorbar(ObservedTime, ObservedOutput, ObservedDataSets(j).y{2}.*errorfactor, 'HandleVisibility','off');
                    setCurveOptions(p_handle2, p_handle);
                end
            end
            % Check for LLOQ to be plotted
            if ~isempty(ObservedDataSets(j).LLOQ)
                LLOQfactor = getUnitFactor(ObservedDataSets(j).LLOQUnit,yAxesOptions.Unit,YDimension, 'MW',MW);
                yLLOQ = LLOQfactor.*ObservedDataSets(j).LLOQ*[1 1];
                xLLOQ = get(gca, 'xlim');
                LLOQ_handle = plot(xLLOQ, yLLOQ, 'k-', 'HandleVisibility','off');
            end
            
        end
        break
        
    end
    
end


function [p_handle, legendLabels] = plotPopulationStatistics(time, Y, Curves, Scaling)
p_handle=[];

time = reshape(time, 1, []);
ll = size(Y,1);
if ll == length(time)
    Y=Y';
end

for i=1:length(Curves.Statistics)
    if strcmp(Curves.Statistics(i).Id, 'ArithmeticMean')
        p_handle(i) = plot(time, mean(Y));
        legendLabels{i}=sprintf('%s-Arithmetic Mean', Curves.Y);
        setCurveOptions(p_handle(i), Curves.Statistics(i));
    end
    if strcmp(Curves.Statistics(i).Id, 'ArithmeticStandardDeviation')
        p_handle(i) = plot([time NaN time], [mean(Y)-std(Y) NaN mean(Y)+std(Y)]);
        legendLabels{i}=sprintf('%s-Arithmetic Standard Deviation', Curves.Y);
        setCurveOptions(p_handle(i), Curves.Statistics(i));
    end
    if strcmp(Curves.Statistics(i).Id, 'GeometricMean')
        p_handle(i) = plot(time, geomean(Y));
        legendLabels{i}=sprintf('%s-Geometric Mean', Curves.Y);
        setCurveOptions(p_handle(i), Curves.Statistics(i));
    end
    if strcmp(Curves.Statistics(i).Id, 'GeometricStandardDeviation')
        p_handle(i) = plot([time NaN time], [exp(mean(log(Y))-std(log(Y))) NaN exp(mean(log(Y))+std(log(Y)))]);
        legendLabels{i}=sprintf('%s-Geometric Standard Deviation', Curves.Y);
        setCurveOptions(p_handle(i) , Curves.Statistics(i));
    end
    if strcmp(Curves.Statistics(i).Id, 'Median')
        p_handle(i) = plot(time, median(Y));
        legendLabels{i}=sprintf('%s-Median', Curves.Y);
        setCurveOptions(p_handle(i), Curves.Statistics(i));
    end
    if strcmp(Curves.Statistics(i).Id, 'Min')
        p_handle(i) = plot(time, min(Y));
        legendLabels{i}=sprintf('%s-Min', Curves.Y);
        setCurveOptions(p_handle(i), Curves.Statistics(i));
    end
    if strcmp(Curves.Statistics(i).Id, 'Max')
        p_handle(i) = plot(time, max(Y));
        legendLabels{i}=sprintf('%s-Max', Curves.Y);
        setCurveOptions(p_handle(i), Curves.Statistics(i));
    end
    if contains(Curves.Statistics(i).Id, 'Percentile')
        perc=sscanf(Curves.Statistics(i).Id, 'Percentile_%d')/100;
        p_handle(i) = plot(time, quantile(Y, perc));
        legendLabels{i}=sprintf('%s-Percentile %d%%', Curves.Y, perc*100);
        setCurveOptions(p_handle(i), Curves.Statistics(i));
    end
    if contains(Curves.Statistics(i).Id, 'Range')
        perc = sscanf(Curves.Statistics(i).Id, 'Range%d')/100;
        ran = quantile(Y, [(1-perc)/2 (1+perc)/2]);
        if strcmpi(Scaling, 'Log')
            % Remove 0s from plot in log
            zeros2remove = min(ran==0);
            time(zeros2remove) = [];
            ran(:,zeros2remove) = [];
        end
        p_handle(i) = patch([time time(end:-1:1)], [ran(1,:) ran(2,end:-1:1)], [perc perc perc]);
        legendLabels{i}=sprintf('%s-Range %d%% to %d%%', Curves.Y, round(100*(1-perc)/2), round(100*(1+perc)/2));
        setCurveOptions(p_handle(i), Curves.Statistics(i));
        set(p_handle(i), 'FaceAlpha', 0.5);
    end
end
