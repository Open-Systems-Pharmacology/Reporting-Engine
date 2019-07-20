function  [fig_handle, GMFE]=plotQualificationGOFMerged(WSettings,plotIndex,Groups,ObservedDataSets,SimulationMappings, AxesOptions, PlotSettings, REInputPath)
%PLOTQUALIFICATIONGOFMERGED Plots Goodness of fit from Configuration Plan
%
% [fig_handle, GMFE]=plotQualificationGOFMerged(WSettings,plotIndex,
%   Groups,ObservedDataSets,SimulationMappings, AxesOptions, PlotSettings, REInputPath)
%
% Inputs:
%   WSettings (structure)    definition of properties used in all
%                   workflow functions see GETDEFAULTWORKFLOWSETTINGS
%   plotIndex (integer) index of plot
%   Groups (structure) Goodness of fit plot information
%   ObservedDataSets (structure) Observed data
%   SimulationMappings (structure) Map simulation results to project
%   AxesOptions (structure) to set axes options
%   PlotSettings (structure) to set plot options
%   REInputPath (string) path of RE input files and folders
% Output
%   fig_handle (handle) handle of output figures
%   GMFE (cells)  Geometric Mean Fold Error
%

% Open Systems Pharmacology Suite;  http://open-systems-pharmacology.org

%---------------------------------------------

% Initialize the unity line and the time line of the plot
minX=NaN; maxX=NaN; minY=NaN; maxY=NaN; maxtime=NaN; maxRes=NaN;

% create figure for Obs vs Pred
[ax1, fig_handle.PredictedVsObserved] = getReportFigureQP(WSettings,1,1,[],PlotSettings);
[xAxesOptions, yAxesOptions] = setFigureOptions(AxesOptions.GOFMergedPlotsPredictedVsObserved);

% create figure for Residuals vs time
[ax2, fig_handle.ResidualsOverTime] = getReportFigureQP(WSettings,1,1,[],PlotSettings);
[TimeAxesOptions, ResAxesOptions] = setFigureOptions(AxesOptions.GOFMergedPlotsResidualsOverTime);

% Map for each group the corresponding observations within a same structure
for i=1:length(Groups)
    % Load simulation according to mapping
    for j=1:length(Groups(i).OutputMappings)
        
        % Load the mapped GOF Simulation Results
        Simulations = Groups(i).OutputMappings(j);
        [csvSimFile, xmlfile] = getSimFile(Simulations, SimulationMappings, REInputPath);
        if isempty(csvSimFile)
            ME = MException('plotQualificationGOFMerged:notFoundInPath', ...
                'In GOF Merged Plot %d, Group %d, mapping %d, Project "%s" or Simulation "%s" were not found in SimulationMappings', plotIndex, i, j, Simulations.Project, Simulations.Simulation);
            throw(ME);
        end
        SimResult = loadSimResultcsv(csvSimFile, Simulations);
        
        if isempty(SimResult.outputPathList)
            ME = MException('plotQualificationGOFMerged:emptyOutputPathInSimulation', ...
                'In GOF Merged Plot %d, Group %d, mapping %d, OutputPath is empty in Project "%s" Simulation "%s"', plotIndex, i, j, Simulations.Project, Simulations.Simulation);
            throw(ME);
        end
        
        % Initialize simulation, and get Molecular Weight in g/mol for correct use of getUnitFactor
        initSimulation(xmlfile,'none');
        % Get Molecular Weight for Conversion
        MW = getMolecularWeightForPath(Simulations.Output);
        
        % Get the right simulation output to be compared
        [predictedTime, predicted] = testSimResults(Simulations, SimResult, MW, TimeAxesOptions, yAxesOptions);
        
        if isempty(predicted)
            ME = MException('plotQualificationGOFMerged:notFoundInPath', ...
                'In GOF Merged Plot %d, Group %d SubGroup %d : %s not found', plotIndex, i, j, Simulations.Output);
            throw(ME);
        end
        
        % Get the right simulation output to be compared
        [ObsTime, Obs] = testObservations(Simulations, ObservedDataSets, MW, TimeAxesOptions, yAxesOptions);
        
        if isempty(Obs)
            ME = MException('plotQualificationGOFMerged:notFoundInPath', ...
                'In GOF Merged Plot %d, Group %d SubGroup %d : %s not found', plotIndex, i, j, Simulations.ObservedData);
            throw(ME);
        end
        
        % Get points comparable between obs and pred
        comparable_index= getComparablePredicted(ObsTime , predictedTime);
        
        % Group and save comparable points
        Group(i).dataTP(j).yobs = reshape(Obs, [], 1);
        Group(i).dataTP(j).ypred = reshape(predicted(comparable_index), [], 1);
        Yres = predicted(comparable_index)-Obs;
        if strcmp(ResAxesOptions.Dimension,'Fraction')
            Yres=Yres./Obs;
        else
            Resfactor=getUnitFactor(yAxesOptions.Unit,ResAxesOptions.Unit,ResAxesOptions.Dimension, 'MW',MW);
            Yres=Yres.*Resfactor;
        end
        Group(i).dataTP(j).yres= reshape(Yres, [], 1);
        Group(i).dataTP(j).time = reshape(ObsTime, [], 1);
        
        % Set the min/max of the axis
        PosObs = Group(i).dataTP(j).yobs>0;
        PosPred = Group(i).dataTP(j).ypred>0;
        minX=nanmin(nanmin(Group(i).dataTP(j).yobs(PosObs)), minX);
        maxX=nanmax(nanmax(Group(i).dataTP(j).yobs(PosObs)), maxX);
        minY=nanmin(nanmin(Group(i).dataTP(j).ypred(PosPred)), minY);
        maxY=nanmax(nanmax(Group(i).dataTP(j).ypred(PosPred)), maxY);
        maxtime=nanmax(nanmax(ObsTime), maxtime);
        maxRes=nanmax(nanmax(abs(Yres)), maxRes);
        
        
    end
end

% -------------------------------------------------------------
% Plot the figures
% Observation vs Prediction Figure
set(0, 'CurrentFigure', fig_handle.PredictedVsObserved);
plot([0 0.8*nanmin(minX, minY) 1.2*nanmax(maxX, maxY)], [0 0.8*nanmin(minX, minY) 1.2*nanmax(maxX, maxY)], '--k', 'Linewidth', 1,'HandleVisibility','off');
axis([0.8*nanmin(minX, minY) 1.2*nanmax(maxX, maxY) 0.8*nanmin(minX, minY) 1.2*nanmax(maxX, maxY)]);
legendLabels={};

% Initialize error for computing GMFE
Error=[];

for i=1:length(Groups)
    
    for j=1:length(Groups(i).OutputMappings)
        CurveOptions.Color=Groups(i).OutputMappings(j).Color;
        CurveOptions.Symbol=Groups(i).Symbol;
        CurveOptions.LineStyle='none';
        
        % One legend caption per group
        if j>1
            pp=plot(Group(i).dataTP(j).yobs, Group(i).dataTP(j).ypred, 'HandleVisibility','off');
        else
            pp=plot(Group(i).dataTP(j).yobs, Group(i).dataTP(j).ypred);
        end
        setCurveOptions(pp, CurveOptions);
        PosObs = Group(i).dataTP(j).yobs>0;
        PosPred = Group(i).dataTP(j).ypred>0;
        Error = [Error; log10(Group(i).dataTP(j).ypred(PosObs & PosPred))-log10(Group(i).dataTP(j).yobs((PosObs & PosPred)))];
    end
    legendLabels=[legendLabels Groups(i).Caption];
end

xLabelFinal = getLabelWithUnit(sprintf('Observed %s', findDimensionfromUnit(xAxesOptions.Unit)),xAxesOptions.Unit);
yLabelFinal = getLabelWithUnit(sprintf('Simulated %s', findDimensionfromUnit(yAxesOptions.Unit)),yAxesOptions.Unit);
xlabel(xLabelFinal); ylabel(yLabelFinal);

GMFE = 10.^(sum(abs(Error))/length(Error));

if ~isempty(legendLabels)
    lgd = legend(legendLabels, 'Location', 'northoutside');
    reshapeQualificationLegend(lgd);
else
    legend('off');
end

% Residuals vs Time Figure
% create figure for Residuals vs time
set(0, 'CurrentFigure', fig_handle.ResidualsOverTime);
plot([0 1.2*maxtime], [0 0], '--k', 'Linewidth', 1, 'HandleVisibility','off');
axis([0 1.2*maxtime -1.2*abs(maxRes) 1.2*abs(maxRes)]);

legendLabels={};

for i=1:length(Groups)
    for j=1:length(Groups(i).OutputMappings)
        CurveOptions.Color=Groups(i).OutputMappings(j).Color;
        CurveOptions.Symbol=Groups(i).Symbol;
        CurveOptions.LineStyle='none';
        
        % One legend caption per group
        if j>1
            pp=plot(Group(i).dataTP(j).time, Group(i).dataTP(j).yres, 'HandleVisibility','off');
        else
            pp=plot(Group(i).dataTP(j).time, Group(i).dataTP(j).yres);
        end
        setCurveOptions(pp, CurveOptions);
    end
    legendLabels=[legendLabels Groups(i).Caption];
end
xLabelFinal = getLabelWithUnit('Time',TimeAxesOptions.Unit);
yLabelFinal = getLabelWithUnit('Residuals',ResAxesOptions.Unit);
xlabel(xLabelFinal); ylabel(yLabelFinal);

if ~isempty(legendLabels)
    lgd = legend(legendLabels, 'Location', 'northoutside');
    reshapeQualificationLegend(lgd);
else
    legend('off');
end

% ---------------- Auxiliary function ------------------------------------

% get comparable time points between observations and simulations
function Index = getComparablePredicted(timeObs, timePred)

timeObs = reshape(timeObs, 1, []);
timePred = reshape(timePred, [], 1);

timeObs2 = repmat(timeObs, length(timePred), 1);
timePred2 = repmat(timePred, 1, length(timeObs));

error = (timeObs2-timePred2).*(timeObs2-timePred2);

[Y, Index] = min(error);


% For simulations: Get the right simulation curve with right unit
function [predictedTime, predicted] = testSimResults(Simulations, SimResult, MW, TimeAxesOptions, yAxesOptions)

% If no path is matched ObsTime will be emtpy
predictedTime=[];
predicted=[];

% Get dimensions for scaling plots
YDimension=findDimensionfromUnit(yAxesOptions.Unit);

TimeDimension=findDimensionfromUnit(TimeAxesOptions.Unit);

for j = 1:length(SimResult.outputPathList)
    
    findPathOutput = contains(SimResult.outputPathList{j}, Simulations.Output);
    
    if findPathOutput
        
        % Get the data from simulations
        predictedTime=SimResult.time;
        predicted=SimResult.y{j};
        
        % Get unit factor to convert to reference unit
        Yfactor=getUnitFactor(SimResult.outputUnit{j},yAxesOptions.Unit,YDimension, 'MW',MW);
        Timefactor=getUnitFactor(SimResult.timeUnit,TimeAxesOptions.Unit,TimeDimension);
        
        % Convert units to
        predicted = predicted.*Yfactor;
        predictedTime = predictedTime.*Timefactor;
        break
    end
end

% For observation: Get the right observation curve with right unit
function [ObsTime, Obs] = testObservations(Simulations, ObservedDataSets, MW, TimeAxesOptions, yAxesOptions)

% If no path is matched ObsTime will be emtpy
ObsTime=[];
Obs=[];

YDimension=findDimensionfromUnit(yAxesOptions.Unit);

TimeDimension=findDimensionfromUnit(TimeAxesOptions.Unit);

for j = 1:length(ObservedDataSets)
    
    % Get the right observed data
    findPathOutput = strcmp(ObservedDataSets(j).Id, Simulations.ObservedData);
    if findPathOutput
        
        % Get the first outputPathList, may be modified if more than one
        % Observed data is in the file
        Obs=ObservedDataSets(j).y{1};
        ObsTime=ObservedDataSets(j).time;
        
        % Get unit factor to convert to reference unit
        Yfactor=getUnitFactor(ObservedDataSets(j).outputUnit{1},yAxesOptions.Unit,YDimension, 'MW',MW);
        Timefactor=getUnitFactor(ObservedDataSets(j).timeUnit,TimeAxesOptions.Unit,TimeDimension);
        
        % Convert units to
        Obs = Obs.*Yfactor;
        ObsTime = ObsTime.*Timefactor;
        break
    end
end
