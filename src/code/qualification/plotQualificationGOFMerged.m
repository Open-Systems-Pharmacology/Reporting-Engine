function  GMFE=plotQualificationGOFMerged(WSettings,figureHandle,Groups,ObservedDataSets,SimulationMappings, AxesOptions, PlotSettings)
% PLOTQUALIFICATIONGOFMERGED plot predcited vs observed and residuals vs time
%
% GMFE=plotQualificationGOFMerged(WSettings,figureHandle,Groups,ObservedDataSets,SimulationMappings, PlotSettings)
%
% Inputs:
%       WSettings (structure)    definition of properties used in all
%                   workflow functions see GETDEFAULTWORKFLOWSETTINGS
%   figureHandle ( handle) handle of figure
%   Groups (structure) with simulated results
%   ObservedDataSets (structure) of all loaded observed data
%   SimulationMappings (structure) to map simulations with observations
%   AxesOptions (structure) to set plot options
%   PlotSettings (structure) to set plot options
%
% Output
%   GMFE (numeric) measure of goodness of fit

% Open Systems Pharmacology Suite;  http://open-systems-pharmacology.org

% Initialize the unity line and the time line of the plot
minX=NaN; maxX=NaN; minY=NaN; maxY=NaN; maxtime=NaN;

% create figure for Obs vs Pred
[ax, fig_handle1] = getReportFigureQP(WSettings,1,1,figureHandle,PlotSettings);
[xAxesOptions, yAxesOptions] = setFigureOptions(AxesOptions.GOFMergedPlotsPredictedVsObserved);

% create figure for Residuals vs time
[ax, fig_handle2] = getReportFigureQP(WSettings,1,1,figureHandle+1,PlotSettings);
[TimeAxesOptions, ResAxesOptions] = setFigureOptions(AxesOptions.GOFMergedPlotsResidualsOverTime);

% Get dimensions for scaling plots
dimensionList=getDimensions;
if ~strcmp(yAxesOptions.Dimension, 'Dimensionless')
    YDimension=dimensionList{strContains(yAxesOptions.Dimension, dimensionList)};
else
    YDimension='';
end
if ~strcmp(ResAxesOptions.Dimension, 'Dimensionless')
    ResDimension=dimensionList{strContains(ResAxesOptions.Dimension, dimensionList)};
else
    ResDimension='';
end
TimeDimension=dimensionList{strContains(TimeAxesOptions.Dimension, dimensionList)};

% Map for each group the corresponding observations within a same structure
for i=1:length(Groups)
    % Load simulation according to mapping
    for j=1:length(Groups(i))
        Simulations = Groups(i).OutputMappings(j);
        [csvSimFile, xmlfile] = getSimFile(Simulations, SimulationMappings);
        SimResult = loadSimResultcsv(csvSimFile, Simulations);
        
        initSimulation(xmlfile,'none');
        Compound=Simulations.Project;
        MW = getParameter(sprintf('*|%s|Molecular weight',Compound),1,'parametertype','readonly');
        
        % Get the right simulation output to be compared
        for k=1:length(SimResult.outputPathList)
            % Find first integer of string pattern found else empty
            findPathOutput = strfind(SimResult.outputPathList{k}, Simulations.Output);
            if ~isempty(findPathOutput)
                % Get the simulation results
                predicted=SimResult.y{k};
                predictedTime=SimResult.time;
                
                % Convert units to reference unit
                % Caution: Code to be updated ? what is MW unit ?
                % Some concentrations are massic whereas others are molar
                Yfactor=getUnitFactor(SimResult.outputUnit{j},yAxesOptions.Unit,YDimension, 'MW',MW*1e10);
                if strcmp(ResDimension,'')
                    Resfactor=Yfactor;
                else
                    Resfactor=getUnitFactor(SimResult.outputUnit{j},yAxesOptions.Unit,ResDimension,'MW',MW*1e10);
                end
                Timefactor=getUnitFactor(SimResult.timeUnit,TimeAxesOptions.Unit,TimeDimension);
                
                predicted = predicted.*Yfactor;
                predictedRes = predicted.*Resfactor;
                predictedTime = predictedTime.*Timefactor;
                
                break
            end
        end
        % If the output was not found
        if ~exist('predicted')
            ME = MException('GOFMergedPlot:notFoundInPath', ...
                ['In Groups %d Mappings %d, ' Simulations.Output ' not found within simulation file ' csvSimFile], i, j);
            throw(ME);
        end
        
        % Get the right observations to be compared
        for k=1:length(ObservedDataSets)
            findPathOutput = strfind(ObservedDataSets(k).Id,Simulations.ObservedData);
            if ~isempty(findPathOutput)
                % Get the Observation results
                Obs=ObservedDataSets(k).y{1};
                ObsTime=ObservedDataSets(k).time;
                
                % Convert units to reference unit
                % Caution: Code to be updated
                % Some concentrations are massic whereas others are molar
                Yfactor=getUnitFactor(ObservedDataSets(k).outputUnit{1},yAxesOptions.Unit,YDimension, 'MW',MW*1e10);
                if strcmp(ResDimension,'')
                    Resfactor=Yfactor;
                else
                    Resfactor=getUnitFactor(ObservedDataSets(k).outputUnit{1},yAxesOptions.Unit,ResDimension,'MW',MW);
                end
                Timefactor=getUnitFactor(ObservedDataSets(k).timeUnit,TimeAxesOptions.Unit,TimeDimension);
                
                Obs = Obs.*Yfactor;
                ObsRes = Obs.*Resfactor;
                ObsTime = ObsTime.*Timefactor;
                
                % Get points comparable between obs and pred
                comparable_index= getComparablePredicted(ObsTime , predictedTime);
                
                % Group and save comparable points
                Group(i).dataTP(j).yobs = Obs;
                Group(i).dataTP(j).ypred = predicted(comparable_index);
                Group(i).dataTP(j).yres = predictedRes(comparable_index)-ObsRes;
                if strcmp(ResDimension,'')
                    Group(i).dataTP(j).yres=Group(i).dataTP(j).yres./ObsRes;
                end
                Group(i).dataTP(j).time = ObsTime;
                
                % Set the min/max of the axis
                minX=nanmin(min(Obs), minX);
                maxX=nanmax(max(Obs), maxX);
                minY=nanmax(min(predicted(comparable_index)), minY);
                maxY=nanmax(max(predicted(comparable_index)), maxY);
                maxtime=nanmax(max(ObsTime), maxtime);
                break
            end
        end
        % If the output was not found
        if ~exist('Obs')
            ME = MException('GOFMergedPlot:notFoundInPath', ...
                ['In Groups %d Mappings %d, ' Simulations.ObservedData ' not found within Observed Dataset'], i,j);
            throw(ME);
        end
    end
    clear Obs
    clear predicted
end

% -------------------------------------------------------------
% Plot the figures
% Observation vs Prediction Figure
figure(fig_handle1);
plot([0.8*min(minX, minY) 1.2*max(maxX, maxY)], [0.8*min(minX, minY) 1.2*max(maxX, maxY)], '--k', 'Linewidth', 1,'HandleVisibility','off');
axis([0.8*min(minX, minY) 1.2*max(maxX, maxY) 0.8*min(minX, minY) 1.2*max(maxX, maxY)]);
legendLabels={};

% Initialize error for computing GMFE
Error=[];

for i=1:length(Groups)
    
    for j=1:length(Groups(i))
        CurveOptions.Color=Groups(i).OutputMappings(j).Color;
        CurveOptions.Symbol=Groups(i).Symbol;
        CurveOptions.LineStyle='none';
        legendLabels{length(legendLabels)+1}=Groups(i).Caption;
        
        pp=plot(Group(i).dataTP(j).yobs, Group(i).dataTP(j).ypred);
        setCurveOptions(pp, CurveOptions);
        
        Error = [Error log10(Group(i).dataTP(j).ypred)-log10(Group(i).dataTP(j).yobs)];
    end
end
xLabelFinal = getLabelWithUnit('Observations',xAxesOptions.Unit);
yLabelFinal = getLabelWithUnit('Predictions',yAxesOptions.Unit);
xlabel(xLabelFinal); ylabel(yLabelFinal);

GMFE = 10.^(sum(abs(Error))/length(Error));

%xLabelFinal = getLabelWithUnit('Observed',timeUnit);
%yLabelFinal = getLabelWithUnit('Predicted',[]);
%xlabel(xLabelFinal); ylabel(yLabelFinal);
%legend(legendLabels);
legend('off')

% Residuals vs Time Figure
% create figure for Residuals vs time
figure(fig_handle2);
plot([0 1.2*maxtime], [0 0], '--k', 'Linewidth', 1, 'HandleVisibility','off');
legendLabels={};

for i=1:length(Groups)
    for j=1:length(Groups(i))
        CurveOptions.Color=Groups(i).OutputMappings(j).Color;
        CurveOptions.Symbol=Groups(i).Symbol;
        CurveOptions.LineStyle='none';
        legendLabels{length(legendLabels)+1}=Groups(i).Caption;
        
        pp=plot(Group(i).dataTP(j).time, Group(i).dataTP(j).yres);
        setCurveOptions(pp, CurveOptions);
    end
end
xLabelFinal = getLabelWithUnit('Time',TimeAxesOptions.Unit);
yLabelFinal = getLabelWithUnit('Residuals',ResAxesOptions.Unit);
xlabel(xLabelFinal); ylabel(yLabelFinal);
%legend(legendLabels);
legend('off')

% -------------------------------------------------------------
% Auxiliary function:
% get comparable time points between observations and simulations
function Index = getComparablePredicted(timeObs, timePred)

timeObs = reshape(timeObs, 1, []);
timePred = reshape(timePred, [], 1);

timeObs2 = repmat(timeObs, length(timePred), 1);
timePred2 = repmat(timePred, 1, length(timeObs));

error = (timeObs2-timePred2).*(timeObs2-timePred2);

[Y, Index] = min(error);