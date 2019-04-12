function  GMFE=plotQualificationGOFMerged(WSettings,figureHandle,Groups,ObservedDataSets,SimulationMappings, PlotSettings)
% PLOTQUALIFICATIONGOFMERGED plot predcited vs observed and residuals vs
% time
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
%   PlotSettings (structure) to set plot options
%
% Output
%   GMFE (numeric) measure of goodness of fit

% Open Systems Pharmacology Suite;  http://open-systems-pharmacology.org

% Initialize the unity line and the time line of the plot
minX=NaN; maxX=NaN; minY=NaN; maxY=NaN; maxtime=NaN;

% Map for each group the corresponding observations within a same structure
for i=1:length(Groups)
    % Load simulation according to mapping
    for j=1:length(Groups(i))
        Simulations = Groups(i).OutputMappings(j);
        csvSimFile = getSimFile(Simulations, SimulationMappings);
        SimResult = loadSimResultcsv(csvSimFile, Simulations);
        
        % Get the right simulation output to be compared
        for k=1:length(SimResult.outputPathList)
            % Find first integer of string pattern found else empty
            findPathOutput = strfind(SimResult.outputPathList{k}, Simulations.Output);
            if ~isempty(findPathOutput)
                predicted=SimResult.y{k};
                predictedTime=SimResult.time;
                predictedTimeUnit=SimResult.timeUnit;
                break
            end
        end
        % If the output was not found
        if ~exist('predicted')
            ME = MException('GOFMergedPlot:notFoundInPath', ...
                ['In Groups %d Mappings %d, ' Simulations.Output ' not found within simulation file ' csvSimFile], i, j);
            throw(ME);
        end
        for k=1:length(ObservedDataSets)
            findPathOutput = strfind(ObservedDataSets(k).Id,Simulations.ObservedData);
            if ~isempty(findPathOutput)
                Group(i).dataTP(j).y=ObservedDataSets(k).y{1};
                Group(i).dataTP(j).time=ObservedDataSets(k).time;
                
                % Get points comparable between obs and pred
                timeUnit = ObservedDataSets(k).timeUnit;
                predictedTime = ConvertTimeUnit(predictedTime, predictedTimeUnit, timeUnit);
                comparable_index= getComparablePredicted(Group(i).dataTP(j).time, predictedTime);
                Group(i).dataTP(j).predicted=predicted(comparable_index);
                
                minX=nanmin(min(Group(i).dataTP(j).y), minX);
                maxX=nanmax(max(Group(i).dataTP(j).predicted), maxX);
                minY=nanmax(min(Group(i).dataTP(j).y), minY);
                maxY=nanmax(max(Group(i).dataTP(j).predicted), maxY);
                
                maxtime=nanmax(max(Group(i).dataTP(j).time), maxtime);
                break
            end
        end
    end
end

% create figure for Obs vs Pred
ax = getReportFigure(WSettings,1,1,figureHandle,'figureformat','square');
plot([0.8*min(minX, minY) 1.2*max(maxX, maxY)], [0.8*min(minX, minY) 1.2*max(maxX, maxY)], '--k', 'Linewidth', 1,'HandleVisibility','off');
legendLabels={};

for i=1:length(Groups)
    for j=1:length(Groups(i))
        CurveOptions.Color=Groups(i).OutputMappings(j).Color;
        CurveOptions.Symbol=Groups(i).Symbol;
        CurveOptions.LineStyle='none';
        legendLabels{length(legendLabels)+1}=Groups(i).Caption;
        
        pp=plot(Group(i).dataTP(j).y, Group(i).dataTP(j).predicted);
        setCurveOptions(pp, CurveOptions);
    end
end
xlabel('Observed'); ylabel('Predicted');
legend(legendLabels);

% create figure for Residuals vs time
ax = getReportFigure(WSettings,1,1,figureHandle+1,'figureformat','square');
plot([0 1.2*maxtime], [0 0], '--k', 'Linewidth', 1, 'HandleVisibility','off');
legendLabels={};

for i=1:length(Groups)
    for j=1:length(Groups(i))
        CurveOptions.Color=Groups(i).OutputMappings(j).Color;
        CurveOptions.Symbol=Groups(i).Symbol;
        CurveOptions.LineStyle='none';
        legendLabels{length(legendLabels)+1}=Groups(i).Caption;
        
        pp=plot(Group(i).dataTP(j).time, Group(i).dataTP(j).y-Group(i).dataTP(j).predicted);
        setCurveOptions(pp, CurveOptions);
    end
end
xlabel(['Time ' timeUnit]); ylabel('Residuals');
legend(legendLabels);

function Index = getComparablePredicted(timeObs, timePred)

timeObs = reshape(timeObs, 1, []);
timePred = reshape(timePred, [], 1);

timeObs2 = repmat(timeObs, length(timePred), 1);
timePred2 = repmat(timePred, 1, length(timeObs));

error = (timeObs2-timePred2).*(timeObs2-timePred2);

[Y, Index] = min(error);