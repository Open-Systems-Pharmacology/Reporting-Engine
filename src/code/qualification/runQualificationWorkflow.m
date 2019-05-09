function runQualificationWorkflow(WSettings, ConfigurationPlan, TaskList, ObservedDataSets)
%RUNQUALIFICATIONWORKFLOW run the Worklfow tasklist of ConfigurationPlan
%
%   runQualifiactionWorkflow(WSettings, ConfigurationPlan, TaskList, ObservedDataSets)
%       WSettings (structure)
%       ConfigurationPlan (structures array) information and simulation mappings
%       for plotting the plots in task list
%       TaskList (cells array) list of tasks/plots to be performed
%       ObservedDataSets (structures array) observations
%

% Open Systems Pharmacology Suite;  http://open-systems-pharmacology.org

%---------------------------------------------------

[WSettings] = initializeWorkflow(WSettings);

%---------------------------------------------------
for i=1:length(TaskList)
    % Implement Plot Settings
    if strcmp(TaskList{i}, 'PlotSettings')
        PlotSettings=ConfigurationPlan.Plots.PlotSettings;
        break
    else
        PlotSettings={};
    end
end

%---------------------------------------------------
for i=1:length(TaskList)
    
    % Implement Axes Settings
    if strcmp(TaskList{i}, 'AxesSettings')
        AxesSettings=ConfigurationPlan.Plots.AxesSettings;
        break
    else
        AxesSettings=[];
    end
end


%---------------------------------------------------

% Plot Time Profile
for i=1:length(TaskList)
    
    if strcmp(TaskList{i}, 'TimeProfile')
        
        for j=1:length(ConfigurationPlan.Plots.TimeProfile)
            
            TimeProfile=ConfigurationPlan.Plots.TimeProfile(j);
            
            % Update plot settings if necessary (to be performed)
            % Currently not handled
            if isfield(TimeProfile, 'PlotSettings')
                nPlotSettings=TimeProfile.PlotSettings;
            else
                nPlotSettings=PlotSettings;
            end
            nPlotSettings.title = TimeProfile.Plot.Name;
            
            % Check if Individual or Population Time Profile
            if isfield(TimeProfile.Plot, 'Type')
                % If loop subject to change:
                % So far, Type is TimeProfile for population, and no field
                % Type is entered for individual.
                if ~iscell(TimeProfile.Plot.Analysis.Fields)
                    TimeProfile.Plot.Analysis.Fields=num2cell(TimeProfile.Plot.Analysis.Fields);
                end
                if ~iscell(TimeProfile.Plot.ObservedDataCollection.ObservedData)
                    TimeProfile.Plot.ObservedDataCollection.ObservedData=num2cell(TimeProfile.Plot.ObservedDataCollection.ObservedData);
                end
                
                % Get analysis plan passed on for each field
                for jj=1:length(TimeProfile.Plot.Analysis.Fields)
                    % Curves and Axis may contain all the analyis in a specific
                    % field. But the path will be called similarly.
                    PopulationCurves = TimeProfile.Plot.Analysis.Fields{jj};
                    
                    PopulationAxes(1).Type='X';
                    PopulationAxes(1).Dimension='Time';
                    PopulationAxes(1).Unit='h';
                    PopulationAxes(1).Scaling='Linear';
                    
                    PopulationAxes(2).Type='Y';
                    PopulationAxes(2).Dimension=PopulationCurves.Dimension;
                    PopulationAxes(2).Unit=PopulationCurves.Unit;
                    PopulationAxes(2).Scaling=PopulationCurves.Scaling;
                    
                    Curves(1).Name=PopulationCurves.Name;
                    Curves(1).X='Time';
                    Curves(1).Y=PopulationCurves.QuantityPath;
                    Curves(1).Type='Population';
                    Curves(1).Statistics = TimeProfile.Plot.Analysis.Statistics;
                    
                    % Add observed data into curves structures
                    for kk=1:length(TimeProfile.Plot.ObservedDataCollection.CurveOptions)
                        Curves(kk+1).Name=TimeProfile.Plot.ObservedDataCollection.CurveOptions(kk).Caption;
                        Curves(kk+1).X='Time';
                        Curves(kk+1).Y=TimeProfile.Plot.ObservedDataCollection.ObservedData{kk};
                        Curves(kk+1).CurveOptions=TimeProfile.Plot.ObservedDataCollection.CurveOptions(kk).CurveOptions;
                    end
                    
                    try
                        plotQualificationTimeProfile(WSettings, jj, TimeProfile, ObservedDataSets, ConfigurationPlan.SimulationMappings, Curves, PopulationAxes, nPlotSettings);
                        saveQualificationFigure(gcf, ConfigurationPlan.Sections, TimeProfile.SectionId, 'PopulationTimeProfile')
                        clear Curves PopulationAxes
                    catch exception
                        writeToReportLog('ERROR', sprintf('Error in TimeProfile plot %d. \n %s \n', j, exception.message), 'true', exception);
                        warning('Error in TimeProfile plot %d. \n %s \n', j, exception.message);
                        close all;
                    end
                end
                
            else
                
                % Plot the Time Profile results
                try
                    plotQualificationTimeProfile(WSettings, j, TimeProfile, ObservedDataSets,ConfigurationPlan.SimulationMappings, TimeProfile.Plot.Curves, TimeProfile.Plot.Axes, nPlotSettings);
                    % Pause option for debugging
                    % pause()
                    saveQualificationFigure(gcf, ConfigurationPlan.Sections, TimeProfile.SectionId, 'TimeProfile')
                catch exception
                    writeToReportLog('ERROR', sprintf('Error in TimeProfile plot %d. \n %s \n', j, exception.message), 'true', exception);
                    warning('Error in TimeProfile plot %d. \n %s \n', j, exception.message);
                    close all;
                    
                end
            end
        end
        break
    end
    
end

%---------------------------------------------------
% Plot GOFMerged
for i=1:length(TaskList)
    if strcmp(TaskList{i}, 'GOFMergedPlots')
        
        for j=1:length(ConfigurationPlan.Plots.GOFMergedPlots)
            
            GOFMerged=ConfigurationPlan.Plots.GOFMergedPlots(j);
            
            for k=1:length(GOFMerged)
                
                % Update plot settings if necessary (to be performed)
                % Currently not handled
                if isfield(GOFMerged(k), 'PlotSettings')
                    nPlotSettings=GOFMerged(k).PlotSettings;
                else
                    nPlotSettings=PlotSettings;
                end
                nPlotSettings.title = GOFMerged(k).Caption;
                for l=1:length(AxesSettings)
                    if isfield(AxesSettings(l), 'GOFMergedPlotsPredictedVsObserved')
                        AxesOptions.GOFMergedPlotsPredictedVsObserved=AxesSettings(l).GOFMergedPlotsPredictedVsObserved;
                        break
                    else
                        AxesOptions.GOFMergedPlotsPredictedVsObserved=[];
                    end
                end
                for l=1:length(AxesSettings)
                    if isfield(AxesSettings(l), 'GOFMergedPlotsResidualsOverTime')
                        AxesOptions.GOFMergedPlotsResidualsOverTime=AxesSettings(l).GOFMergedPlotsResidualsOverTime;
                        break
                    else
                        AxesOptions.GOFMergedPlotsResidualsOverTime=[];
                    end
                end
                % Get all the groups for one GOF merged plot
                Groups = GOFMerged.Groups;
                
                % Plot the Goodness of fit as obs vs pred and residuals
                try
                    GMFE = plotQualificationGOFMerged(WSettings,j,Groups,ObservedDataSets,ConfigurationPlan.SimulationMappings, AxesOptions, nPlotSettings);
                    
                    % Pause option for debugging
                    % pause()
                    % Check plot type to perform predictedVsObserved, residualsOverTime or both
                    if ~isempty(strfind(GOFMerged.PlotType, 'residualsOverTime'))
                        saveQualificationFigure(gcf, ConfigurationPlan.Sections, GOFMerged.SectionId, 'GOFMergedResiduals');
                    end
                    if ~isempty(strfind(GOFMerged.PlotType, 'predictedVsObserved'))
                        saveQualificationFigure(gcf, ConfigurationPlan.Sections, GOFMerged.SectionId, 'GOFMergedPredictedVsObserved');
                    end
                    [SectionPath, indexed_item] = getSection(ConfigurationPlan.Sections, GOFMerged.SectionId);
                    % Create GMFE markdown
                    GMFEfile = fullfile(SectionPath, sprintf('%0.3d_GMFE%s', indexed_item+1, '.md'));
                    fileID = fopen(GMFEfile,'wt');
                    fprintf(fileID,'GMFE = %f \n',GMFE);
                    fclose(fileID);
                catch exception
                    writeToReportLog('ERROR', sprintf('Error in GOFMerged plot %d, Group %d. \n %s \n', j, k, exception.message), 'true', exception);
                    warning('Error in GOFMerged plot %d, Group %d. \n %s \n', j, k, exception.message);
                    % Close open figures
                    close all
                end
                
            end
            
        end
        break
    end
end


%---------------------------------------------------
% Plot Comparison of Time profiles
% TO BE COMPLETED
for i=1:length(TaskList)
    if strcmp(TaskList{i}, 'ComparisonTimeProfilePlots')
        %{
        for j=1:length(ConfigurationPlan.Plots.ComparisonTimeProfilePlots)
            
            ComparisonTimeProfile=ConfigurationPlan.Plots.ComparisonTimeProfilePlots(j);
            
            %            for k=1:length(ComparisonTimeProfile.OutputMappings)
            %                csvSimFile = getSimFile(TimeProfile, ComparisonTimeProfile.OutputMappings(k));
            %                SimResult(k) = loadSimResultcsv(csvSimFile, ComparisonTimeProfile.OutputMappings(k).Simulation);
            %            end
            
            % Update plot settings if necessary (to be performed)
            % Currently not handled
            if isfield(ComparisonTimeProfile, 'PlotSettings')
                nPlotSettings=ComparisonTimeProfile.PlotSettings;
            else
                nPlotSettings=PlotSettings;
            end
            
            % Load the mapped Time Profile Simulation Results
            csvSimFile = getSimFile(TimeProfile, ConfigurationPlan.SimulationMappings);
            SimResult = loadSimResultcsv(csvSimFile, TimeProfile.Simulation);
        
        
            ax = getReportFigure(WSettings,1,1,[],'figureformat','landscape');
            
            % Plot the results
            plotQualificationComparisonTimeProfile(WSettings, j, SimResult, ComparisonTimeProfile.OutputMapping(k).Plot.Curves, ComparisonTimeProfile.OutputMapping(k).Plot.Axes);
            
            % Pause option for debugging
            % pause()
            saveQualificationFigure(gcf, ConfigurationPlan.Sections, TimeProfile.SectionId, 'ComparisonTimeProfile')
        end
        %}
        break
    end
end

%---------------------------------------------------
% Plot of PK Ratio
for i=1:length(TaskList)
    if strcmp(TaskList{i}, 'PKRatioPlots')
        
        for j=1:length(ConfigurationPlan.Plots.PKRatioPlots)
            
            PKRatioPlots=ConfigurationPlan.Plots.PKRatioPlots(j);
            
            % Update plot settings if necessary (to be performed)
            % Currently not handled
            if isfield(PKRatioPlots, 'PlotSettings')
                nPlotSettings=PKRatioPlots.PlotSettings;
            else
                nPlotSettings=PlotSettings;
            end
            nPlotSettings.title = PKRatioPlots.Caption;
            if isfield(PKRatioPlots, 'Color')
                CurveOptions.Color=PKRatioPlots.Color;
            else
                CurveOptions.Color='#000000'; %Black
            end
            if isfield(PKRatioPlots, 'Symbol')
                CurveOptions.Symbol=PKRatioPlots.Symbol;
            else
                CurveOptions.Symbol='Circle';
            end
            if isfield(AxesSettings, 'PKRatioPlots')
                AxesOptions=AxesSettings.PKRatioPlots;
            else
                AxesOptions=[];
            end
            
            try
                % Plot the results
                [PKRatioTable, GMFE] = plotQualificationPKRatio(WSettings,j,PKRatioPlots.PKParameter, PKRatioPlots.PKRatios,ObservedDataSets, ConfigurationPlan.SimulationMappings, AxesOptions, nPlotSettings, CurveOptions);
                
                saveQualificationFigure(gcf, ConfigurationPlan.Sections, PKRatioPlots.SectionId, 'PKRatio');
                saveQualificationTable(PKRatioTable, ConfigurationPlan.Sections, PKRatioPlots.SectionId, 'PKRatio');
                
                [SectionPath, indexed_item] = getSection(ConfigurationPlan.Sections, PKRatioPlots.SectionId);
                % Create GMFE markdown
                GMFEfile = fullfile(SectionPath, sprintf('%0.3d_GMFE%s', indexed_item+1, '.md'));
                fileID = fopen(GMFEfile,'wt');
                fprintf(fileID,'GMFE = %f \n',GMFE);
                fclose(fileID);
            catch exception
                writeToReportLog('ERROR', sprintf('Error in PKRatio plot %d. \n %s \n', j, exception.message), 'true', exception);
                warning('Error in PKRatio plot %d. \n %s \n', j, exception.message);
                % Close open figures
                close all
            end
        end
        break
    end
end

%---------------------------------------------------
% Auxiliary function to save the figure in the right orientation
% using a count of the indexed files

function saveQualificationFigure(figureHandle, Sections, SectionId, PlotType)

[SectionPath, indexed_item] = getSection(Sections, SectionId);

set(figureHandle,'PaperOrientation','portrait');
saveas(figureHandle,fullfile(SectionPath, sprintf('%0.3d_plot%s', indexed_item+1, PlotType)), 'png');
close(figureHandle);

function saveQualificationTable(QualificationTable, Sections, SectionId, Type)

[SectionPath, indexed_item] = getSection(Sections, SectionId);
fileName = fullfile(SectionPath, sprintf('%0.3d_table%s.md', indexed_item+1, Type));

writeCell2md(QualificationTable, 'outfile', fileName, 'alignment', 'right');


