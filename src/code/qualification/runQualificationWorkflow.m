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
for TaskListIndex=1:length(TaskList)
    % Implement Plot Settings
    if strcmp(TaskList{TaskListIndex}, 'PlotSettings')
        PlotSettings=ConfigurationPlan.Plots.PlotSettings;
        break
    else
        PlotSettings={};
    end
end

%---------------------------------------------------
for TaskListIndex=1:length(TaskList)
    
    % Implement Axes Settings
    if strcmp(TaskList{TaskListIndex}, 'AxesSettings')
        AxesSettings=ConfigurationPlan.Plots.AxesSettings;
        break
    else
        AxesSettings=[];
    end
end

%---------------------------------------------------
% Plot Time Profile
for TaskListIndex=1:length(TaskList)
    
    if strcmp(TaskList{TaskListIndex}, 'TimeProfile')
        
        for TimeProfileIndex=1:length(ConfigurationPlan.Plots.TimeProfile)
            
            TimeProfile=ConfigurationPlan.Plots.TimeProfile(TimeProfileIndex);
            
            % Update plot settings if necessary
            nPlotSettings = setPlotSettings(PlotSettings, TimeProfile);
            
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
                        Curves(kk+1).Y=TimeProfile.Plot.ObservedDataCollection.CurveOptions(kk).Path;
                        Curves(kk+1).CurveOptions=TimeProfile.Plot.ObservedDataCollection.CurveOptions(kk).CurveOptions;
                    end
                    
                    try
                        plotQualificationTimeProfile(WSettings, TimeProfileIndex, TimeProfile, ObservedDataSets, ConfigurationPlan.SimulationMappings, ...
                            Curves, PopulationAxes, nPlotSettings, ConfigurationPlan.REInput_path);
                        saveQualificationFigure(gcf, ConfigurationPlan.Sections, TimeProfile.SectionId, 'PopulationTimeProfile')
                        clear Curves PopulationAxes
                    catch exception
                        writeToReportLog('ERROR', sprintf('Error in TimeProfile plot %d "%s". \n %s \n', TimeProfileIndex, nPlotSettings.title, exception.message), 'true', exception);
                        warning('Error in TimeProfile plot %d. \n %s \n', TimeProfileIndex, exception.message);
                        close all;
                    end
                end
                
            else
                
                % Plot the Time Profile results
                try
                    plotQualificationTimeProfile(WSettings, TimeProfileIndex, TimeProfile, ObservedDataSets,ConfigurationPlan.SimulationMappings, TimeProfile.Plot.Curves, ...
                        TimeProfile.Plot.Axes, nPlotSettings, ConfigurationPlan.REInput_path);
                    % Pause option for debugging
                    % pause()
                    saveQualificationFigure(gcf, ConfigurationPlan.Sections, TimeProfile.SectionId, 'TimeProfile')
                catch exception
                    writeToReportLog('ERROR', sprintf('Error in TimeProfile plot %d "%s". \n %s \n', TimeProfileIndex, nPlotSettings.title, exception.message), 'true', exception);
                    warning('Error in TimeProfile plot %d. \n %s \n', TimeProfileIndex, exception.message);
                    close all;
                    
                end
            end
        end
        break
    end
    
end

%---------------------------------------------------
% Plot GOFMerged
for TaskListIndex=1:length(TaskList)
    if strcmp(TaskList{TaskListIndex}, 'GOFMergedPlots')
        
        for GOFMergedIndex=1:length(ConfigurationPlan.Plots.GOFMergedPlots)
            % Check if misread of GOFMergedPlots as cells
            if iscell(ConfigurationPlan.Plots.GOFMergedPlots(GOFMergedIndex))
                GOFMerged=ConfigurationPlan.Plots.GOFMergedPlots{GOFMergedIndex};
            else
                GOFMerged=ConfigurationPlan.Plots.GOFMergedPlots(GOFMergedIndex);
            end
            
            % Update plot settings if necessary
            nPlotSettings = setPlotSettings(PlotSettings, GOFMerged);
            
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
                [GOF_handle, GOFMergedGMFE] = plotQualificationGOFMerged(WSettings,GOFMergedIndex,Groups,ObservedDataSets,...
                    ConfigurationPlan.SimulationMappings, AxesOptions, nPlotSettings, ConfigurationPlan.REInput_path);
                
                % Save the results into a structure that can be called and
                % saved using saveArtifacts
                GOFMergedPlotsArtifacts(GOFMergedIndex).Plot = GOF_handle;
                GOFMergedPlotsArtifacts(GOFMergedIndex).GMFE = GOFMergedGMFE;
                
                % Get Plot Types as elements
                GOFMerged.PlotType=getElementsfromPath(GOFMerged.PlotType);
                
                saveArtifacts(GOFMergedPlotsArtifacts(GOFMergedIndex), GOFMerged, ConfigurationPlan, 'GOFMerged');
                
                clear AxesOptions nPlotSettings
                
            catch exception
                writeToReportLog('ERROR', sprintf('Error in GOFMerged plot %d "%s". \n %s \n', GOFMergedIndex, nPlotSettings.title, exception.message), 'true', exception);
                warning('Error in GOFMerged plot %d. \n %s \n', GOFMergedIndex, exception.message);
                % Close open figures
                close all
                clear AxesOptions nPlotSettings
            end
            
        end
        break
    end
end

%---------------------------------------------------
% Plot Comparison of Time profiles
for TaskListIndex=1:length(TaskList)
    if strcmp(TaskList{TaskListIndex}, 'ComparisonTimeProfilePlots')
        
        for ComparisonTimeProfileIndex=1:length(ConfigurationPlan.Plots.ComparisonTimeProfilePlots)
            ComparisonTimeProfile=ConfigurationPlan.Plots.ComparisonTimeProfilePlots(ComparisonTimeProfileIndex);
            
            % Update plot settings if necessary
            nPlotSettings = setPlotSettings(PlotSettings, ComparisonTimeProfile);
            
            for l=1:length(AxesSettings)
                if isfield(AxesSettings(l), 'ComparisonTimeProfile')
                    AxesOptions=AxesSettings(l).ComparisonTimeProfile;
                    break
                else
                    AxesOptions=[];
                end
            end
            
            % These lines are to check and convert if the OutputMapping is a cell array
            % whereas it should be a structure array
            if ~iscell(ComparisonTimeProfile.OutputMappings)
                ComparisonTimeProfile.OutputMappings = num2cell(ComparisonTimeProfile.OutputMappings);
            end
            
            
            try
                plotQualificationComparisonTimeProfile(WSettings, ComparisonTimeProfileIndex, ComparisonTimeProfile, ObservedDataSets, ConfigurationPlan.SimulationMappings, ComparisonTimeProfile.OutputMappings, ...
                    AxesOptions, nPlotSettings, ConfigurationPlan.REInput_path);
                %pause()
                saveQualificationFigure(gcf, ConfigurationPlan.Sections, ComparisonTimeProfile.SectionId, 'ComparisonTimeProfile');
                clear AxesOptions nPlotSettings
            catch exception
                writeToReportLog('ERROR', sprintf('Error in ComparisonTimeProfile plot %d "%s". \n %s \n', ComparisonTimeProfileIndex, nPlotSettings.title, exception.message), 'true', exception);
                warning('Error in ComparisonTimeProfile plot %d. \n %s \n', ComparisonTimeProfileIndex, exception.message);
                close all;
                clear AxesOptions nPlotSettings
            end
        end
        break
    end
end

%---------------------------------------------------
% Plot of PK Ratio
for TaskListIndex=1:length(TaskList)
    if strcmp(TaskList{TaskListIndex}, 'PKRatioPlots')
        
        for PKRatioPlotIndex=1:length(ConfigurationPlan.Plots.PKRatioPlots)
            % Artifacts optional field may change switch of PK Ratios
            % struct to cell format
            if iscell(ConfigurationPlan.Plots.PKRatioPlots(PKRatioPlotIndex))
                PKRatioPlots=ConfigurationPlan.Plots.PKRatioPlots{PKRatioPlotIndex};
            else
                PKRatioPlots=ConfigurationPlan.Plots.PKRatioPlots(PKRatioPlotIndex);
            end
            
            % Update plot settings if necessary
            nPlotSettings = setPlotSettings(PlotSettings, PKRatioPlots);
            
            for l=1:length(AxesSettings)
                if isfield(AxesSettings(l), 'PKRatioPlots')
                    AxesOptions=AxesSettings(l).PKRatioPlots;
                    break
                else
                    AxesOptions=[];
                end
            end
            
            % Get PK Parameters as elements
            PKRatioPlots.PKParameter=getElementsfromPath(PKRatioPlots.PKParameter);
            
            try
                % Plot the results
                [fig_handle, PKRatioTable, PKRatioQuali, PKRatioGMFE] = plotQualificationPKRatio(WSettings,PKRatioPlotIndex,PKRatioPlots.PKParameter, PKRatioPlots.Groups, ObservedDataSets, ...
                    ConfigurationPlan.SimulationMappings, AxesOptions, nPlotSettings, ConfigurationPlan.REInput_path);
                fig_handle.PKRatio.CurrentAxes.YTick=[0.1,0.25,0.5,1,2,4,10];
                
                % Save the results into a structure that can be called and
                % saved using saveArtifacts
                PKRatioPlotsArtifacts(PKRatioPlotIndex).Plot = fig_handle;
                PKRatioPlotsArtifacts(PKRatioPlotIndex).Measure = PKRatioQuali;
                PKRatioPlotsArtifacts(PKRatioPlotIndex).GMFE = PKRatioGMFE;
                PKRatioPlotsArtifacts(PKRatioPlotIndex).Table = PKRatioTable;
                
                saveArtifacts(PKRatioPlotsArtifacts(PKRatioPlotIndex), PKRatioPlots, ConfigurationPlan, 'PKRatio');
                
                clear AxesOptions nPlotSettings
                
            catch exception
                writeToReportLog('ERROR', sprintf('Error in PKRatio plot %d "%s". \n %s \n', PKRatioPlotIndex, nPlotSettings.title, exception.message), 'true', exception);
                warning('Error in PKRatio plot %d. \n %s \n', PKRatioPlotIndex, exception.message);
                % Close open figures
                close all
                clear AxesOptions nPlotSettings
            end
            
        end
        break
    end
end

%---------------------------------------------------
% Plot of DDI Ratio
for TaskListIndex=1:length(TaskList)
    if strcmp(TaskList{TaskListIndex}, 'DDIRatioPlots')
        
        for DDIRatioPlotIndex=1:length(ConfigurationPlan.Plots.DDIRatioPlots)
            % Artifacts optional field may change switch of DDI Ratios
            % struct to cell format
            if iscell(ConfigurationPlan.Plots.DDIRatioPlots(DDIRatioPlotIndex))
                DDIRatioPlots=ConfigurationPlan.Plots.DDIRatioPlots{DDIRatioPlotIndex};
            else
                DDIRatioPlots=ConfigurationPlan.Plots.DDIRatioPlots(DDIRatioPlotIndex);
            end
            
            % Update plot settings if necessary
            nPlotSettings = setPlotSettings(PlotSettings, DDIRatioPlots);
            
            for l=1:length(AxesSettings)
                if isfield(AxesSettings(l), 'DDIRatioPlotsPredictedVsObserved')
                    AxesOptions.DDIRatioPlotsPredictedVsObserved=AxesSettings(l).DDIRatioPlotsPredictedVsObserved;
                    break
                else
                    AxesOptions.DDIRatioPlotsPredictedVsObserved=[];
                end
            end
            for l=1:length(AxesSettings)
                if isfield(AxesSettings(l), 'DDIRatioPlotsResidualsVsObserved')
                    AxesOptions.DDIRatioPlotsResidualsVsObserved=AxesSettings(l).DDIRatioPlotsResidualsVsObserved;
                    break
                else
                    AxesOptions.DDIRatioPlotsResidualsVsObserved=[];
                end
            end
            
            % Get PK Parameters as elements
            DDIRatioPlots.PKParameter=getElementsfromPath(DDIRatioPlots.PKParameter);
            
            try
                % Plot the results
                [fig_handle, DDIRatioTable, DDIRatioQuali, DDIRatioGMFE] = plotQualificationDDIRatio(WSettings,DDIRatioPlotIndex,DDIRatioPlots.PKParameter, ...
                    DDIRatioPlots.Groups, ObservedDataSets, ConfigurationPlan.SimulationMappings, ...
                    AxesOptions, nPlotSettings, ConfigurationPlan.REInput_path);
                
                % Save the results into a structure that can be called and
                % saved using saveArtifacts
                DDIRatioPlotsArtifacts(DDIRatioPlotIndex).Plot = fig_handle;
                DDIRatioPlotsArtifacts(DDIRatioPlotIndex).Measure = DDIRatioQuali;
                DDIRatioPlotsArtifacts(DDIRatioPlotIndex).GMFE = DDIRatioGMFE;
                DDIRatioPlotsArtifacts(DDIRatioPlotIndex).Table = DDIRatioTable;
                
                % get Elements of Plot Type to save Artifacts
                DDIRatioPlots.PlotType=getElementsfromPath(DDIRatioPlots.PlotType);
                
                saveArtifacts(DDIRatioPlotsArtifacts(DDIRatioPlotIndex), DDIRatioPlots, ConfigurationPlan, 'DDIRatio');
                
                clear AxesOptions nPlotSettings
                
            catch exception
                %pause()
                writeToReportLog('ERROR', sprintf('Error in DDIRatio plot %d "%s". \n %s \n', DDIRatioPlotIndex, nPlotSettings.title, exception.message), 'true', exception);
                warning('Error in DDIRatio plot %d %s. \n %s \n', DDIRatioPlotIndex, nPlotSettings.title, exception.message);
                % Close open figures
                close all
                clear AxesOptions nPlotSettings
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

set(0, 'CurrentFigure', figureHandle);
set(figureHandle,'visible','off');
set(figureHandle,'PaperOrientation','portrait');

% check paper units
if strcmpi(get(figureHandle,'PaperUnits'),'centimeters')
    p = get(figureHandle,'PaperPosition');
    set(figureHandle,'PaperPosition',p*2.5)
end

% Remove Legend and Save Figure
legend('off');
saveas(figureHandle,fullfile(SectionPath, sprintf('%0.3d_plot%s', indexed_item+2, PlotType)), 'png');

% Re-add, reshape legend and save legend only
lgd = legend();
ImageLegend = Legend2Image(lgd);

imwrite(ImageLegend,sprintf('%s.png', fullfile(SectionPath, sprintf('%0.3d_legend%s', indexed_item+1, PlotType))));
close(figureHandle);

function saveQualificationTable(QualificationTable, Sections, SectionId, Type)

[SectionPath, indexed_item] = getSection(Sections, SectionId);
fileName = fullfile(SectionPath, sprintf('%0.3d_table%s.md', indexed_item+1, Type));

writeCell2md(QualificationTable, 'outfile', fileName, 'alignment', 'right');

function nPlotSettings = setPlotSettings(PlotSettings, StructureInput)

if isfield(StructureInput, 'PlotSettings')
    nPlotSettings=StructureInput.PlotSettings;
else
    nPlotSettings=PlotSettings;
end

if isfield(StructureInput, 'Title')
    nPlotSettings.title = StructureInput.Title;
elseif isfield(StructureInput, 'Caption')
    nPlotSettings.title = StructureInput.Caption;
elseif isfield(StructureInput, 'Plot')
    if isfield(StructureInput.Plot, 'Name')
        nPlotSettings.title = StructureInput.Plot.Name;
    end
else
    nPlotSettings.title = [];
end

function saveGMFE(GMFE, Sections, SectionId, PlotTitle)

[SectionPath, indexed_item] = getSection(Sections, SectionId);

% Create GMFE markdown
GMFEfile = fullfile(SectionPath, sprintf('%0.3d_%sGMFE%s', indexed_item+1, PlotTitle, '.md'));
fileID = fopen(GMFEfile,'wt');
fprintf(fileID,'GMFE = %f \n', GMFE);
fclose(fileID);

function saveArtifacts(Artifacts, PlotConfiguration, ConfigurationPlan, PlotType)

% Check Configuration of Artifacts
if ~isfield(PlotConfiguration, 'Artifacts') || isempty(PlotConfiguration.Artifacts)
    % If field missing or empty, plot everything as Plot, GMFE, Measure and Table
    PlotConfiguration.Artifacts = {'Plot', 'GMFE', 'Measure', 'Table'};
    if strcmp(PlotType, 'GOFMerged')
        % GOF Merged does not have Measure or Table output
        PlotConfiguration.Artifacts = {'Plot', 'GMFE'};
    end
end


for indexArtifacts=1:length(PlotConfiguration.Artifacts)
    % DDI and PK Ratios have Plots, GMFEs and Measures per PK parameter
    if isfield(PlotConfiguration, 'PKParameter')
        for savedPKparameters=1:length(PlotConfiguration.PKParameter)
            if strcmp(PlotConfiguration.Artifacts{indexArtifacts}, 'Plot')
                % DDI Ratio
                if isfield(PlotConfiguration, 'PlotType')
                    for plottedTypes=1:length(PlotConfiguration.PlotType)
                        if strcmp(PlotConfiguration.PlotType{plottedTypes}, 'predictedVsObserved')
                            saveQualificationFigure(Artifacts.Plot(savedPKparameters).predictedVsObserved, ConfigurationPlan.Sections, ...
                                PlotConfiguration.SectionId, sprintf('%s%sPredictedVsObserved', PlotType, PlotConfiguration.PKParameter{savedPKparameters}));
                        end
                        if strcmp(PlotConfiguration.PlotType{plottedTypes}, 'residualsVsObserved')
                            saveQualificationFigure(Artifacts.Plot(savedPKparameters).residualsVsObserved, ConfigurationPlan.Sections, ...
                                PlotConfiguration.SectionId, sprintf('%s%sResidualsVsObserved', PlotType, PlotConfiguration.PKParameter{savedPKparameters}));
                        end
                    end
                    % PK Ratio
                else
                    saveQualificationFigure(Artifacts.Plot(savedPKparameters).PKRatio, ConfigurationPlan.Sections, ...
                        PlotConfiguration.SectionId, sprintf('%s%s', PlotType, PlotConfiguration.PKParameter{savedPKparameters}));
                end
                
            end
            if strcmp(PlotConfiguration.Artifacts{indexArtifacts}, 'GMFE')
                UpdatedPlotType = sprintf('%s%s', PlotType, PlotConfiguration.PKParameter{savedPKparameters});
                saveGMFE(Artifacts.GMFE(savedPKparameters), ConfigurationPlan.Sections, PlotConfiguration.SectionId, UpdatedPlotType);
            end
            if strcmp(PlotConfiguration.Artifacts{indexArtifacts}, 'Measure')
                saveQualificationTable(Artifacts.Measure(savedPKparameters).Output, ConfigurationPlan.Sections, PlotConfiguration.SectionId, ...
                    sprintf('%s%sQualificationMeasure', PlotType, PlotConfiguration.PKParameter{savedPKparameters}));
            end
        end
        if strcmp(PlotConfiguration.Artifacts{indexArtifacts}, 'Table')
            saveQualificationTable(Artifacts.Table, ConfigurationPlan.Sections, PlotConfiguration.SectionId, sprintf('%sQualificationTable', PlotType));
        end
        % Else GOF Merged plot (Only Plot and GMFE are possible options)
    else
        if strcmp(PlotConfiguration.Artifacts{indexArtifacts}, 'Plot')
            if isfield(PlotConfiguration, 'PlotType')
                for plottedTypes=1:length(PlotConfiguration.PlotType)
                    if strcmp(PlotConfiguration.PlotType{plottedTypes}, 'predictedVsObserved')
                        saveQualificationFigure(Artifacts.Plot.PredictedVsObserved, ConfigurationPlan.Sections, ...
                            PlotConfiguration.SectionId, sprintf('%sPredictedVsObserved', PlotType));
                    end
                    if strcmp(PlotConfiguration.PlotType{plottedTypes}, 'residualsOverTime')
                        saveQualificationFigure(Artifacts.Plot.ResidualsOverTime, ConfigurationPlan.Sections, ...
                            PlotConfiguration.SectionId, sprintf('%sResidualsOverTime', PlotType));
                    end
                end
            end
        end
        if strcmp(PlotConfiguration.Artifacts{indexArtifacts}, 'GMFE')
            saveGMFE(Artifacts.GMFE, ConfigurationPlan.Sections, PlotConfiguration.SectionId, PlotType);
        end
    end
end

