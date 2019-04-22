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
            
            % Load the mapped Time Profile Simulation Results
            csvSimFile = getSimFile(TimeProfile, ConfigurationPlan.SimulationMappings);
            SimResult = loadSimResultcsv(csvSimFile, TimeProfile.Simulation);
            
            % Plot the results
            plotQualificationTimeProfile(WSettings, j, SimResult, ObservedDataSets, TimeProfile.Plot.Curves, TimeProfile.Plot.Axes, PlotSettings);
            % Pause option for debugging
            % pause()
            saveQualificationFigure(gcf, ConfigurationPlan.Sections, TimeProfile.SectionId, 'TimeProfile')
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
                if isfield(GOFMerged(j), 'PlotSettings')
                    nPlotSettings=GOFMerged(j).PlotSettings;
                else
                    nPlotSettings=PlotSettings;
                end
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
                % TO BE MODELED: Output GMFE
                plotQualificationGOFMerged(WSettings,j,Groups,ObservedDataSets,ConfigurationPlan.SimulationMappings, AxesOptions, nPlotSettings);
                
                % Pause option for debugging
                % pause()
                % Check plot type to perform predictedVsObserved, residualsOverTime or both
                if ~isempty(strfind(GOFMerged.PlotType, 'residualsOverTime'))
                    saveQualificationFigure(gcf, ConfigurationPlan.Sections, GOFMerged.SectionId, 'GOFMergedResiduals');
                end
                if ~isempty(strfind(GOFMerged.PlotType, 'predictedVsObserved'))
                    saveQualificationFigure(gcf, ConfigurationPlan.Sections, GOFMerged.SectionId, 'GOFMergedPredictedVsObserved')
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
            if isfield(AxesSettings, 'PKRatioPlots')
                AxesOptions=AxesSettings.PKRatioPlots;
            else
                AxesOptions=[];
            end
            
            % Plot the results
            plotQualificationPKRatio(WSettings,j,PKRatioPlots.PKParameter, PKRatioPlots.PKRatios,ObservedDataSets, ConfigurationPlan.SimulationMappings, AxesOptions, PlotSettings);
            
        end
    end
    break
end

%---------------------------------------------------
% Auxiliary function to save the figure in the right orientation
% using a count of the indexed files

function saveQualificationFigure(figureHandle, Sections, SectionId, PlotType)

[SectionPath, indexed_item] = getSection(Sections, SectionId);

set(figureHandle,'PaperOrientation','portrait');
saveas(figureHandle,fullfile(SectionPath, sprintf('%0.3d_plot%s', indexed_item+1, PlotType)), 'png');
close(figureHandle);
