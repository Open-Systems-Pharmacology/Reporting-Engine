function runQualificationWorkflow(WSettings, ConfigurationPlan, TaskList, ObservedDataSets)
%RUNQUALIFICATIONWORKFLOW run the TaskList from the ConfigurationPlan
%
%   runQualifiactionWorkflow(ConfigurationPlan,TaskList)
%       ConfigurationPlan (array of structures) information of the
%       Qualification Plan
%       TaskList (array of cells) list of plots to be performed
%

% Open Systems Pharmacology Suite;  http://open-systems-pharmacology.org

%---------------------------------------------------

% Implement Plot Settings
if isfield(ConfigurationPlan.Plots, 'PlotSettings')
    PlotSettings=ConfigurationPlan.Plots.PlotSettings;
else
    PlotSettings={};
end


for i=1:length(TaskList)
    
    % Plot Time Profile
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
            SectionPath = getSection(ConfigurationPlan.Sections, TimeProfile.SectionId);
            saveas(gcf,fullfile(SectionPath, ['OutputTimeProfile' num2str(j)]), 'pdf');
            close;
        end
    end
    
    
    % Plot GOFMerged
    % Caution, function not finished, to be completed
    % So far save empty plots
    if strcmp(TaskList{i}, 'GOFMergedPlots')
        
        for j=1:length(ConfigurationPlan.Plots.GOFMergedPlots)
            
            GOFMerged=ConfigurationPlan.Plots.GOFMergedPlots(j);
            
            % Update plot settings if necessary (to be performed)
            % Currently not handled
            WSettings = getDefaultWorkflowSettings('Qualification','default'); close;
            
            for k=1:length(GOFMerged)
                
                % Update plot settings if necessary (to be performed)
                % Currently not handled
                if isfield(GOFMerged(j), 'PlotSettings')
                    nPlotSettings=GOFMerged(j).PlotSettings;
                else
                    nPlotSettings=PlotSettings;
                end
                
                % Get all the groups for one GOF merged plot
                Groups = GOFMerged.Groups;
                
                % Plot the Goodness of fit and output GMFE
                GOFResults=plotQualificationGOFMerged(WSettings,j,Groups,ObservedDataSets,ConfigurationPlan.SimulationMappings, nPlotSettings);
                
                % Pause option for debugging
                % pause()
                
                % Check plot type to perform predictedVsObserved, residualsOverTime or both
                if ~isempty(strfind(GOFMerged.PlotType, 'residualsOverTime'))
                    SectionPath = getSection(ConfigurationPlan.Sections, GOFMerged.SectionId);
                    saveas(gcf,fullfile(SectionPath, ['OutputGOFMergedResiduals' num2str(k)]), 'png');
                    close(gcf);
                end
                if ~isempty(strfind(GOFMerged.PlotType, 'predictedVsObserved'))
                    SectionPath = getSection(ConfigurationPlan.Sections, GOFMerged.SectionId);
                    saveas(gcf,fullfile(SectionPath, ['OutputGOFMergedPredvsObs' num2str(k)]), 'png');
                    close(gcf);
                end
            end
            
        end
    end
    
    % Plot Comparison of Time profiles
    % TO BE IMPLEMENTED
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
            WSettings = getDefaultWorkflowSettings('Qualification','default'); close;
            
            ax = getReportFigure(WSettings,1,1,[],'figureformat','landscape');
            
            %plotReportPredictedVsObserved(WSettings,[],DataTP,yLabel,yUnit,legendEntries,yscale,lloq);
            
            % Plot the results
            %plotQualificationComparisonTimeProfile(WSettings, j, SimResult, ComparisonTimeProfile.OutputMapping(k).Plot.Curves, ComparisonTimeProfile.OutputMapping(k).Plot.Axes);
            
            SectionPath = getSection(ConfigurationPlan.Sections, ComparisonTimeProfile.SectionId);
            saveas(gcf,fullfile(SectionPath, ['ComparisonTimeProfilePlots' num2str(j)]), 'png');
            close;
        end
        %}
    end
    
end