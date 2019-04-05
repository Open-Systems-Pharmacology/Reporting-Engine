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

for i=1:length(TaskList)
    
    % Plot Time Profile
    if strcmp(TaskList{i}, 'TimeProfile')
        
        for j=1:2%length(ConfigurationPlan.Plots.TimeProfile)
            
            TimeProfile=ConfigurationPlan.Plots.TimeProfile(j);
            
            % Update plot settings if necessary (to be performed)
            % WSettings = getDefaultWorkflowSettings('',''); close;
            
            % Load the mapped Time Profile Simulation Results
            csvSimFile = getSimFile(TimeProfile, ConfigurationPlan.SimulationMappings);
            SimResult = loadSimResultcsv(csvSimFile, TimeProfile.Simulation);
            
            % Plot the results
            plotQualificationTimeProfile(WSettings, j, SimResult, ObservedDataSets, TimeProfile.Plot.Curves, TimeProfile.Plot.Axes);
            
            SectionPath = getSection(ConfigurationPlan.Sections, TimeProfile.SectionId);
            saveas(gcf,fullfile(SectionPath, TimeProfile.Plot.Name), 'fig');
            close;
        end
    end
    
    
    % Plot GOFMerged
    if strcmp(TaskList{i}, 'GOFMerged')
        
    end
end