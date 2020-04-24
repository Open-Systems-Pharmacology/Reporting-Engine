function saveSubArtifacts(sub,DDIRatioPlots,ConfigurationPlan)

DDIRatioPlotsSub = DDIRatioPlots;
DDIRatioPlotsSub.Artifacts = DDIRatioPlots.Artifacts(ismember({'Plot','GMFE','Measure'},DDIRatioPlots.Artifacts));
 
subTypes = fieldnames(sub);
for eachSubType = 2:length(subTypes)
    curSubType = subTypes{eachSubType};
    subUnits = fieldnames(sub.(curSubType).fig_handle);
    for eachSubunit = 1:length(subUnits)
        curSubUnit = subUnits{eachSubunit};
        DDIRatioPlotsArtifacts = struct();
        
        % Save the results into a structure that can be called and
        % saved using saveArtifacts
        DDIRatioPlotsArtifacts.Plot = sub.(curSubType).fig_handle.(curSubUnit);
        DDIRatioPlotsArtifacts.Measure = sub.(curSubType).DDIRatioQuali.(curSubUnit);
        DDIRatioPlotsArtifacts.GMFE = sub.(curSubType).GMFE.(curSubUnit);

        saveArtifacts(DDIRatioPlotsArtifacts, DDIRatioPlotsSub, ConfigurationPlan, 'DDIRatio');
    end
end
