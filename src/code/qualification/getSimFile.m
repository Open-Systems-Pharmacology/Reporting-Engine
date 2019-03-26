function csvSimFile = getSimFile(PlotType, SimulationMapping)

csvSimFile=[];

Project=PlotType.Project;
Simulation=PlotType.Simulation;

for i=1:length(SimulationMapping)
    if strcmp(SimulationMapping(i).Project, Project) && ...
            strcmp(SimulationMapping(i).Simulation, Simulation)
        
        path = SimulationMapping(i).Path;
        csvSimFile = fullfile(path, [Simulation '-Results.csv']);
    end
end
