function csvSimFile = getSimFile(PlotType, SimulationMapping)
% GETSIMFILE get the simulated result file specified in the
% SimulationMapping
%
% csvSimFile = getSimFile(PlotType, SimulationMapping)
%
%
%       csvSimFile (structure) contains all relevant information of
%       Simulation results csv file
%       PlotType (structure) Project and Simulation to be mapped
%       SimulationMapping (structure) contains all the mapping of file to
%       be loaded in function of Project and Simulation
%

% Open Systems Pharmacology Suite;  http://open-systems-pharmacology.org

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
