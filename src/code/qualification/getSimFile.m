function [csvSimFile, xmlSimFile] = getSimFile(PlotType, SimulationMapping, REInputPath)
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
xmlSimFile=[];

Project=PlotType.Project;
Simulation=PlotType.Simulation;

for i=1:length(SimulationMapping)
    if strcmp(SimulationMapping(i).Project, Project) && ...
            strcmp(SimulationMapping(i).Simulation, Simulation)
        
        path = SimulationMapping(i).Path;
        
        % The name of the Simulation file is the same as folder name
        % and may mismatch from the Simulation Id
        % Consequently, the name of the folder is used instead
        SimulationFile = getElementsfromPath(fullfile(path), '\');
        
        csvSimFile = fullfile(REInputPath, path, [SimulationFile{2} '-Results.csv']);
        xmlSimFile = fullfile(REInputPath, path, [SimulationFile{2} '.xml']);
        break
    end
end