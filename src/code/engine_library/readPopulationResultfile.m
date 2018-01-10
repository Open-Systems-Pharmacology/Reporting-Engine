function nBunch = readPopulationResultfile(simulationName,simulationType,simResultPrefix)
%READPOPULATIONRESULTFILE reads a simulation result file exported by PK-SIM
% and converts it to the temporary matfile
%
% Inputs: 
%   - simulationName (string) name of the simulation 
%   - simulationType (string) "Results" for population, 
%                             "MeanModel" for a MeanModel simulation
%   - simResultPrefix (string) prefix for matfile
% 
% Outputs  nBunch (double) number of files, simulation is result is splitted into different bunches
%
%   temporary result is a structure
%   - SimResult structure with following fields
%       name (string)   name of simulation
%       time (doublevector)     timevector
%       values (cellarray):  the ith entry contains a double matrix with the
%           concentration values of the ith quantity ( modeloutput)
%       individualIdVector (double vector):  vector with the individual ids
%       outputPathList (cellarray):  the ith entry contains a string with the
%          pathname of the ith quantity ( modeloutput)
%       outputUnit (cellarray):  the ith entry contains a string with the
%          unit of the ith quantity ( modeloutput)

% Open Systems Pharmacology Suite;  http://open-systems-pharmacology.org


% getName of the csvfile
csvfile = fullfile('simulations',sprintf('%s-%s.csv',simulationName,simulationType)); 

iBunch = 1;
while exist(csvfile,'file')
   
    % read data
    t = readtable(csvfile);
    % outputs
    varNames = strrep(strrep(t.Properties.VariableDescriptions,'Original column heading: ',''),'''','');
    jj = cellfun(@isempty,varNames);
    varNames(jj) = t.Properties.VariableNames(jj);
    outputPathList = varNames(3:end);
        
    
    % split path and unit
    for iP=1:length(outputPathList)
        tmp=outputPathList{iP};
        ji=strfind(tmp,'[');
        if ~isempty(ji)
            outputPathList{iP}=strtrim(tmp(1:ji-1));
            outputUnit{iP}=strtrim(tmp(ji:end)); %#ok<AGROW>
        else
            outputPathList{iP}=tmp;
            outputUnit{iP}=''; %#ok<AGROW>
        end
    end
    
    % read numeric data
    individualIdVector = unique(t.IndividualId);
    time = unique(t{:,2});
    
    % get values for each output
    for iO=3:size(t,2)
        values{iO-2}=reshape(t{:,iO},length(time),length(individualIdVector)); %#ok<AGROW>
    end
    
    % collect infos in Structure
    SimResult.name = simulationName;
    SimResult.time = time';
    SimResult.values = values';
    SimResult.individualIdVector = individualIdVector;
    SimResult.outputPathList = outputPathList;
    SimResult.outputUnit = outputUnit;
    
    % save as temporary file
    save(fullfile('tmp',simulationName,sprintf('%ssimResult_%d.mat',simResultPrefix,iBunch)),'SimResult');
    
    
    % get name of resultfile
    iBunch = 1 + iBunch;
    csvfile = fullfile('simulations',sprintf('%s-%d-%s.csv',simulationName,iBunch,simResultPrefix));

end

nBunch = iBunch-1;

return
