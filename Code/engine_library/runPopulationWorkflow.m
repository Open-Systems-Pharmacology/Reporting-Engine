function runPopulationWorkflow(TaskList,Settings,PopRunSet)
% master routine for population workflow
%
% Inputs:
%       TaskList (structure)    list of task which should be executed see GETDEFAULTTASKLISTPOPULATIONWORKFLOW
%       Settings (structure)    definition of properties used in all
%                   workflow functions see GETDEFAULTWORKFLOWSETTINGS
%       PopRunSet (structure)   list of population simulations see GETDEFAULTPOPRUNSET

% Open Systems Pharmacology Suite;  support@systems-biology.com
% Date: 14-July-2017


% try
    %% initialize workflow
    [Settings] = initializeWorkflow('Population',Settings);
    
    %% check if necessary inputs for task are available
    
    successInputCheck = true;

    if TaskList.simulatePopulation
        for iSet = 1:length(PopRunSet)
            successInputCheck = checkPopulationSimulationInput(Settings,PopRunSet(iSet)) & successInputCheck;
        end
    end

    if TaskList.calculatePKParameter
        if ~ TaskList.simulatePopulation
            for iSet = 1:length(PopRunSet)
                if ~exist(fullfile('Simulations',[PopRunSet(iSet).name '-Results.csv']),'file')
                    writeToLog(sprintf('ERROR: results for "%s" does not exist, please set Task simulatePopulation to true',PopRunSet(iSet).name),Settings.logfile,true,false);
                    successInputCheck = false;
                end
                successInputCheck = checkPopulationSimulationInput(Settings,PopRunSet(iSet)) & successInputCheck;
            end
        end
    end
            
    
    if ~successInputCheck
        return
    end
    
    %% Tasks processing
    
    % simulate the poulationssimulations
    if TaskList.simulatePopulation
        
        for iSet = 1:length(PopRunSet)
            runPopulationSimulation(Settings,PopRunSet(iSet));
        end
    end

    % calculate PK Parameter
    if TaskList.calculatePKParameter
        for iSet = 1:length(PopRunSet)
            calculatePopulationPKParameter(Settings,PopRunSet(iSet));
        end
    end
    
% catch exception
%     
%     save('exception.mat','exception')
%     writeToLog(exception.message,Settings.logfile,true,false);
%     
% end



return


function successInputCheck = checkPopulationSimulationInput(Settings,PopRunSet)

successInputCheck = true;

% check name on special letters
if any(ismember(PopRunSet.name,'./ ?§$%&()[]{}+~*#'))
    writeToLog(sprintf('ERROR: Popsetname "%s" contains special signs, do not use them.',PopRunSet.name),Settings.logfile,true,false);
    successInputCheck = false;
end

% mandatory inputs
if ~exist(PopRunSet.xml,'file')
    writeToLog(sprintf('ERROR: "%s" does not exist',PopRunSet.xml),Settings.logfile,true,false);
    successInputCheck = false;
end
if ~exist(PopRunSet.pop,'file')
    writeToLog(sprintf('ERROR: "%s" does not exist',PopRunSet.pop),Settings.logfile,true,false);
    successInputCheck = false;
end
if ~exist(PopRunSet.outputList,'file')
    writeToLog(sprintf('ERROR: "%s" does not exist',PopRunSet.outputList),Settings.logfile,true,false);
    successInputCheck = false;
end
                
% optional inputs
if  ~isempty(PopRunSet.studyDesign) && ~exist(PopRunSet.studyDesign,'file') 
    writeToLog(sprintf('ERROR: "%s" is given, but does not exist',PopRunSet.studyDesign),Settings.logfile,true,false);
    successInputCheck = false;
end
    
