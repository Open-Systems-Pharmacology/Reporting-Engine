function [ApplicationProtocol,isValid] = getApplicationProtocollFromXML(Settings,xml)
% GETAPPLICATIONPROTOCOLLFROMXML get properties of the applicationsprotocol
%       defined within an exported simulationfile (.xml)
%
% Inputs:
%       Settings (structure)    definition of properties used in all
%                   workflow functions see GETDEFAULTWORKFLOWSETTINGS
%   xml (string) name of xmlfile
% 
% Outputs;
%   ApplicationProtocol structure with following fields
%       name                (string)    name of application
%       startTime           (double)    [min] starting time of application
%       dose                (double)    [kg]  applicated drug amount
%       dosePerBodyWeight   (double)    [kg/kg] applicated drug amount per body weight
%       drugMass            (double)    [µmol] applicated drug amount
%       isDosePerBodyweight (boolean)   if true, dose is adjusted to bodyweight
%                                       false dose is absolute  
%       infusionTime        (double)    [min] time of infusion, 
%                                       zero if the application is no infusion
%  isValid  (boolean)  if true,  all necessary field were available in the xml file, 
%                      if false, some fields could not be found, may be the xml is a 
%                       MoBi file with a user defined application. Special attention is needed 


% Open Systems Pharmacology Suite;  support@systems-biology.com
% Date: 14-July-2017

initSimulation(xml,'none');
simulationIndex=1;

% select application parameter 
% initialize loop parameters

iApp = 1;
ApplicationProtocol = getDefaultApplicationProtocol('');
isValid = true;

while  existsParameter(sprintf('*Application_%d|ProtocolSchemaItem*',iApp),simulationIndex,'parametertype','readonly')

    % initialize
    ApplicationProtocol(iApp) = getDefaultApplicationProtocol(sprintf('Application_%d',iApp));
    
    % get the relevant Parameter
    if existsParameter(sprintf('*Application_%d|ProtocolSchemaItem|Start time',iApp),simulationIndex,'parametertype','readonly');
        ApplicationProtocol(iApp).startTime = getParameter(sprintf('*Application_%d|ProtocolSchemaItem|Start time',iApp),simulationIndex,'parametertype','readonly');
    end
    if existsParameter(sprintf('*Application_%d|ProtocolSchemaItem|Dose',iApp),simulationIndex,'parametertype','readonly');
        ApplicationProtocol(iApp).dose = getParameter(sprintf('*Application_%d|ProtocolSchemaItem|Dose',iApp),simulationIndex,'parametertype','readonly');
        ApplicationProtocol(iApp).dose =  ApplicationProtocol(iApp).dose; 
        ApplicationProtocol(iApp).isDosePerBodyweight = getParameter(sprintf('*Application_%d|ProtocolSchemaItem|Dose',iApp),simulationIndex,'parametertype','readonly','property','isFormula');
    end
    if existsParameter(sprintf('*Application_%d|ProtocolSchemaItem|DrugMass',iApp),simulationIndex,'parametertype','readonly');
        ApplicationProtocol(iApp).drugMass = getParameter(sprintf('*Application_%d|ProtocolSchemaItem|DrugMass',iApp),simulationIndex,'parametertype','readonly');
    end
    if existsParameter(sprintf('*Application_%d|ProtocolSchemaItem|DosePerBodyWeight',iApp),simulationIndex,'parametertype','readonly');
        ApplicationProtocol(iApp).dosePerBodyWeight = getParameter(sprintf('*Application_%d|ProtocolSchemaItem|DosePerBodyWeight',iApp),simulationIndex,'parametertype','readonly');
    end
    if existsParameter(sprintf('*Application_%d|ProtocolSchemaItem|Infusion time',iApp),simulationIndex,'parametertype','readonly');
        ApplicationProtocol(iApp).infusionTime =   getParameter(sprintf('*Application_%d|ProtocolSchemaItem|Infusion time',iApp),simulationIndex,'parametertype','readonly');
    end
        
       
    iApp = iApp+1;
end

% take only application within simulation range
simulationTime = getSimulationTime(simulationIndex);
jj = [ApplicationProtocol.startTime] < simulationTime(end);
ApplicationProtocol = ApplicationProtocol(jj);

% take only application with a Drugmass > 0
jj = [ApplicationProtocol.drugMass] > 0;
ApplicationProtocol = ApplicationProtocol(jj);

% do consistency checks
for iApp = 1:length(ApplicationProtocol)
    tmp = struct2cell(ApplicationProtocol(iApp));
    jj = cellfun(@isnumeric,tmp);
    if any(isnan(cell2mat(tmp(jj))))
        isValid = false;
    
        writeToLog('Applicationprotocoll is a non default protocol, which is not possible  to interpret',...
            Settings.logfile,true,false);
        return
    end
end


return


function ApplicationProtocol = getDefaultApplicationProtocol(name)

 ApplicationProtocol =    struct('name',name,'startTime',nan,...
        'dose',nan,'dosePerBodyWeight',nan,'drugMass',nan,'isDosePerBodyweight',false,...
        'infusionTime',0);
