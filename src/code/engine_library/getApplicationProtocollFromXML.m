function [ApplicationProtocol,isValid] = getApplicationProtocollFromXML(WSettings,xmlOrIndex,parPathsStudyDesign,parValuesStudyDesign)
% GETAPPLICATIONPROTOCOLLFROMXML get properties of the applicationsprotocol
%       defined within an exported simulationfile (.xml)
%
% [ApplicationProtocol,isValid] = getApplicationProtocollFromXML(WSettings,xml,parPathsStudyDesign,parValuesStudyDesign)
% 
% Inputs:
%       WSettings (structure)    definition of properties used in all
%                   workflow functions see GETDEFAULTWORKFLOWSETTINGS
%   xmlOrIndex (string/double) name of xmlfile or if already initialized simulatin Index
%   parPathsStudyDesign (cellarray of strings) optional, list of studedsign parameters
%   parValuesStudyDesign (double matrix) optional, values of studedsign parameters
% 
% Outputs:
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
%       compound        (string)        name of applicated compound
%  isValid  (boolean)  if true,  all necessary field were available in the xml file, 
%                      if false, some fields could not be found, may be the xml is a 
%                       MoBi file with a user defined application. Special attention is needed 


% Open Systems Pharmacology Suite;  http://open-systems-pharmacology.org

% initialize return parameter
isValid = true;
ApplicationProtocol = getDefaultApplicationProtocol('');




% initialize simulation file
if isnumeric(xmlOrIndex)
    simulationIndex = xmlOrIndex;
else
    initSimulation(xmlOrIndex,'none');
    simulationIndex=1;
end

% get list of paths added by studydesign
if exist('parPathsStudyDesign','var') && ~isempty(parPathsStudyDesign)
    studyDesignPathId = cell(length(parPathsStudyDesign),1);
    for iPar = 1:length(parPathsStudyDesign)
        studyDesignPathId{iPar} = getParameter(['*' parPathsStudyDesign{iPar}],simulationIndex,'parametertype','readonly','property','Path');
    end
else
    studyDesignPathId = {};
    parValuesStudyDesign = [];
end

% select application parameter 
[ise,desc] = existsParameter(sprintf('*Application_*|ProtocolSchemaItem|Start time'),simulationIndex,'parametertype','readonly');
if ise
    isIndividualized = false;
    for iApp = 1:size(desc,1)-1
        
        tmp = regexp(desc{iApp+1,2},'\|','split');
        pathPrefix = strjoin(tmp(1:end-1),'|');
        pathCompound = [strjoin(tmp(1:end-2),'|') '|*|Concentration'];
        
        % initialize
        ApplicationProtocol(iApp) = getDefaultApplicationProtocol( strjoin(tmp(3:(end-2)),'|')); 
        
        % get the relevant Parameter
        [ApplicationProtocol(iApp),isValid,isIndividualized] = addParameter(WSettings,ApplicationProtocol(iApp),[pathPrefix,'|Start time'],...
            'startTime',simulationIndex,studyDesignPathId,parValuesStudyDesign,isValid,isIndividualized);
        [ApplicationProtocol(iApp),isValid,isIndividualized] = addParameter(WSettings,ApplicationProtocol(iApp),[pathPrefix,'|Dose'],...
            'dose',simulationIndex,studyDesignPathId,parValuesStudyDesign,isValid,isIndividualized);
        [ApplicationProtocol(iApp),isValid,isIndividualized] = addParameter(WSettings,ApplicationProtocol(iApp),[pathPrefix,'|DosePerBodyWeight'],...
            'dosePerBodyWeight',simulationIndex,studyDesignPathId,parValuesStudyDesign,isValid,isIndividualized);
        [ApplicationProtocol(iApp),isValid,isIndividualized] = addParameter(WSettings,ApplicationProtocol(iApp),[pathPrefix,'|DosePerBodySurfaceArea'],...
            'dosePerBodySurfaceArea',simulationIndex,studyDesignPathId,parValuesStudyDesign,isValid,isIndividualized);
        
        % drugmass may overwrite formulacheck of dose, call it after the dose call
        [ApplicationProtocol(iApp),isValid,isIndividualized] = addParameter(WSettings,ApplicationProtocol(iApp),[pathPrefix,'|DrugMass'],...
            'drugMass',simulationIndex,studyDesignPathId,parValuesStudyDesign,isValid,isIndividualized);
        
        [ApplicationProtocol(iApp),isValid,isIndividualized] = addParameter(WSettings,ApplicationProtocol(iApp),[pathPrefix,'|Infusion time'],...
            'infusionTime',simulationIndex,studyDesignPathId,parValuesStudyDesign,isValid,isIndividualized);
                
        % get CompoundName
        tmp = getParameter(pathCompound,simulationIndex,'parametertype','readonly','property','path');
        if iscell(tmp) && length(tmp)>1
            jj = cellfun(@(x) ~contains(x,'InsolubleDrug'),tmp);
            tmp = tmp{jj};
        end
        tmp = regexp(tmp,'\|','split');
        ApplicationProtocol(iApp).compound = tmp{end-1};
        
    end
else
    writeToReportLog('WARNING',sprintf('There is no valid application within the xml %s',xml),false);
    isValid = false;
end    

if ~isValid
    return
end

% take only application within simulation range
simulationTime = getSimulationTime(simulationIndex);
jj = [ApplicationProtocol.startTime] < simulationTime(end);
ApplicationProtocol = ApplicationProtocol(jj);

% take only application with a Drugmass > 0
jj = any([ApplicationProtocol.drugMass] > 0,1);
ApplicationProtocol = ApplicationProtocol(jj);

% check if there is a application protocol
if isempty(ApplicationProtocol)
    writeToReportLog('WARNING',sprintf('There is no valid application within the xml %s',xml),false);
    isValid = false;
    return
end


% list application protocol in logfile
 writeToReportLog('INFO',sprintf('Applicationprotocol'),false);
for iApp = 1:length(ApplicationProtocol)
    
    if isIndividualized
        dosetxt = sprintf('individualized dose'); 
    else
        if ApplicationProtocol(iApp).isDosePerBodyweight
            dosetxt = sprintf('dose %g mg/kg',ApplicationProtocol(iApp).dosePerBodyWeight*1e6);
        elseif ApplicationProtocol(iApp).isDosePerBodySurfaceArea
            dosetxt = sprintf('dose %g mg/m^2',ApplicationProtocol(iApp).dosePerBodySurfaceArea*1e8);
        else
            dosetxt = sprintf('dose %g mg',ApplicationProtocol(iApp).dose*1e6);
        end
    end
    
    msg = sprintf('  %s: start time %g min; %s; infusion time %g min',ApplicationProtocol(iApp).name,ApplicationProtocol(iApp).startTime,...
        dosetxt,ApplicationProtocol(iApp).infusionTime);
    writeToReportLog('INFO',msg,false);
end

return


function ApplicationProtocol = getDefaultApplicationProtocol(name)

 ApplicationProtocol =    struct('name',name,'startTime',nan,...
        'dose',nan,'dosePerBodyWeight',nan,'dosePerBodySurfaceArea',nan,'drugMass',nan,...
        'isDosePerBodyweight',false,'isDosePerBodySurfaceArea',false,...
        'infusionTime',0,'compound','');
    
return    
    
 function [ApplicationProtocol,isValid,isIndividualized] = addParameter(~,ApplicationProtocol,pathID,...
         fieldNameStructure,simulationIndex,studyDesignPathId,parValuesStudyDesign,isValid,isIndividualized)
         
if existsParameter(pathID,simulationIndex,'parametertype','readonly')
    ApplicationProtocol.(fieldNameStructure) = getParameter(pathID,simulationIndex,'parametertype','readonly');

    % check formulas
    if strcmp(fieldNameStructure,'dose')
            formula = getParameter(pathID,simulationIndex,'parametertype','readonly','property','Formula');
            ApplicationProtocol.isDosePerBodyweight = contains(lower(formula),'doseperbodyweight');
            ApplicationProtocol.isDosePerBodySurfaceArea = contains(lower(formula),'doseperbodysurfacearea');
                
    end

    
    % check if pathID is part of studyDesign
    if ~isempty(studyDesignPathId)
        pathID = getParameter(pathID,simulationIndex,'parametertype','readonly','property','Path');
        for iPar = 1:length(studyDesignPathId)
            jj = strcmp(pathID,studyDesignPathId{iPar});
            if any(jj)
                isIndividualized = true;
                ApplicationProtocol.(fieldNameStructure) = parValuesStudyDesign(:,iPar);
        
                switch fieldNameStructure
                    case 'startTime'
                        writeToReportLog('WARNING','startTime is individualized by study design, applicationprotocol of xml is set to invalid',false);
                        isValid = false;
                    case 'drugMass'
                        % if drugmass is set individualized, dose ist absolute
                        ApplicationProtocol.isDosePerBodyweight = false;
                    otherwise
                        % no special action needed
                end
            end
        end
                
    end
else
    % check if mandatory field are missing
    if ismember(fieldNameStructure,{'startTime','dose','drugMass'})
                writeToReportLog('WARNING',sprintf('mandatory field %s is missing, applicationprotocol of xml is set to invalid',fieldNameStructure),false);
                isValid = false;
    end
            
end


return    
    
