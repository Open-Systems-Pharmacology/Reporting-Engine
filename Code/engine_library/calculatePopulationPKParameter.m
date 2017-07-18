function calculatePopulationPKParameter(Settings,PopRunSet)
% CALCULATEPOPULATIONPKPARAMETER calculates PK-Parameter for simulations defined within a PopRunSet 
% see GETDEFAULTPOPRUNSET
%
% Inputs:
%   Settings  structure containing global settings see GETDEFAULTWORKFLOWSETTINGS
%   PopRunSet  defines properties and corresponding files of a simulation
%           see GETDEFAULTPOPRUNSET
 

% Open Systems Pharmacology Suite;  support@systems-biology.com
% Date: 14-July-2017

writeToLog(sprintf('Calculate PK Parameter of %s',PopRunSet.name),Settings.logfile,true,false);

% analyse applicationProtocol
[ApplicationProtocol,isValid] = getApplicationProtocollFromXML(Settings,PopRunSet.xml);  

if isValid
    writeToLog(sprintf('Applicationprotocol'),Settings.logfile,true,false);
    for iApp = 1:length(ApplicationProtocol)
        if ApplicationProtocol(1).isDosePerBodyweight
            msg = sprintf('  start time %g min; dose %g g/kg; infusion time %g min',ApplicationProtocol(iApp).startTime,...
                ApplicationProtocol(iApp).dosePerBodyWeight*1e6,ApplicationProtocol(iApp).infusionTime);
        else
            msg = sprintf('  start time %g min; dose %g g; infusion time %g min',ApplicationProtocol(iApp).startTime,...
                ApplicationProtocol(iApp).dose*1e6,ApplicationProtocol(iApp).infusionTime);
        end
        writeToLog(msg,Settings.logfile,true,false);
    end
end

% read population csv
[parPathes,parValues] = readPopulationCSV(PopRunSet.pop);
% reduce to demographic
[jj,ix] = ismember({'IndividualId','Organism|Weight','Organism|Height','Organism|BMI'},parPathes);
parPathes = parPathes(ix(jj));
parValues = parValues(:,ix(jj));

% load time profiles
SimResult = readPopulationResultfile(PopRunSet.name);


% add studyDesign if one is given
if ~isempty(PopRunSet.studyDesign)
    error('not implemented yet') %ToDo
end

% initialize PK List
PKPList = cell(length(SimResult.outputList) ,1);

% decide on the PKParameter function and calculate PKParameter
% if given take the projectspecific function
if ~isempty(PopRunSet.calculatePKParameter_fh)
    
elseif isValid
    
    jj_BW = strcmp('Organism|Weight',parPathes);
    weight =  parValues(:,jj_BW);
    
    % add studyDesign if one is given
    if ~isempty(PopRunSet.studyDesign)
        error('not implemented yet') %ToDo
    end
    
    for iO = 1:length(SimResult.outputList)
        PKPList{iO} = calculatePKParameter_forApplicationProtocol(Settings,SimResult.time,SimResult.values{iO},weight,ApplicationProtocol);
    end

else            
    error('application protocol could not be interpreted. Please use project specific function for PK Parameter calculation');
end

% export PKParameter
exportPKParameter(PKPList,SimResult);

writeToLog(sprintf('Calculation finished \n'),Settings.logfile,true,false);

return

function exportPKParameter(PKPList,SimResult)
% get name of resultfile
resultfile =fullfile('Simulations',[SimResult.name '-PK-Analyses.csv']); 

% check if siluation directory already exist
if ~exist('Simulations','dir')
    mkdir('Simulations')
end

% Write header
fid = fopen(resultfile,'w');
fprintf(fid,'%s','IndividualId;Quantity Path;Parameter;Value;Unit');
fprintf(fid,'\r\n');


% Write values
for iO = 1:length(PKPList)
    PKParameter = PKPList{iO};
    for iPKP = 1:length(PKParameter)
        for iInd = 1:length(SimResult.individualIdVector)
            linetxt = sprintf('%d;%s;%s;%g;%s',SimResult.individualIdVector(iInd),SimResult.outputList{iO},...
                PKParameter(iPKP).name,PKParameter(iPKP).value(iInd),PKParameter(iPKP).unit);
            fprintf(fid,'%s',linetxt);
            fprintf(fid,'\r\n');
        end
    end
end

fclose(fid);

return