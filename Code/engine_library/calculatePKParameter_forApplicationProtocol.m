function PKParameterTemplate = calculatePKParameter_forApplicationProtocol(Settings,sim_time,sim_values,weight,ApplicationProtocol)
%CALCULATEPKPARAMETER_FORAPPLICATIONPROTOCOL calculates PK Parameter, same List as in PK-Sim, for multi or single application
%
%Inputs: 
%   Settings (structure)    definition of properties used in all
%                   workflow functions see GETDEFAULTWORKFLOWSETTINGS
%   sim_time (double vector)    [min] vector of time      
%   sim_values (double matrix)  [MoBi Units, e.g. �mol/l for
%               concentration],  values of time profile
%   weight (double vector)     [kg] vector body weight
%   ApplicationProtocol (structure with hold the information on the applicationprotocoll
%                   see GETAPPLICATIONPROTOCOLLFROMXML

% Open Systems Pharmacology Suite;  support@systems-biology.com
% Date: 14-July-2017

% get start times of applications
[startingTimes] =  unique([ApplicationProtocol.startTime]);

% get flag
if length(startingTimes) >1
    flag = 'multiDose';
else
    flag = 'singleDose';
end

% get PKParameter template corresponding to flag
switch flag
    case 'singleDose'
        PKParameterTemplate = iniPKParameterTemplate('C_max','�mol/l','1','cMax','total');
        PKParameterTemplate(end+1) = iniPKParameterTemplate('C_max_norm','mg/l','1./dose_per_weight*1e6', 'cMax','total');
        PKParameterTemplate(end+1) = iniPKParameterTemplate('t_max','h','1/60','tMax','total');
        PKParameterTemplate(end+1) = iniPKParameterTemplate('C_tEnd','�mol/l','1','cend','total');
        PKParameterTemplate(end+1) = iniPKParameterTemplate('AUC','�mol*min/l','1','AUC_0last','total');
        PKParameterTemplate(end+1) = iniPKParameterTemplate('AUC_norm','�g*min/l','1./dose_per_weight.*1e9','AUC_0last','total');
        PKParameterTemplate(end+1) = iniPKParameterTemplate('AUC_inf','�mol*min/l','1','AUC_0inf','total');
        PKParameterTemplate(end+1) = iniPKParameterTemplate('AUC_inf_norm','�g*min/l','1./dose_per_weight.*1e9','AUC_0inf','total');
        PKParameterTemplate(end+1) = iniPKParameterTemplate('MRT','h','1/60','MRT','total');
        PKParameterTemplate(end+1) = iniPKParameterTemplate('Thalf','h','1/60','tHalf','total');
        PKParameterTemplate(end+1) = iniPKParameterTemplate('FractionAucLastToInf','','1/100','perc_AUC_lastToInf','total');
        PKParameterTemplate(end+1) = iniPKParameterTemplate('CL','ml/min/kg','1000./weight','CL','total');
        PKParameterTemplate(end+1) = iniPKParameterTemplate('Vss','ml/kg','1000./weight','VSS','total');
        PKParameterTemplate(end+1) = iniPKParameterTemplate('Vd','ml/kg','1000./weight','Vd','total');
        
    case 'multiDose';
        PKParameterTemplate = iniPKParameterTemplate('C_max','�mol/l','1','cMax','total');
        PKParameterTemplate(end+1) = iniPKParameterTemplate('C_max_norm','mg/l','1./dose_per_weight*1e6', 'cMax','total');
        PKParameterTemplate(end+1) = iniPKParameterTemplate('C_max_t1_t2','�mol/l','1','cMax','t1_t2');
        PKParameterTemplate(end+1) = iniPKParameterTemplate('C_max_t1_t2_norm','mg/l','1./dose_per_weight*1e6', 'cMax','t1_t2');
        PKParameterTemplate(end+1) = iniPKParameterTemplate('C_max_tLast_tEnd','�mol/l','1','cMax','tLast_tEnd');
        PKParameterTemplate(end+1) = iniPKParameterTemplate('C_max_tLast_tEnd_norm','mg/l','1./dose_per_weight*1e6', 'cMax','tLast_tEnd');
        PKParameterTemplate(end+1) = iniPKParameterTemplate('t_max','h','1/60','tMax','total');
        PKParameterTemplate(end+1) = iniPKParameterTemplate('t_max_t1_t2','h','1/60','tMax','t1_t2');
        PKParameterTemplate(end+1) = iniPKParameterTemplate('t_max_tLast_tEnd','h','1/60','tMax','tLast_tEnd');
        PKParameterTemplate(end+1) = iniPKParameterTemplate('C_trough_t2','�mol/l','1','cend','t1_t2'); 
        PKParameterTemplate(end+1) = iniPKParameterTemplate('C_trough_tLast','�mol/l','1','cend','tLast_tEnd');
        PKParameterTemplate(end+1) = iniPKParameterTemplate('AUC_t1_t2','�mol*min/l','1','AUC_0last','t1_t2');
        PKParameterTemplate(end+1) = iniPKParameterTemplate('AUC_t1_t2_norm','�g*min/l','1./dose_per_weight.*1e9','AUC_0last','t1_t2');
        PKParameterTemplate(end+1) = iniPKParameterTemplate('AUC_tLast_minus_1_tLast','�mol*min/l','1','AUC_0last','tLast_minus_1_tLast');
        PKParameterTemplate(end+1) = iniPKParameterTemplate('AUC_tLast_minus_1_tLast_norm','�g*min/l','1./dose_per_weight.*1e9','AUC_0last','tLast_minus_1_tLast');
        PKParameterTemplate(end+1) = iniPKParameterTemplate('AUC_inf_t1','�mol*min/l','1','AUC_0inf','t1_t2');
        PKParameterTemplate(end+1) = iniPKParameterTemplate('AUC_inf_t1_norm','�g*min/l','1./dose_per_weight.*1e9','AUC_0inf','t1_t2');
        PKParameterTemplate(end+1) = iniPKParameterTemplate('AUC_inf_tLast','�mol*min/l','1','AUC_0inf','tLast_tEnd');
        PKParameterTemplate(end+1) = iniPKParameterTemplate('AUC_inf_tLast_norm','�g*min/l','1./dose_per_weight.*1e9','AUC_0inf','tLast_tEnd');
        PKParameterTemplate(end+1) = iniPKParameterTemplate('MRT','h','1/60','MRT','t1_t2');
        PKParameterTemplate(end+1) = iniPKParameterTemplate('Thalf','h','1/60','tHalf','t1_t2');
        PKParameterTemplate(end+1) = iniPKParameterTemplate('Thalf_tLast_tEnd','h','1/60','tHalf','tLast_tEnd');
        
end

% get List of Intervals
[intervalList,~,ix_PK] = unique({PKParameterTemplate.timeRange});



% loop on Intervals
for iInt = 1:length(intervalList)
    
    switch intervalList{iInt}
        case 'total'
            timeRange = [startingTimes(1) sim_time(end)];
            
        case 't1_t2'
            timeRange = [startingTimes(1) startingTimes(2)];

        case 'tLast_minus_1_tLast'
            timeRange = [startingTimes(end-1) startingTimes(end)];
            
        case 'tLast_tEnd'
            timeRange = [startingTimes(end) sim_time(end)];
    end

    % get application of the range
%     jjApp = timeRange(1) == [ApplicationProtocol.startTime]; % toDo
    jjApp =sim_time(end) > [ApplicationProtocol.startTime];

    
    % sum over drugmass
    drugMass = zeros(size(weight));
    for iApp = find(jjApp)
        
        if ApplicationProtocol(iApp).isDosePerBodyweight
            BW_ref = ApplicationProtocol(iApp).dose./ApplicationProtocol(iApp).dosePerBodyWeight;
            drugMass = ApplicationProtocol(iApp).drugMass./BW_ref.*weight + drugMass;
        else
            drugMass = drugMass +  ApplicationProtocol(iApp).drugMass.*ones(size(weight));
        end
    end

    % get dose per body weight used for factor calculation
    dose_per_weight = drugMass./weight; %#ok<NASGU>

    
    infusionTime = unique([ApplicationProtocol(jjApp).infusionTime]);
    if length(infusionTime)>1
        writeToLog('WARNING for the time slot %d - %d infusion time was not unique MRT will be nan',...
            timeRange(1),timeRange(end),Settings.logfile,true,false);
        infusionTime = nan;
    end

    

    % calculate pKparameter
    PK = getPKParametersForConcentration(sim_time,sim_values,'Dose',drugMass,'infusionTime',infusionTime,'timeRange',timeRange);
    ix_tend = find( sim_time <= timeRange(2),1,'last');
    PK.cend = sim_values(ix_tend,:);
    
    % get List of PK Parameter
    jj_PK = ix_PK == iInt;
    
    % transfer values
    for iPKP = find(jj_PK)'
        % get unitfactor
        f = eval(PKParameterTemplate(iPKP).unitFactor);
        
        PKParameterTemplate(iPKP).value = PK.(PKParameterTemplate(iPKP).fnPK).*f';
    end
end

% remove temporyr field
PKParameterTemplate = rmfield(PKParameterTemplate,'unitFactor');
PKParameterTemplate = rmfield(PKParameterTemplate,'fnPK');
PKParameterTemplate = rmfield(PKParameterTemplate,'timeRange');

return



function  PKParameterTemplate = iniPKParameterTemplate(name,unit,unitFactor,fnPK,timeRange)
% 1: name of Parameter (as in PKsim)
% 2. Unit as in PK-sim
% 3. unitfactor from MoBi Units to PK-Sim Units
% 4. field name of structure PK ( return value of getPKParametersForConcentration)
% 5. time range used for calculation
% 6. placeholder for Values
PKParameterTemplate = struct('name',name,'unit',unit,'unitFactor',unitFactor,...
    'fnPK',fnPK,'timeRange',timeRange,'values',[]);
    
return
