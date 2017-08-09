function PKParameterTemplate = calculatePKParameterForApplicationProtocol(WSettings,ApplicationProtocol,simTime,simValues,weight,height)
%CALCULATEPKPARAMETERFORAPPLICATIONPROTOCOL calculates PK Parameter, same List as in PK-Sim, for multi or single application
%  for initialsiation of the availabel PK Parameters call function 
%
% PKParameterTemplate = calculatePKParameterForApplicationProtocol(WSettings,ApplicationProtocol);
%
% for calculation call
% 
% PKParameterTemplate = calculatePKParameterForApplicationProtocol(WSettings,simTime,simValues,weight,ApplicationProtocol)
% 
%Inputs: 
%   WSettings (structure)    definition of properties used in all
%                   workflow functions see GETDEFAULTWORKFLOWSETTINGS
%   ApplicationProtocol (structure with hold the information on the applicationprotocoll
%                   see GETAPPLICATIONPROTOCOLLFROMXML
%   simTime (double vector)    [min] vector of time      
%   simValues (double matrix)  [MoBi Units, e.g. µmol/l for
%               concentration],  values of time profile
%   weight (double vector)     [kg] vector body weight

% Open Systems Pharmacology Suite;  http://open-systems-pharmacology.org


% get start times of applications
[startingTimes] =  unique([ApplicationProtocol.startTime]);

% get flag
if length(startingTimes) >1
    flag = 'multiDose';
else
    flag = 'singleDose';
end

% get PKParameter template corresponding to flag
PKParameterTemplate = iniPKParameterTemplateList(flag);

% check if  initialiation
if nargin ==2
    return
end

%% start calculation

% get List of Intervals
[intervalList,~,ix_PK] = unique({PKParameterTemplate.timeRange});



% loop on Intervals
for iInt = 1:length(intervalList)
    
    switch intervalList{iInt}
        case 'total'
            timeRange = [startingTimes(1) simTime(end)];
            
        case 't1_t2'
            timeRange = [startingTimes(1) startingTimes(2)];

        case 'tLast_minus_1_tLast'
            timeRange = [startingTimes(end-1) startingTimes(end)];
            
        case 'tLast_tEnd'
            timeRange = [startingTimes(end) simTime(end)];
    end

    % get List of PK Parameter
    jj_PK = ix_PK == iInt;

    
    % get application of the range
    jjApp = timeRange(1) == [ApplicationProtocol.startTime]; 
%     jjApp =simTime(end) > [ApplicationProtocol.startTime]; toDo

    
    % sum over drugmass
    drugMass = zeros(size(weight));
    for iApp = find(jjApp)
        
        if ApplicationProtocol(iApp).isDosePerBodyweight
            BW_ref = ApplicationProtocol(iApp).dose./ApplicationProtocol(iApp).dosePerBodyWeight;
            drugMass = ApplicationProtocol(iApp).drugMass./BW_ref.*weight + drugMass;
        elseif ApplicationProtocol(iApp).isDosePerBodySurfaceArea
            BSA = calculateBodySurfaceArea(WSettings,weight,height);
            BSA_ref = ApplicationProtocol(iApp).dose./ApplicationProtocol(iApp).dosePerBodySurfaceArea;
            drugMass = ApplicationProtocol(iApp).drugMass./BSA_ref.*BSA + drugMass;
        else
            drugMass = drugMass +  ApplicationProtocol(iApp).drugMass.*ones(size(weight));
        end
    end

    % get dose per body weight used for factor calculation
    dose_per_weight = drugMass./weight; %#ok<NASGU>

    % get infusion time
    infusionTime = unique([ApplicationProtocol(jjApp).infusionTime]);
    if length(infusionTime)>1 
        if ismember('MRT',{PKParameterTemplate(jj_PK).fnPK});
            writeToLog(sprintf('WARNING for the time slot %d - %d infusion time was not unique MRT will be nan',...
                timeRange(1),timeRange(end)),WSettings.logfile,true,false);
        end
        infusionTime = nan;
    end

    % calculate pKparameter
    if isnan(infusionTime)
        PK = getPKParametersForConcentration(simTime,simValues,'Dose',drugMass,'timeRange',timeRange);
        PK.MRT = nan(size(PK.MRT));
    else
        PK = getPKParametersForConcentration(simTime,simValues,'Dose',drugMass,'infusionTime',infusionTime,'timeRange',timeRange);
    end
    ix_tend = find( simTime <= timeRange(2),1,'last');
    PK.cend = simValues(ix_tend,:);
    
    
    % transfer values
    for iPKP = find(jj_PK)'
        % get unitfactor
        f = eval(PKParameterTemplate(iPKP).unitFactor);
        
        PKParameterTemplate(iPKP).value = PK.(PKParameterTemplate(iPKP).fnPK).*f';
    end
end

% remove temporary field
PKParameterTemplate = rmfield(PKParameterTemplate,'unitFactor');
PKParameterTemplate = rmfield(PKParameterTemplate,'fnPK');
PKParameterTemplate = rmfield(PKParameterTemplate,'timeRange');

return

function PKParameterTemplate = iniPKParameterTemplateList(flag)
% get PKParameter template corresponding to flag
switch flag
    case 'singleDose'
        % PK Parameter for concentrations
        PKParameterTemplate = iniPKParameterTemplate('C_max','C_{max}','µmol/l','1','cMax','total');
        PKParameterTemplate(end+1) = iniPKParameterTemplate('C_max_norm','C_{max} norm','mg/l','1./dose_per_weight*1e6', 'cMax','total');
        PKParameterTemplate(end+1) = iniPKParameterTemplate('t_max','t_{max}','h','1/60','tMax','total');
        PKParameterTemplate(end+1) = iniPKParameterTemplate('C_tEnd','C_{tEnd}','µmol/l','1','cend','total');
        PKParameterTemplate(end+1) = iniPKParameterTemplate('AUC','AUC','µmol*min/l','1','AUC_0last','total');
        PKParameterTemplate(end+1) = iniPKParameterTemplate('AUC_norm','AUC_{last}','µg*min/l','1./dose_per_weight.*1e9','AUC_0last','total');
        PKParameterTemplate(end+1) = iniPKParameterTemplate('AUC_inf','AUC_{inf}','µmol*min/l','1','AUC_0inf','total');
        PKParameterTemplate(end+1) = iniPKParameterTemplate('AUC_inf_norm','AUC_{inf} norm','µg*min/l','1./dose_per_weight.*1e9','AUC_0inf','total');
        PKParameterTemplate(end+1) = iniPKParameterTemplate('MRT','MRT','h','1/60','MRT','total');
        PKParameterTemplate(end+1) = iniPKParameterTemplate('Thalf','Thalf','h','1/60','tHalf','total');
        PKParameterTemplate(end+1) = iniPKParameterTemplate('FractionAucLastToInf','Fraction AUC last to inf','','1/100','perc_AUC_lastToInf','total');
        PKParameterTemplate(end+1) = iniPKParameterTemplate('CL','CL','ml/min/kg','1000./weight','CL','total');
        PKParameterTemplate(end+1) = iniPKParameterTemplate('Vss','Vss','ml/kg','1000./weight','VSS','total');
        PKParameterTemplate(end+1) = iniPKParameterTemplate('Vd','Vd','ml/kg','1000./weight','Vd','total');
        
        % PK Parameter for fractions
        PKParameterTemplate(end+1) = iniPKParameterTemplate('F_max','F_{max}','','1','cMax','total');
        PKParameterTemplate(end+1) = iniPKParameterTemplate('F_tEnd','F_{max}','','1','cend','total');
        
        
    case 'multiDose';
        % PK Parameter for concentrations
        PKParameterTemplate = iniPKParameterTemplate('C_max','C_{max}','µmol/l','1','cMax','total');
        PKParameterTemplate(end+1) = iniPKParameterTemplate('C_max_norm','C_{max} norm','mg/l','1./dose_per_weight*1e6', 'cMax','total');
        PKParameterTemplate(end+1) = iniPKParameterTemplate('C_max_t1_t2','C_{max} (first dosing interval)','µmol/l','1','cMax','t1_t2');
        PKParameterTemplate(end+1) = iniPKParameterTemplate('C_max_t1_t2_norm','C_{max} norm (first dosing interval)',...
            'mg/l','1./dose_per_weight*1e6', 'cMax','t1_t2');
        PKParameterTemplate(end+1) = iniPKParameterTemplate('C_max_tLast_tEnd','C_{max} (last dosing interval)','µmol/l','1','cMax','tLast_tEnd');
        PKParameterTemplate(end+1) = iniPKParameterTemplate('C_max_tLast_tEnd_norm','C_{max} norm (last dosing interval',...
            'mg/l','1./dose_per_weight*1e6', 'cMax','tLast_tEnd');
        PKParameterTemplate(end+1) = iniPKParameterTemplate('t_max','t_{max}','h','1/60','tMax','total');
        PKParameterTemplate(end+1) = iniPKParameterTemplate('t_max_t1_t2','t_{max} (first dosing interval)','h','1/60','tMax','t1_t2');
        PKParameterTemplate(end+1) = iniPKParameterTemplate('t_max_tLast_tEnd','t_{max} (last dosing interval)','h','1/60','tMax','tLast_tEnd');
        PKParameterTemplate(end+1) = iniPKParameterTemplate('C_trough_t2','C_{trough} (first dosing interval)','µmol/l','1','cend','t1_t2'); 
        PKParameterTemplate(end+1) = iniPKParameterTemplate('C_trough_tLast','C_{trough} (last dosing interval)','µmol/l','1','cend','tLast_tEnd');
        PKParameterTemplate(end+1) = iniPKParameterTemplate('AUC_t1_t2','AUC_{last} (first dosing interval)','µmol*min/l','1','AUC_0last','t1_t2');
        PKParameterTemplate(end+1) = iniPKParameterTemplate('AUC_t1_t2_norm','AUC_{last} norm (first dosing interval)',...
            'µg*min/l','1./dose_per_weight.*1e9','AUC_0last','t1_t2');
        PKParameterTemplate(end+1) = iniPKParameterTemplate('AUC_tLast_minus_1_tLast','AUC_{last} (last dosing interval)',...
            'µmol*min/l','1','AUC_0last','tLast_minus_1_tLast');
        PKParameterTemplate(end+1) = iniPKParameterTemplate('AUC_tLast_minus_1_tLast_norm','AUC_{last} norm (first dosing interval)',...
            'µg*min/l','1./dose_per_weight.*1e9','AUC_0last','tLast_minus_1_tLast');
        PKParameterTemplate(end+1) = iniPKParameterTemplate('AUC_inf_t1','AUC_{inf} (first dosing interval)','µmol*min/l','1','AUC_0inf','t1_t2');
        PKParameterTemplate(end+1) = iniPKParameterTemplate('AUC_inf_t1_norm','AUC_{inf} norm (first dosing interval)',...
            'µg*min/l','1./dose_per_weight.*1e9','AUC_0inf','t1_t2');
        PKParameterTemplate(end+1) = iniPKParameterTemplate('AUC_inf_tLast','AUC_{inf} (last dosing interval)',...
            'µmol*min/l','1','AUC_0inf','tLast_tEnd');
        PKParameterTemplate(end+1) = iniPKParameterTemplate('AUC_inf_tLast_norm','AUC_{inf} norm (last dosing interval)',...
            'µg*min/l','1./dose_per_weight.*1e9','AUC_0inf','tLast_tEnd');
        PKParameterTemplate(end+1) = iniPKParameterTemplate('MRT','MRT (first dosing interval)','h','1/60','MRT','t1_t2');
        PKParameterTemplate(end+1) = iniPKParameterTemplate('Thalf','t_{1/2} (first dosing interval)','h','1/60','tHalf','t1_t2');
        PKParameterTemplate(end+1) = iniPKParameterTemplate('Thalf_tLast_tEnd','t_{1/2} (last dosing interval)','h','1/60','tHalf','tLast_tEnd');
        
        % PK Parameter for fractions
        PKParameterTemplate(end+1) = iniPKParameterTemplate('F_max','F_{max}','','1','cMax','total');
        PKParameterTemplate(end+1) = iniPKParameterTemplate('F_max_t1_t2','F_{max} (first dosing interval)','','1','cMax','t1_t2');
        PKParameterTemplate(end+1) = iniPKParameterTemplate('F_max_tLast_tEnd','F_{max} (last dosing interval)','','1','cMax','tLast_tEnd');
        PKParameterTemplate(end+1) = iniPKParameterTemplate('F_at_t2','F_{t end} (first dosing interval)','','1','cend','t1_t2'); 
        PKParameterTemplate(end+1) = iniPKParameterTemplate('F_at_tLast','F_{t end} (last dosing interval)','','1','cend','tLast_tEnd');

end

return

function  PKParameterTemplate = iniPKParameterTemplate(name,reportName,unit,unitFactor,fnPK,timeRange)
% name of Parameter (as in PKsim)
% report name used for display
% Unit as in PK-sim
% unitfactor from MoBi Units to PK-Sim Units (temporary field)
% field name of structure PK ( return value of
%       getPKParametersForConcentration) (temporary field)
% time range used for calculation (temporary field)
% placeholder for Values
PKParameterTemplate = struct('name',name,'reportName',reportName,'unit',unit,'unitFactor',unitFactor,...
    'fnPK',fnPK,'timeRange',timeRange,'value',[]);
    
return
