function [fig_handle, DDIRatioTable, DDIRatioQuali] = plotQualificationDDIRatio(WSettings,figureHandle,PKParameter,DDIRatioGroups,ObservedDataSets, SimulationMappings, AxesOptions, PlotSettings, REInputPath)
%PLOTQUALIFICATIONPKRATIO Plots PK ratio from qualification workflow
%
% plotQualificationDDIRatio(WSettings,figureHandle,PKParameter,DDIRatioGroups,ObservedDataSets, SimulationMappings, AxesOptions, PlotSettings)
%
% Inputs:
%       WSettings (structure)    definition of properties used in all
%                   workflow functions see GETDEFAULTWORKFLOWSETTINGS
%   figureHandle ( handle) handle of figure
%   PKParameter (string) name of the PK parameter to be evaluated
%   ObservedDataSets (structure) of all loaded observed data
%   SimulationMappings (structure) to map simulations with observations
%   AxesOptions (structure) to set axes options
%   PlotSettings (structure) to set plot options
% Output
%   csv (cellarray) table with numeric information to the plot

% Open Systems Pharmacology Suite;  http://open-systems-pharmacology.org

%---------------------------------------------
% Create figure with first setting from WSettings using getReportFigure
% To be updated using the Configuration plan Settings as optional arguments

% WARNING: FUNCTION STILL AS DRAFT: MANY FEATURES NOT AVAILABLE/HANDLED
% TO BE ADDED: LEGENDS, LABELS, GMFE, TABLE, ENSURE ROBUSTNESS

% xRatio and yRatio vectors for Guest et al. equation
% Load the template for the plot/ Formula not correct at the moment
xRatio=10.^(-2:0.01:2);
[yRatio1, yRatio2] = DDIRatioGuestEquation(xRatio);

for k=1:length(PKParameter)
    % create figure for Obs vs Pred
    [ax, fig_handle(k).predictedVsObserved] = getReportFigureQP(WSettings,1,1,2*k+figureHandle,PlotSettings);
    setFigureOptions(AxesOptions.DDIRatioPlotsPredictedVsObserved);
    figure(fig_handle(k).predictedVsObserved);
    plot(xRatio, xRatio, '--k', 'Linewidth', 1, 'HandleVisibility','off');
    plot(xRatio, xRatio/2, '--k', 'Linewidth', 1, 'HandleVisibility','off');
    plot(xRatio, xRatio*2, '--k', 'Linewidth', 1, 'HandleVisibility','off');
    plot(xRatio, yRatio1, '--k', 'Linewidth', 1, 'HandleVisibility','off');
    plot(xRatio, yRatio2, '--k', 'Linewidth', 1, 'HandleVisibility','off');
    xlabel('Observed Ratio'); ylabel('Predicted Ratio');
    
    % create figure for Residuals
    [ax, fig_handle(k).residualsVsObserved] = getReportFigureQP(WSettings,1,1,2*k+figureHandle+1,PlotSettings);
    setFigureOptions(AxesOptions.DDIRatioPlotsResidualsVsObserved);
    figure(fig_handle(k).residualsVsObserved);
    plot(xRatio, ones(size(xRatio)), '--k', 'Linewidth', 1, 'HandleVisibility','off');
    plot(xRatio, ones(size(xRatio))/2, '--k', 'Linewidth', 1, 'HandleVisibility','off');
    plot(xRatio, ones(size(xRatio))*2, '--k', 'Linewidth', 1, 'HandleVisibility','off');
    plot(xRatio, yRatio1./xRatio, '--k', 'Linewidth', 1, 'HandleVisibility','off');
    plot(xRatio, yRatio2./xRatio, '--k', 'Linewidth', 1, 'HandleVisibility','off');
    xlabel('Observed Ratio'); ylabel('Residuals');
end

% Initialize error for computing GMFE
leg_labels={};
Error=[];

% Loop on the Ratios to be plotted by PK Ratio plot
for i=1:length(DDIRatioGroups)
    
    DDIRatios=DDIRatioGroups(i).DDIRatios;
    
    for j=1:length(DDIRatios)
        
        % Get the observed data
        for k=1:length(ObservedDataSets)
            if strcmp(DDIRatios(j).ObservedData, ObservedDataSets(k).Id)
                ObservedData = ObservedDataSets(k).y;
                break
            end
        end
        Results(i).Study(j)=table2array(ObservedData(ObservedData.ID==DDIRatios(j).ObservedDataRecordId,'StudyID'));
        
        % Load the mapped Time Profile Simulation Results
        [csvSimFileControl, xmlfileControl] = getSimFile(DDIRatios(j).SimulationControl, SimulationMappings, REInputPath);
        if isempty(csvSimFileControl)
            ME = MException('plotQualificationDDIRatio:notFoundInPath', ...
                'In DDI Ratio plot %d, Group %d, Ratio %d, Project "%s" or Simulation "%s" for Control was not found in SimulationMappings', figureHandle, i, j, DDIRatios(j).SimulationControl.Project, DDIRatios(j).SimulationControl.Simulation);
            throw(ME);
        end
        SimResultControl = loadSimResultcsv(csvSimFileControl, DDIRatios(j).SimulationControl.Simulation);
        
        % All the output are kept so far, may be removed if not necessary
        [AGEControl, BWControl, MWControl, drugmassControl] = getInfofromSimulation(xmlfileControl, DDIRatios(j).Output);
        
        Result(i).AGEControl(j)=AGEControl;
        Result(i).BWControl(j)=BWControl;
        Result(i).MWControl(j)=MWControl;
        Result(i).drugmassControl(j)=drugmassControl(1);
        
        [csvSimFileDDI, xmlfileDDI] = getSimFile(DDIRatios(j).SimulationDDI, SimulationMappings, REInputPath);
        if isempty(csvSimFileDDI)
            ME = MException('plotQualificationDDIRatio:notFoundInPath', ...
                'In DDI Ratio plot %d, Group %d, Ratio %d, Project "%s" or Simulation "%s" for DDI was not found in SimulationMappings', figureHandle, i, j, DDIRatios(j).SimulationDDI.Project, DDIRatios(j).SimulationDDI.Simulation);
            throw(ME);
        end
        SimResultDDI = loadSimResultcsv(csvSimFileDDI, DDIRatios(j).SimulationDDI.Simulation);
        
        [AGEDDI, BWDDI, MWDDI, drugmassDDI] = getInfofromSimulation(xmlfileDDI, DDIRatios(j).Output);
        
        Result(i).AGEDDI(j)=AGEDDI;
        Result(i).BWDDI(j)=BWDDI;
        Result(i).MWDDI(j)=MWDDI;
        Result(i).drugmassDDI(j)=drugmassControl(1);
        
        % Use the right Output for Control
        for k=1:length(SimResultControl.outputPathList)
            findPathOutput = strfind(SimResultControl.outputPathList{k}, DDIRatios(j).Output);
            if ~isempty(findPathOutput)
                Controlpred=SimResultControl.y{k};
                ControlTime=SimResultControl.time;
                ControlTimeUnit=SimResultControl.timeUnit;
                
                % Only get output in a time range defined by user
                Xfactor=getUnitFactor(ControlTimeUnit,DDIRatios(j).SimulationControl.TimeUnit,'time');
                
                SimTime = (ControlTime.*Xfactor >= DDIRatios(j).SimulationControl.StartTime &  ControlTime.*Xfactor <= DDIRatios(j).SimulationControl.EndTime);
                ControlTime = ControlTime(SimTime).*Xfactor;
                Controlpred = Controlpred(SimTime);
                
                break
            end
        end
        
        % Use the right Output for DDI
        for k=1:length(SimResultDDI.outputPathList)
            findPathOutput = strfind(SimResultDDI.outputPathList{k}, DDIRatios(j).Output);
            if ~isempty(findPathOutput)
                DDIpred=SimResultDDI.y{k};
                DDITime=SimResultDDI.time;
                DDITimeUnit=SimResultDDI.timeUnit;
                
                % Only get output in a time range defined by user
                Xfactor=getUnitFactor(DDITimeUnit,DDIRatios(j).SimulationDDI.TimeUnit,'time');
                
                SimTime = (DDITime.*Xfactor >= DDIRatios(j).SimulationDDI.StartTime &  DDITime.*Xfactor <= DDIRatios(j).SimulationDDI.EndTime);
                DDITime = DDITime(SimTime).*Xfactor;
                DDIpred = DDIpred(SimTime);
                break
            end
        end
        
        % Get the PK parameters out of the simulation
        allPKpredControl=getPKParametersForConcentration(ControlTime,Controlpred,'Dose', drugmassControl);
        
        % Get the PK parameters out of the simulation
        allPKpredDDI=getPKParametersForConcentration(DDITime,DDIpred,'Dose', drugmassDDI);
        
        
        for k=1:length(PKParameter)
            % Get the PK parameters requested in PKParameter
            if strcmpi(PKParameter{k}, 'AUC')
                PKpredField{k}= 'AUC_last';
            elseif strcmpi(PKParameter{k}, 'CMAX')
                PKpredField{k}= 'cMax';
            else
                PKpredField{k}=PKParameter{k};
            end
            if isfield(allPKpredControl, PKpredField{k}) && isfield(allPKpredDDI, PKpredField{k})
                
                
                PKRobs(j,k) = table2array(ObservedData(ObservedData.ID==DDIRatios(j).ObservedDataRecordId,...
                    (strcmpi(ObservedData.Properties.VariableNames, [PKParameter{k} 'RAvg']))));
                
                Result(i).DDIPK(j,k)=getfield(allPKpredDDI, PKpredField{k});
                
                Result(i).ControlPK(j,k)=getfield(allPKpredControl, PKpredField{k});
                
                % Assumes currently that the output from PKSim has same unit for Control
                % and DDI
                Result(i).RatioPK(j,k)=Result(i).DDIPK(j,k)./Result(i).ControlPK(j,k);
                
                figure(fig_handle(k).predictedVsObserved);
                pp=plot(PKRobs(j,k), Result(i).RatioPK(j,k), 'o', 'Linewidth',1);
                setCurveOptions(pp, DDIRatioGroups(i));
                
                figure(fig_handle(k).residualsVsObserved);
                pp=plot(PKRobs(j,k), Result(i).RatioPK(j,k)./PKRobs(j,k), 'o', 'Linewidth',1);
                setCurveOptions(pp, DDIRatioGroups(i));
                
                if isfield(DDIRatioGroups(i), 'Caption')
                leg_labels=[leg_labels DDIRatioGroups(i).Caption];
                end
                
            else
                ME = MException('DDIRatio:notFoundInField', ...
                    'Requested PK Parameter %s not found in parameters extracted using getPKParametersForConcentration', PKParameter{k});
                throw(ME);
            end
        end
    end
end

for k=1:length(PKParameter)
    figure(fig_handle(k).predictedVsObserved);
    legend(leg_labels);
    figure(fig_handle(k).predictedVsObserved);
    legend(residualsVsObserved);
end

% TABLES TO BE CALCULATED AND UPDATED
% Get the DDI Ratio Table
DDIRatioHeader = {'Perpetrator', 'Victim', 'Dose gap', 'Males', ...
    'Predicted AUC Ratio', 'Observed AUC Ratio', 'Predicted/Obs AUC Ratio', ...
    'Predicted Cmax Ratio', 'Observed Cmax Ratio', 'Predicted/Obs AUC Ratio', 'Reference'};

DDIRatioTable = [DDIRatioHeader];

disp(DDIRatioTable);
    
% Get the DDI Ratio Qualification
DDIRatioQualiHeader = {'', 'Number', 'Ratio'};
DDIRatioQuali_1st_Column = {'Points total'; 'Points within Guest et al.'; 'Points within 2-fold'};

DDIRatioQuali = [DDIRatioQualiHeader  ; ...
    DDIRatioQuali_1st_Column num2cell(zeros(3,2))];


disp(DDIRatioQuali);


function [AGE, BW, MW, drugmass] = getInfofromSimulation(xmlfile, Output)

initSimulation(xmlfile,'none');

MW = getMolecularWeightForPath(Output);

drugmass = getParameter('*Application_*|ProtocolSchemaItem|DrugMass',1,'parametertype','readonly');

AGE = getParameter('*|Organism|Age',1,'parametertype','readonly');
BW = getParameter('*|Organism|Weight',1,'parametertype','readonly');

function [UpperLimit, LowerLimit] = DDIRatioGuestEquation(Robs, delta)
% Upper and Lower limits for DDI Ratio plots as proposed by Guest et al.


% Delta might be different from 1 to use different limits
if ~exist('delta')
    delta=1;
end

symRobs=Robs;

symRobs(Robs<1)=1./symRobs(symRobs<1);

Limit = (delta + 2.*(symRobs-1))./symRobs;

UpperLimit = symRobs.*Limit;
LowerLimit = symRobs./Limit;

UpperLimit(Robs<1)=1./UpperLimit(Robs<1);
LowerLimit(Robs<1)=1./LowerLimit(Robs<1);