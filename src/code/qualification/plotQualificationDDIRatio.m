function fig_handle = plotQualificationDDIRatio(WSettings,figureHandle,PKParameter,DDIRatioGroups,ObservedDataSets, SimulationMappings, AxesOptions, PlotSettings)
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
    [ax, fig_handle(k).predictedVsObserved] = getReportFigureQP(WSettings,1,1,figureHandle,PlotSettings);
    setFigureOptions(AxesOptions.DDIRatioPlotsPredictedVsObserved);
    figure(fig_handle(k).predictedVsObserved);
    plot(xRatio, xRatio, '--k', 'Linewidth', 1, 'HandleVisibility','off');
    plot(xRatio, xRatio/2, '--k', 'Linewidth', 1, 'HandleVisibility','off');
    plot(xRatio, xRatio*2, '--k', 'Linewidth', 1, 'HandleVisibility','off');
    plot(xRatio, yRatio1, '--k', 'Linewidth', 1, 'HandleVisibility','off');
    plot(xRatio, yRatio2, '--k', 'Linewidth', 1, 'HandleVisibility','off');
    xlabel('Observed Ratio'); ylabel('Predicted Ratio');
    
    % create figure for Residuals
    [ax, fig_handle(k).residualsVsObserved] = getReportFigureQP(WSettings,1,1,figureHandle+1,PlotSettings);
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
        [csvSimFileControl, xmlfileControl] = getSimFile(DDIRatios(j).SimulationControl, SimulationMappings);
        SimResultControl = loadSimResultcsv(csvSimFileControl, DDIRatios(j).SimulationControl.Simulation);
        
        % All the output are kept so far, may be removed if not necessary
        [AGEControl, BWControl, MWControl, drugmassControl] = getInfofromSimulation(xmlfileControl, DDIRatios(j).Output);
        
        Result(i).AGEControl(j)=AGEControl;
        Result(i).BWControl(j)=BWControl;
        Result(i).MWControl(j)=MWControl;
        Result(i).drugmassControl(j)=drugmassControl(1);
        
        [csvSimFileDDI, xmlfileDDI] = getSimFile(DDIRatios(j).SimulationDDI, SimulationMappings);
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
                
                % To be implmented:
                % only get output in a time range defined by user
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
                
                Result(i).RatioPK(j,k)=Result(i).DDIPK(j,k)./Result(i).ControlPK(j,k);
                
                figure(fig_handle(k).predictedVsObserved);
                pp=plot(PKRobs(j,k), Result(i).RatioPK(j,k), 'o', 'Linewidth',1);
                
                figure(fig_handle(k).residualsVsObserved);
                pp=plot(PKRobs(j,k), Result(i).RatioPK(j,k)./PKRobs(j,k), 'o', 'Linewidth',1);
                %setCurveOptions(pp, CurveOptions);
                
            else
                ME = MException('DDIRatio:notFoundInField', ...
                    'Requested PK Parameter %s not found in parameters extracted using getPKParametersForConcentration', PKParameter{k});
                throw(ME);
            end
        end
    end
end

%{
% Definition of the PKRatio table
PKRatioHeader = {'StudyID', 'Age (y)', 'BodyWeight (kg)', ['Report' PKParameter], ['PKSim' PKParameter], 'Ratio'};
PKRatioTable = [PKRatioHeader ; ...
    Study' num2cell([AGE', BW', PKobs', PKpred', PKpred'./PKobs'])];

disp(PKRatioTable);

% Calculation of GMFE
GMFE = 10.^(sum(abs(log(PKpred./PKobs)))/length(PKobs));
fprintf('GMFE = %f \n', GMFE);

% Perform the plot based on Curves indications
Xrange=[0.8*min(AGE) 1.2*max(AGE)]; Yrange=[1 1];

plot(Xrange, Yrange, '-k', 'LineWidth', 1);
plot(Xrange, Yrange/2, '--r', 'LineWidth', 1); plot(Xrange, Yrange*2, '--r', 'LineWidth', 1);
plot(Xrange, Yrange/1.5, '--b', 'LineWidth', 1); plot(Xrange, Yrange*1.5, '--b', 'LineWidth', 1);

ylim([min([0.1 PKpred./PKobs]) max([10 PKpred./PKobs])]);
ylabel(sprintf('PK ratio %s Pred / %s Obs', PKParameter, PKParameter));
% Placement of legend without masking the plot has not been figure out
% Update to be perform
legend('off')
%}

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