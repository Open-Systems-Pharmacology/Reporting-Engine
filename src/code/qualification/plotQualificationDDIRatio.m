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

% Initialize error for computing GMFE
leg_labels={};
DDIRatioTableContent={};

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
        Result(i).drugmassDDI(j)=drugmassDDI(1);
        
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
                
                % Get the observation Ratios
                Observations(i).RatioPK(k,j) = table2array(ObservedData(ObservedData.ID==DDIRatios(j).ObservedDataRecordId,...
                    (strcmpi(ObservedData.Properties.VariableNames, [PKParameter{k} 'RAvg']))));
                
                Result(i).DDIPK(k,j)=getfield(allPKpredDDI, PKpredField{k});
                
                Result(i).ControlPK(k,j)=getfield(allPKpredControl, PKpredField{k});
                
                % Assumes currently that the output from PKSim has same unit for Control and DDI
                Result(i).RatioPK(k,j)=Result(i).DDIPK(k,j)./Result(i).ControlPK(k,j);
                
            else
                ME = MException('DDIRatio:notFoundInField', ...
                    'Requested PK Parameter %s not found in parameters extracted using getPKParametersForConcentration', PKParameter{k});
                throw(ME);
            end
        end
        % Build the DDI Ratio Table:
        Perpetrator = table2cell(ObservedData(ObservedData.ID==DDIRatios(j).ObservedDataRecordId,...
            {'Perpetrator', 'Dose', 'DoseUnit', 'RoutePerpetrator'}));
        Victim = table2cell(ObservedData(ObservedData.ID==DDIRatios(j).ObservedDataRecordId, {'Victim', 'RouteVictim'}));
        Reference = table2cell(ObservedData(ObservedData.ID==DDIRatios(j).ObservedDataRecordId, {'StudyID'}));
        
        % Reshape the ratio table as a line Pred Obs Pred/Obs
        DDIRatioLinePK=[];
        for k=1:length(PKParameter)
            DDIRatioLinePK = [DDIRatioLinePK Result(i).RatioPK(k,j), Observations(i).RatioPK(k,j), Result(i).RatioPK(k,j)./Observations(i).RatioPK(k,j)];
        end
        
        DDIRatioLine = [{sprintf('%s, %s %s, %s', Perpetrator{:}), sprintf('%s, %s', Victim{:})},...
            num2cell(DDIRatioLinePK), {sprintf('%s', Reference{:})}];
        
        DDIRatioTableContent = [DDIRatioTableContent; DDIRatioLine];
        
    end
end

% xRatio and yRatio vectors for Guest et al. equation
% Load the template for the plot/ Formula not correct at the moment
xRatio=10.^(-1.5:0.01:1.5);
[yRatio1, yRatio2] = DDIRatioGuestEquation(xRatio);

for k=1:length(PKParameter)
    % Initialize the Qualification Measures
    QualiMeasure(k).PointsTotal = 0;
    QualiMeasure(k).PointsGuest= 0;
    QualiMeasure(k).Points2fold= 0;
end

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
    
    for i=1:length(DDIRatioGroups)
        % Update the number of points for each group
        conditionPoints = ~isnan(Observations(i).RatioPK(k,:));
        QualiMeasure(k).PointsTotal = QualiMeasure(k).PointsTotal + length(Observations(i).RatioPK(k,conditionPoints));
        [UpperBound, LowerBound] = DDIRatioGuestEquation(Observations(i).RatioPK(k,:));
        conditionGuest = (Result(i).RatioPK(k,:) >= LowerBound & Result(i).RatioPK(k,:) <= UpperBound);
        QualiMeasure(k).PointsGuest= QualiMeasure(k).PointsGuest + length(Result(i).RatioPK(k,conditionGuest));
        condition2fold = (Result(i).RatioPK(k,:) >= Observations(i).RatioPK(k,:)/2 & Result(i).RatioPK(k,:) <= Observations(i).RatioPK(k,:)*2);
        QualiMeasure(k).Points2fold= QualiMeasure(k).Points2fold+ length(Result(i).RatioPK(k,condition2fold));
        
        % Plot part
        figure(fig_handle(k).predictedVsObserved);
        pp=plot(Observations(i).RatioPK(k,:), Result(i).RatioPK(k,:), 'o', 'Linewidth',1);
        setCurveOptions(pp, DDIRatioGroups(i));
        
        figure(fig_handle(k).residualsVsObserved);
        pp=plot(Observations(i).RatioPK(k,:), Result(i).RatioPK(k,:)./Observations(i).RatioPK(k,:), 'o', 'Linewidth',1);
        setCurveOptions(pp, DDIRatioGroups(i));
        
        if isfield(DDIRatioGroups(i), 'Caption')
            leg_labels=[leg_labels DDIRatioGroups(i).Caption];
        end
    end
    figure(fig_handle(k).predictedVsObserved);
    legend(leg_labels, 'Location', 'northoutside');
    figure(fig_handle(k).residualsVsObserved);
    legend(leg_labels, 'Location', 'northoutside');
end

% Get the DDI Ratio Table
DDIRatioHeaderPK={};
for k=1:length(PKParameter)
    DDIRatioHeaderPK = [DDIRatioHeaderPK {sprintf('Predicted %s Ratio', PKParameter{k}), ...
        sprintf('Observed %s Ratio', PKParameter{k}), sprintf('Pred/Obs %s Ratio', PKParameter{k})}];
end
DDIRatioHeader = [{'Perpetrator', 'Victim'}, DDIRatioHeaderPK, {'Reference'}];

DDIRatioTable = [DDIRatioHeader;
    DDIRatioTableContent];

testNaN=cellfun(@(x) any(isnan(x)),DDIRatioTable);
DDIRatioTable(testNaN)={'-'};

disp(DDIRatioTable);

% Get the DDI Ratio Qualification
for k=1:length(PKParameter)
    DDIRatioQualiHeader = {PKParameter{k}, 'Number', 'Ratio [%]'};
    DDIRatioQuali_1st_Column = {'Points total'; 'Points within Guest et al.'; 'Points within 2-fold'};
    DDIRatioQuali_Other_Columns = num2cell([QualiMeasure(k).PointsTotal NaN; ...
        QualiMeasure(k).PointsGuest 100.*QualiMeasure(k).PointsGuest./QualiMeasure(k).PointsTotal; ...
        QualiMeasure(k).Points2fold 100.*QualiMeasure(k).Points2fold./QualiMeasure(k).PointsTotal]);
    DDIRatioQuali_Other_Columns{1,2}='-';
    
    
    DDIRatioQuali(k).Output = [DDIRatioQualiHeader  ; ...
        DDIRatioQuali_1st_Column DDIRatioQuali_Other_Columns];
    
    fprintf('Qualification Measures for PK parameter %s \n', PKParameter{k});
    disp(DDIRatioQuali(k).Output);
end


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

symUpperLimit = symRobs.*Limit;
symLowerLimit = symRobs./Limit;

UpperLimit(isnan(Robs))=NaN;
LowerLimit(isnan(Robs))=NaN;

UpperLimit(Robs>=1)=symUpperLimit(Robs>=1);
LowerLimit(Robs>=1)=symLowerLimit(Robs>=1);
UpperLimit(Robs<1)=1./symLowerLimit(Robs<1);
LowerLimit(Robs<1)=1./symUpperLimit(Robs<1);