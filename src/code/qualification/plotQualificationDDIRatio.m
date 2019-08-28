function [fig_handle, DDIRatioTable, DDIRatioQuali, GMFE] = plotQualificationDDIRatio(WSettings,plotIndex,PKParameter,DDIRatioGroups,ObservedDataSets, SimulationMappings, AxesOptions, PlotSettings, REInputPath)
%PLOTQUALIFICATIONDDIRATIO Plots DDI Ratios from Configuration Plan
%
% [fig_handle, DDIRatioTable, DDIRatioQuali] = plotQualificationDDIRatio(WSettings,plotIndex,
%   PKParameter,DDIRatioGroups,ObservedDataSets, SimulationMappings, AxesOptions, PlotSettings, REInputPath)
%
% Inputs:
%   WSettings (structure)    definition of properties used in all
%                   workflow functions see GETDEFAULTWORKFLOWSETTINGS
%   plotIndex (integer) index of plot
%   PKParameter (cells) name of the PK parameter to be evaluated
%   DDIRatioGroups (structure) DDI Ratio plot information
%   ObservedDataSets (structure) Observed data
%   SimulationMappings (structure) Map simulation results to project
%   AxesOptions (structure) to set axes options
%   PlotSettings (structure) to set plot options
%   REInputPath (string) path of RE input files and folders
% Output
%   fig_handle (handle) handle of output figures
%   DDIRatioTable (cells) Table of DDI Ratio information
%   DDIRatioQuali (cells)  Table of DDI Ratio qualification
%

% Open Systems Pharmacology Suite;  http://open-systems-pharmacology.org

%---------------------------------------------
% Create figure with first setting from WSettings using getReportFigure
% To be updated using the Configuration plan Settings as optional arguments

% Initialize content of DDI table, and axis
DDIRatioTableContent={};

for k=1:length(PKParameter)
    axisObsVsPred(k).min=NaN;
    axisObsVsPred(k).max=NaN;
    axisResVsObs(k).min=NaN;
    axisResVsObs(k).max=NaN;
end

% Loop on the Ratio Groups to be plotted by DDI Ratio plot
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
        
        % Get the matching Record ID in the Table
        ID = (ObservedData.ID==DDIRatios(j).ObservedDataRecordId);
        
        if max(ID)==0
            ME = MException('plotQualificationDDIRatio:notFoundInPath', ...
                'In DDI Ratio Plot %d, Group %d, Ratio %d, Study ID "%d" was not found in Observed Dataset', plotIndex, i, j, DDIRatios(j).ObservedDataRecordId);
            throw(ME);
        end
        
        % Get the Study
        Result(i).Study(j) = ObservedData.StudyID(ID);
        
        % Load the mapped Time Profile Simulation Results
        [csvSimFileControl, xmlfileControl] = getSimFile(DDIRatios(j).SimulationControl, SimulationMappings, REInputPath);
        if isempty(csvSimFileControl)
            ME = MException('plotQualificationDDIRatio:notFoundInPath', ...
                'In DDI Ratio Plot %d, Group %d, Ratio %d, Project "%s" or Simulation "%s" for Control was not found in SimulationMappings', plotIndex, i, j, DDIRatios(j).SimulationControl.Project, DDIRatios(j).SimulationControl.Simulation);
            throw(ME);
        end
        SimResultControl = loadSimResultcsv(csvSimFileControl, DDIRatios(j).SimulationControl.Simulation);
        
        if isempty(SimResultControl.outputPathList)
            ME = MException('plotQualificationDDIRatio:emptyOutputPathInSimulation', ...
                'In DDI Ratio Plot %d, Group %d, Ratio %d, OutputPath is empty in Project "%s" Simulation "%s" for Control', plotIndex, i, j, DDIRatios(j).SimulationControl.Project, DDIRatios(j).SimulationControl.Simulation);
            throw(ME);
        end
        
        % All the output are kept so far, may be removed if not necessary
        [AGEControl, BWControl, MWControl, drugmassControl] = getInfofromSimulation(xmlfileControl, DDIRatios(j).Output);
        
        Result(i).AGEControl(j)=AGEControl;
        Result(i).BWControl(j)=BWControl;
        Result(i).MWControl(j)=MWControl;
        Result(i).drugmassControl(j)=drugmassControl(1);
        
        [csvSimFileDDI, xmlfileDDI] = getSimFile(DDIRatios(j).SimulationDDI, SimulationMappings, REInputPath);
        if isempty(csvSimFileDDI)
            ME = MException('plotQualificationDDIRatio:notFoundInPath', ...
                'In DDI Ratio Plot %d, Group %d, Ratio %d, Project "%s" or Simulation "%s" for DDI was not found in SimulationMappings', plotIndex, i, j, DDIRatios(j).SimulationDDI.Project, DDIRatios(j).SimulationDDI.Simulation);
            throw(ME);
        end
        SimResultDDI = loadSimResultcsv(csvSimFileDDI, DDIRatios(j).SimulationDDI.Simulation);
        
        if isempty(SimResultDDI.outputPathList)
            ME = MException('plotQualificationDDIRatio:emptyOutputPathInSimulation', ...
                'In DDI Ratio Plot %d, Group %d, Ratio %d, OutputPath is empty in Project "%s" Simulation "%s" for DDI Treatment', plotIndex, i, j, DDIRatios(j).SimulationDDI.Project, DDIRatios(j).SimulationDDI.Simulation);
            throw(ME);
        end
        
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
                
                % If Start Time and End Time are empty, use min and max values of simulation
                DDIRatios(j).SimulationControl.StartTime(isempty(DDIRatios(j).SimulationControl.StartTime)) = min(ControlTime.*Xfactor);
                DDIRatios(j).SimulationControl.EndTime(isempty(DDIRatios(j).SimulationControl.EndTime)) = max(ControlTime.*Xfactor);
                
                % Check if end time was read as a string for Inf value
                if strcmpi(DDIRatios(j).SimulationControl.EndTime, 'Inf')
                    DDIRatios(j).SimulationControl.EndTime = Inf;
                end
                
                SimTime = (ControlTime.*Xfactor >= DDIRatios(j).SimulationControl.StartTime &  ControlTime.*Xfactor <= DDIRatios(j).SimulationControl.EndTime);
                ControlTime = ControlTime(SimTime).*Xfactor;
                Controlpred = Controlpred(SimTime);
                
                % Check that start and end time lead to a non-empty simulation time range
                if isempty(ControlTime)
                    ME = MException('plotQualificationDDIRatio:notFoundInPath', ...
                        ['In DDI Ratio Plot %d, Group %d, Ratio %d, Output "%s" from Control Project "%s" and Simulation "%s" \n',...
                        'Requested time range leads to empty simulation'], plotIndex, i, j, DDIRatios(j).Output, DDIRatios(j).SimulationControl.Project, DDIRatios(j).SimulationControl.Simulation);
                    throw(ME);
                end
                
                break
            end
        end
        if isempty(findPathOutput)
            ME = MException('plotQualificationDDIRatio:notFoundInPath', ...
                'In DDI Ratio Plot %d, Group %d, Ratio %d, Output "%s" for Control was not found in Project "%s" or Simulation "%s"', plotIndex, i, j, DDIRatios(j).Output, DDIRatios(j).SimulationDDI.Project, DDIRatios(j).SimulationDDI.Simulation);
            throw(ME);
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
                
                % If Start Time and End Time are empty, use min and max values of simulation
                DDIRatios(j).SimulationDDI.StartTime(isempty(DDIRatios(j).SimulationDDI.StartTime)) = min(DDITime.*Xfactor);
                DDIRatios(j).SimulationDDI.EndTime(isempty(DDIRatios(j).SimulationDDI.EndTime)) = max(DDITime.*Xfactor);
                
                % Check if end time was read as a string for Inf value
                if strcmpi(DDIRatios(j).SimulationDDI.EndTime, 'Inf')
                    DDIRatios(j).SimulationDDI.EndTime = Inf;
                end
                
                SimTime = (DDITime.*Xfactor >= DDIRatios(j).SimulationDDI.StartTime &  DDITime.*Xfactor <= DDIRatios(j).SimulationDDI.EndTime);
                DDITime = DDITime(SimTime).*Xfactor;
                DDIpred = DDIpred(SimTime);
                
                % Check that start and end time lead to a non-empty simulation time range
                if isempty(DDITime)
                    ME = MException('plotQualificationDDIRatio:notFoundInPath', ...
                        ['In DDI Ratio Plot %d, Group %d, Ratio %d, Output "%s" from DDI Project "%s" and Simulation "%s" \n',...
                        'Requested time range leads to empty simulation'], plotIndex, i, j, DDIRatios(j).Output, DDIRatios(j).SimulationDDI.Project, DDIRatios(j).SimulationDDI.Simulation);
                    throw(ME);
                end
                
                break
            end
        end
        if isempty(findPathOutput)
            ME = MException('plotQualificationDDIRatio:notFoundInPath', ...
                'In DDI Ratio Plot %d, Group %d, Ratio %d, Output "%s" for Control was not found in Project "%s" or Simulation "%s"', plotIndex, i, j, DDIRatios(j).Output, DDIRatios(j).SimulationDDI.Project, DDIRatios(j).SimulationDDI.Simulation);
            throw(ME);
        end
        
        % Get the PK parameters out of the simulation
        allPKpredControl=getPKParametersForConcentration(ControlTime,Controlpred,'Dose', drugmassControl);
        
        % Get the PK parameters out of the simulation
        allPKpredDDI=getPKParametersForConcentration(DDITime,DDIpred,'Dose', drugmassDDI);
        
        for k=1:length(PKParameter)
            % Get the PK parameters requested in PKParameter
            if strcmpi(PKParameter{k}, 'AUC')
                if isinf(DDIRatios(j).SimulationDDI.EndTime)
                    PKpredField(k).DDI = 'AUC_inf';
                else
                    PKpredField(k).DDI = 'AUC_last';
                end
                if isinf(DDIRatios(j).SimulationControl.EndTime)
                    PKpredField(k).Control = 'AUC_inf';
                else
                    PKpredField(k).Control = 'AUC_last';
                end
            elseif strcmpi(PKParameter{k}, 'CMAX')
                PKpredField(k).DDI = 'cMax';
                PKpredField(k).Control = 'cMax';
            else
                PKpredField(k).DDI = PKParameter{k};
                PKpredField(k).Control = PKParameter{k};
            end
            if isfield(allPKpredControl, PKpredField(k).Control) && isfield(allPKpredDDI, PKpredField(k).DDI)
                
                % Get the observation Ratios
                Observations(i).RatioPK(k,j) = table2array(ObservedData(ObservedData.ID==DDIRatios(j).ObservedDataRecordId,...
                    (strcmpi(ObservedData.Properties.VariableNames, [PKParameter{k} 'RAvg']))));
                
                Result(i).DDIPK(k,j)=getfield(allPKpredDDI, PKpredField(k).DDI);
                
                Result(i).ControlPK(k,j)=getfield(allPKpredControl, PKpredField(k).Control);
                
                % Assumes currently that the output from PKSim has same unit for Control and DDI
                Result(i).RatioPK(k,j)=Result(i).DDIPK(k,j)./Result(i).ControlPK(k,j);
                
                % Set the axis window
                axisObsVsPred(k).min=nanmin(nanmin(Result(i).RatioPK(k,j), Observations(i).RatioPK(k,j)),axisObsVsPred(k).min);
                axisObsVsPred(k).max=nanmax(nanmax(Result(i).RatioPK(k,j), Observations(i).RatioPK(k,j)),axisObsVsPred(k).max);
                axisResVsObs(k).min=nanmin(Result(i).RatioPK(k,j)./Observations(i).RatioPK(k,j),axisResVsObs(k).min);
                axisResVsObs(k).max=nanmax(Result(i).RatioPK(k,j)./Observations(i).RatioPK(k,j),axisResVsObs(k).max);
                
            else
                ME = MException('DDIRatio:notFoundInField', ...
                    'In DDI Ratio Plot %d, Requested PK Parameter %s not found in parameters extracted using getPKParametersForConcentration', plotIndex, PKParameter{k});
                throw(ME);
            end
        end
        % Build the DDI Ratio Table:
        Perpetrator = table2cell(ObservedData(ID,{'Perpetrator', 'Dose', 'DoseUnit', 'RoutePerpetrator'}));
        Victim = table2cell(ObservedData(ID, {'Victim', 'RouteVictim'}));
        Reference = table2cell(ObservedData(ID, {'StudyID'}));
        
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

%-------------------------------------------------------------
% Plot Section

% Initialize Quaklification Measure and
% Guest equation for each PK Parameter

for k=1:length(PKParameter)
    
    % GuestRatio vectors for plotting Guest et al. equation
    GuestRatio(k).x=10.^(log10(0.8*axisObsVsPred(k).min):0.01:log10(1.2*axisObsVsPred(k).max));
    [GuestRatio(k).yup, GuestRatio(k).ylo] = DDIRatioGuestEquation(GuestRatio(k).x);
    
    % Initialize the Qualification Measures
    QualiMeasure(k).PointsTotal = 0;
    QualiMeasure(k).PointsGuest= 0;
    QualiMeasure(k).Points2fold= 0;
end

for k=1:length(PKParameter)
    % Initialize legend labels
    leg_labels={};
    
    % create figure for Obs vs Pred
    [ax, fig_handle(k).predictedVsObserved] = getReportFigureQP(WSettings,1,1,[],PlotSettings);
    setFigureOptions(AxesOptions.DDIRatioPlotsPredictedVsObserved);
    plot(GuestRatio(k).x, GuestRatio(k).x, '--k', 'Linewidth', 1, 'HandleVisibility','off');
    plot(GuestRatio(k).x, GuestRatio(k).x/2, '--k', 'Linewidth', 1, 'HandleVisibility','off');
    plot(GuestRatio(k).x, GuestRatio(k).x*2, '--k', 'Linewidth', 1, 'HandleVisibility','off');
    plot(GuestRatio(k).x, GuestRatio(k).yup, '--k', 'Linewidth', 1, 'HandleVisibility','off');
    plot(GuestRatio(k).x, GuestRatio(k).ylo, '--k', 'Linewidth', 1, 'HandleVisibility','off');
    xlabel(sprintf('Observed %s Ratio', PKParameter{k})); ylabel(sprintf('Predicted %s Ratio', PKParameter{k}));
    axis([0.8*axisObsVsPred(k).min 1.2*axisObsVsPred(k).max 0.8*axisObsVsPred(k).min 1.2*axisObsVsPred(k).max]);
    
    % create figure for Residuals
    [ax, fig_handle(k).residualsVsObserved] = getReportFigureQP(WSettings,1,1,[],PlotSettings);
    setFigureOptions(AxesOptions.DDIRatioPlotsResidualsVsObserved);
    plot(GuestRatio(k).x, ones(size(GuestRatio(k).x)), '--k', 'Linewidth', 1, 'HandleVisibility','off');
    plot(GuestRatio(k).x, ones(size(GuestRatio(k).x))/2, '--k', 'Linewidth', 1, 'HandleVisibility','off');
    plot(GuestRatio(k).x, ones(size(GuestRatio(k).x))*2, '--k', 'Linewidth', 1, 'HandleVisibility','off');
    plot(GuestRatio(k).x, GuestRatio(k).yup./GuestRatio(k).x, '--k', 'Linewidth', 1, 'HandleVisibility','off');
    plot(GuestRatio(k).x, GuestRatio(k).ylo./GuestRatio(k).x, '--k', 'Linewidth', 1, 'HandleVisibility','off');
    xlabel(sprintf('Observed %s Ratio', PKParameter{k})); ylabel(sprintf('Predicted %s Ratio / Observed %s Ratio', PKParameter{k}, PKParameter{k}));
    axis([0.8*axisObsVsPred(k).min 1.2*axisObsVsPred(k).max 0.8*axisResVsObs(k).min 1.2*axisResVsObs(k).max]);
    
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
        set(0, 'CurrentFigure', fig_handle(k).predictedVsObserved);
        pp=plot(Observations(i).RatioPK(k,:), Result(i).RatioPK(k,:), 'o', 'Linewidth',1);
        setCurveOptions(pp, DDIRatioGroups(i));
        
        set(0, 'CurrentFigure', fig_handle(k).residualsVsObserved);
        pp=plot(Observations(i).RatioPK(k,:), Result(i).RatioPK(k,:)./Observations(i).RatioPK(k,:), 'o', 'Linewidth',1);
        setCurveOptions(pp, DDIRatioGroups(i));
        
        if isfield(DDIRatioGroups(i), 'Caption')
            leg_labels=[leg_labels DDIRatioGroups(i).Caption];
        end
    end
    set(0, 'CurrentFigure', fig_handle(k).predictedVsObserved);
    legend(leg_labels, 'Location', 'northoutside');
        
    set(0, 'CurrentFigure', fig_handle(k).residualsVsObserved);
    legend(leg_labels, 'Location', 'northoutside');
    
end

% Calculation of GMFE
RatioPK=[];
for i=1:length(DDIRatioGroups)
    RatioPK = [RatioPK Result(i).RatioPK./Observations(i).RatioPK];
end
RatioPK = reshape(RatioPK', [], length(PKParameter));
for k=1:length(PKParameter)
    GMFE(k) = 10.^(sum(abs(log10(RatioPK(~isnan(RatioPK(:,k)),k))))./length(RatioPK(~isnan(RatioPK(:,k)),k)));
    fprintf('%s: GMFE = %f \n', PKParameter{k}, GMFE(k));
end 

%-------------------------------------------------------------
% Tables Section

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
