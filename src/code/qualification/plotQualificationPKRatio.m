function [fig_handle, PKRatioTable, PKRatioQuali, GMFE] = plotQualificationPKRatio(WSettings,plotIndex,PKParameter,PKRatioGroups,ObservedDataSets, SimulationMappings, AxesOptions, PlotSettings, REInputPath)
%PLOTQUALIFICATIONPKRATIO Plots PK Ratios from Configuration Plan
%
% [fig_handle, DDIRatioTable, DDIRatioQuali] = plotQualificationDDIRatio(WSettings,plotIndex,
%   PKParameter,DDIRatioGroups,ObservedDataSets, SimulationMappings, AxesOptions, PlotSettings, REInputPath)
%
% Inputs:
%   WSettings (structure)    definition of properties used in all
%                   workflow functions see GETDEFAULTWORKFLOWSETTINGS
%   plotIndex (integer) index of plot
%   PKParameter (cells) name of the PK parameter to be evaluated
%   PKRatioPlot (structure) PK Ratio plot information
%   ObservedDataSets (structure) Observed data
%   SimulationMappings (structure) Map simulation results to project
%   AxesOptions (structure) to set axes options
%   PlotSettings (structure) to set plot options
%   REInputPath (string) path of RE input files and folders
% Output
%   fig_handle (handle) handle of output figures
%   PKRatioTable (cells) Table of DDI Ratio information
%   GMFE (cells)  Geometric Mean Fold Error
%

% Open Systems Pharmacology Suite;  http://open-systems-pharmacology.org

%---------------------------------------------

% Loop on the Ratio Groups to be plotted by PK Ratio plot
for i=1:length(PKRatioGroups)
    
    PKRatios=PKRatioGroups(i).PKRatios;
    
    for j=1:length(PKRatios)
        
        % Get the observed data
        for k=1:length(ObservedDataSets)
            if strcmp(PKRatios(j).ObservedData, ObservedDataSets(k).Id)
                ObservedData = ObservedDataSets(k).y;
                break
            end
        end
        
        % Get the matching Record ID in the Table
        ID = (ObservedData.ID==PKRatios(j).ObservedDataRecordId);
        
        if max(ID)==0
            ME = MException('plotQualificationPKRatio:notFoundInPath', ...
                'In PK Ratio Plot %d, Group %d, Ratio %d, Study ID "%d" was not found in Observed Dataset', plotIndex, i, j, PKRatio.ObservedDataRecordId);
            throw(ME);
        end
        
        % Get the Study
        Result(i).Study(j,1) = ObservedData.Study(ID);
        
        % Get the requested PK Parameters, their unit and dimension
        for k=1:length(PKParameter)
            
            Result(i).obsPK(j, k) = table2array(ObservedData(ID, sprintf('%sAvg', PKParameter{k})));
            obsPKUnit = table2array(ObservedData(ID, sprintf('%sAvgUnit', PKParameter{k})));
            Result(i).obsPKUnit{j, k} = convertPKSimUnit(obsPKUnit);
            Result(i).obsPKDimension{j, k} = findDimensionfromUnit(Result(i).obsPKUnit{j, k});
            
            if isempty(Result(i).obsPKDimension{j, k})
                ME = MException('plotQualificationPKRatio:unknownUnit', ...
                    'In PK Ratio Plot %d, Group %d, Ratio %d, Dimension unknown for Unit "%s" in Observed Data Record ID "%d" \n', plotIndex, i, j, char(Result(i).obsPKUnit{j, k}), PKRatios(j).ObservedDataRecordId);
                throw(ME);
            end
            
        end
        
        % Load the mapped Time Profile Simulation Results
        [csvSimFile, xmlfile] = getSimFile(PKRatios(j), SimulationMappings, REInputPath);
        if isempty(csvSimFile)
            ME = MException('plotQualificationPKRatio:notFoundInPath', ...
                'In PK Ratio Plot %d, Group %d, Ratio %d, Project "%s" or Simulation "%s" was not found in SimulationMappings', plotIndex, i, j, PKRatios(j).Project, PKRatios(j).Simulation);
            throw(ME);
        end
        SimResult = loadSimResultcsv(csvSimFile, PKRatios(j).Simulation);
        
        % All the output are kept so far, may be removed if not necessary
        [AGE, BW, MW, drugmass, drugmassUnit] = getInfofromSimulation(xmlfile, PKRatios(j).Output);
        
        Result(i).AGE(j, 1)=AGE;
        Result(i).BW(j, 1)=BW;
        Result(i).MW(j, 1)=MW;
        Result(i).drugmass(j, 1)=drugmass(1);
        Result(i).drugmassUnit{j, 1}=drugmassUnit(1);
        
        % Get the right PK Output
        if isempty(SimResult.outputPathList)
            ME = MException('plotQualificationPKRatio:emptyOutputPathInSimulation', ...
                'In PK Ratio Plot %d, Group %d, Ratio %d, OutputPath is empty in Project "%s" or Simulation "%s"', plotIndex, i, j, PKRatios(j).Project, PKRatios(j).Simulation);
            throw(ME);
        end
        for indexPathList=1:length(SimResult.outputPathList)
            findPathOutput = strfind(SimResult.outputPathList{indexPathList}, PKRatios(j).Output);
            if ~isempty(findPathOutput)
                % Get Time and Concentration in PK Sim internal units
                % Concentration in µmol/l and
                % Time in min
                SimTime=SimResult.time;
                Pred=SimResult.y{indexPathList};
                
                break
            end
        end
        if isempty(findPathOutput)
            ME = MException('plotQualificationPKRatio:notFoundInPath', ...
                'In PK Ratio Plot %d, Group %d, Ratio %d, Output "%s" was not found in Project "%s" or Simulation "%s"', plotIndex, i, j, PKRatios(j).Output, PKRatios(j).Project, PKRatios(j).Simulation);
            throw(ME);
        end
        
        % Get PK parameters from the curves
        % For simulation
        allPKpred=getPKParametersForConcentration(SimTime, Pred, 'Dose', drugmass);
        
        for k=1:length(PKParameter)
            % Get the PK parameters requested in PK Parameters
            % Internal Units are assumed for PK parameters
            % according to Obs Unit
            if strcmp(Result(i).obsPKDimension(j, k), 'AUC (mass)')
                % Internal Unit for AUC is µmol*min/l and MW in g/mol
                AUCpred = getfield(allPKpred, 'AUC_last')*MW;
                AUCpredUnitFactor = getUnitFactor('µg*min/l', Result(i).obsPKUnit{j, k}, 'AUC (mass)');
                Result(i).predPK(j, k) = AUCpred.*AUCpredUnitFactor;
                
            elseif strcmp(Result(i).obsPKDimension(j, k), 'AUC (molar)')
                % Internal Unit for AUC is µmol*min/l
                AUCpred = getfield(allPKpred, 'AUC_last');
                AUCpredUnitFactor = getUnitFactor('µmol*min/l', Result(i).obsPKUnit{j, k}, 'AUC (molar)');
                Result(i).predPK(j, k) = AUCpred.*AUCpredUnitFactor;
                
            elseif strcmp(Result(i).obsPKDimension(j, k), 'Concentration')
                % Internal Unit for Cmax is µmol/l
                CMAXpred = getfield(allPKpred, 'cMax');
                CMAXpredUnitFactor = getUnitFactor('µmol/l', Result(i).obsPKUnit{j, k}, 'Concentration', 'MW', MW);
                Result(i).predPK(j, k) = CMAXpred.*CMAXpredUnitFactor;
                
            elseif strcmp(Result(i).obsPKDimension(j, k), 'Flow')
                % Internal Unit for CL is l/min
                CLpred = getfield(allPKpred, 'CL');
                CLpredUnitFactor = getUnitFactor('l/min', Result(i).obsPKUnit{j, k}, 'Flow');
                Result(i).predPK(j, k) = CLpred.*CLpredUnitFactor;
                
            elseif strcmp(Result(i).obsPKDimension(j, k), 'Flow per weight')
                % Internal Units for CL is l/min and for BW is kg
                CLpred = getfield(allPKpred, 'CL')/BW;
                CLpredUnitFactor = getUnitFactor('l/min/kg', Result(i).obsPKUnit{j, k}, 'Flow per weight');
                Result(i).predPK(j, k) = CLpred.*CLpredUnitFactor;
                
            else
                ME = MException('plotQualificationPKRatio:unknownDimension', ...
                    'In PK Ratio Plot %d, Group %d, Ratio %d, Observed Study ID "%d", PK Parameter "%s" \n Unknown dimension for observed unit "%s" ', plotIndex, i, j, PKRatios(j).ObservedDataRecordId, PKParameter{k}, char(Result(i).obsPKUnit{j, k}));
                throw(ME);
            end
            
            % Get the Ratio
            Result(i).RatioPK(j, k) = Result(i).predPK(j, k)./Result(i).obsPK(j, k);
            
        end
    end
end
%--------------------------------------
% Plot Section
% Check that X parameter is AGE
for k=1:length(AxesOptions)
    if strcmp(AxesOptions(k).Type, 'X')
        Xparam = AxesOptions(k).Dimension;
        % Field is saved as AGE, similar process can be
        % implmented for different X parameters
        if strcmpi(Xparam, 'AGE')
            Xparam = 'AGE';
        else
            ME = MException('plotQualificationPKRatio:XaxisDimension', ...
                'In PK Ratio Plot %d, X-axis does not correspond to Age', plotIndex);
            throw(ME);
        end
        break
    end
end

% Get all values in a single array
% To get limits of X and Y axes + GMFE and final tables
XcellValues = arrayfun(@(x)getfield(x, Xparam), Result, 'UniformOutput', false);
YcellValues = arrayfun(@(x)getfield(x, 'RatioPK'), Result, 'UniformOutput', false);
ObscellValues = arrayfun(@(x)getfield(x, 'obsPK'), Result, 'UniformOutput', false);
PredcellValues = arrayfun(@(x)getfield(x, 'RatioPK'), Result, 'UniformOutput', false);
StudycellValues = arrayfun(@(x)getfield(x, 'Study'), Result, 'UniformOutput', false);
AGEcellValues = arrayfun(@(x)getfield(x, 'AGE'), Result, 'UniformOutput', false);
BWcellValues = arrayfun(@(x)getfield(x, 'BW'), Result, 'UniformOutput', false);

XallValues=[]; YallValues= [];
ObsallValues=[]; PredallValues= [];
StudyallValues={}; AGEallValues= []; BWallValues=[];

for i=1:length(PKRatioGroups)
    XallValues = [XallValues; XcellValues{i}];
    YallValues = [YallValues; YcellValues{i}];
    ObsallValues = [ObsallValues; ObscellValues{i}];
    PredallValues = [PredallValues; PredcellValues{i}];
    StudyallValues = [StudyallValues; StudycellValues{i}];
    AGEallValues = [AGEallValues; AGEcellValues{i}];
    BWallValues = [BWallValues; BWcellValues{i}];
end
Xrange = [0.8*min(XallValues) 1.2*max(XallValues)];
Yrange = [1 1];

% One plot per PK parameter
for k=1:length(PKParameter)
    
    % Initialize legend Labels
    legendLabels={};
    
    % create figure for Obs vs Pred
    [ax, fig_handle(k).PKRatio] = getReportFigureQP(WSettings,1,1,[],PlotSettings);
    setFigureOptions(AxesOptions);
    
    % Get limits of X and Y axes
    YrangePerParameter = [0.8*min(0.5, min(YallValues(:,k))) 1.2*max(2, max(YallValues(:,k)))];
    
    % Plot lines of specific PK Ratios
    plot(Xrange, Yrange, '-k', 'LineWidth', 1, 'HandleVisibility','off');
    plot(Xrange, Yrange/2, '--r', 'LineWidth', 1, 'HandleVisibility','off');
    plot(Xrange, Yrange*2, '--r', 'LineWidth', 1, 'HandleVisibility','off');
    plot(Xrange, Yrange/1.5, '--b', 'LineWidth', 1, 'HandleVisibility','off');
    plot(Xrange, Yrange*1.5, '--b', 'LineWidth', 1, 'HandleVisibility','off');
    
    ylabel(sprintf('Predicted %s / Observed %s', PKParameter{k}, PKParameter{k}));
    
    axis([Xrange YrangePerParameter]);
    
    % Plot PK Ratios per Group
    for i=1:length(PKRatioGroups)
        
        Xvalues = getfield(Result(i), Xparam);
        
        pp=plot(Xvalues, Result(i).RatioPK(:, k), 'o', 'Linewidth',1);
        setCurveOptions(pp, PKRatioGroups(i));
        
        if isfield(PKRatioGroups(i), 'Caption')
            legendLabels=[legendLabels PKRatioGroups(i).Caption];
        end
    end
    
    if ~isempty(legendLabels)
        lgd = legend(legendLabels, 'Location', 'northoutside');
        reshapeQualificationLegend(lgd);
    else
        legend('off');
    end
    
end

% --------------------------------------
% Table and Qualification Section

% Calculation of GMFE
GMFE = 10.^(sum(abs(log10(YallValues)))./size(YallValues, 1));

for k=1:length(PKParameter)
    fprintf('%s: GMFE = %f \n', PKParameter{k}, GMFE(k));
end

% Get the PK Ratio Table
PKRatioHeader={};
PKRatioResults=[];

for k=1:length(PKParameter)
    PKRatioHeader = [PKRatioHeader {sprintf('Predicted %s [%s]', PKParameter{k}, Result(1).obsPKUnit{1, k}), ...
        sprintf('Observed %s [%s]', PKParameter{k}, Result(1).obsPKUnit{1, k}), sprintf('Pred/Obs %s Ratio', PKParameter{k})}];
    
    PKRatioResults = [PKRatioResults PredallValues(:,k) ObsallValues(:,k) YallValues(:,k)];
end

% Definition of the PKRatio table
PKRatioHeader = [{'Study ID', 'Age [y]', 'BodyWeight [kg]'}, PKRatioHeader];

PKRatioTable = [PKRatioHeader ; ...
    StudyallValues  num2cell([AGEallValues BWallValues PKRatioResults])];

disp(PKRatioTable);

% Get the PK Ratio Qualification
for k=1:length(PKParameter)
    % Update the number of points for each group
    conditionPoints = ~isnan(YallValues(:,k));
    QualiMeasure(k).PointsTotal = length(YallValues(conditionPoints,k));
    condition15fold = (YallValues(:,k) >= 1/1.5 & YallValues(:,k) <= 1.5);
    QualiMeasure(k).Points15fold= length(YallValues(condition15fold,k));
    condition2fold = (YallValues(:,k) >= 0.5 & YallValues(:,k) <= 2);
    QualiMeasure(k).Points2fold= length(YallValues(condition2fold,k));
    
    PKRatioQualiHeader = {PKParameter{k}, 'Number', 'Ratio [%]'};
    PKRatioQuali_1st_Column = {'Points total'; 'Points within 1.5 fold'; 'Points within 2-fold'};
    PKRatioQuali_Other_Columns = num2cell([QualiMeasure(k).PointsTotal NaN; ...
        QualiMeasure(k).Points15fold 100.*QualiMeasure(k).Points15fold./QualiMeasure(k).PointsTotal; ...
        QualiMeasure(k).Points2fold 100.*QualiMeasure(k).Points2fold./QualiMeasure(k).PointsTotal]);
    PKRatioQuali_Other_Columns{1,2}='-';
    
    
    PKRatioQuali(k).Output = [PKRatioQualiHeader  ; ...
        PKRatioQuali_1st_Column PKRatioQuali_Other_Columns];
    
    fprintf('Qualification Measures for PK parameter %s \n', PKParameter{k});
    disp(PKRatioQuali(k).Output);
    
end

function [AGE, BW, MW, drugmass, drugmassUnit] = getInfofromSimulation(xmlfile, Output)

initSimulation(xmlfile,'none');

MW = getMolecularWeightForPath(Output);

drugmass = getParameter('*Application_*|ProtocolSchemaItem|DrugMass',1,'parametertype','readonly');
drugmassUnit = getParameter('*Application_*|ProtocolSchemaItem|DrugMass',1,'parametertype','readonly', 'property', 'Unit');

AGE = getParameter('*|Organism|Age',1,'parametertype','readonly');
BW = getParameter('*|Organism|Weight',1,'parametertype','readonly');

function UnitOut = convertPKSimUnit(UnitIn)

UnitOut = string(UnitIn);
UnitOut = replace(UnitOut,'L','l');

%UnitOut(1)='µ';