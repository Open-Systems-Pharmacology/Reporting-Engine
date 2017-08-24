function VPC = getDefaultVPCSettings(WSettings,SimulationSet) 
%GETDEFAULTVPCPOPULATIONSETTINGS get Settings for visual predicitve check
%
%  VPC = getDefaultVPCPopulationSettings(WSettings,PopRunSet,flag)
%
% Inputs:
%       SimulationSet (structure)   list of  simulations see GENERATEWORKFLOWINPUTFORPOPULATIONSIMULATION
%       workflowType (string)  defines type of population, default = 'paralellComparison'
%                               other possibilities are 'pediatric', ratioComparison;

% Open Systems Pharmacology Suite;  http://open-systems-pharmacology.org

% add function handle for figure and legend text
VPC.textFunctionHandle = @textVPC;


if strcmp(WSettings.workflowType,'popModel')
    
    % add default structure for demographi plots
    VPC.PhysProperties = addVPCPhysProperties(SimulationSet,WSettings.workflowMode);
end

% add default structure for timeprofile plots
VPC.Timeprofile = addVPCTimeprofile(SimulationSet,WSettings.workflowType);

if strcmp(WSettings.workflowType,'popModel')
    % add default structure for PK-Parameter plots
    VPC.PKParameter = addVPCPKParameter(SimulationSet,WSettings.workflowMode);
end

% get List of optimized Parameters
VPC.optimizedParameters = {};

return


function Timeprofile = addVPCTimeprofile(SimulationSet,workflowType)

switch workflowType
    case 'popModel'
        % list indices of popSet for that the VPC should be done
        jj_ref = [SimulationSet.isReference];
        Timeprofile.ixOfRunSets =  find(~jj_ref);
        
        % name of Reference Population
        Timeprofile.ixRunSetRef = find(jj_ref,1);
    case 'meanModel'
        % list indices of MeanModelSet for that the VPC should be done
        Timeprofile.ixOfRunSets = 1:length(SimulationSet);
        
        % name of Reference Population
        Timeprofile.ixRunSetRef = [];
end

% display Unit for tim
Timeprofile.timeDisplayUnit = 'h';

% zoom int time axes
Timeprofile.timelimit = [];

% set scale of y axes
Timeprofile.yScale = {'lin','log'};

% get Types of Plot
Timeprofile.plotTypes = {'yVsTime','resVsTime','resVsY','predVsObs','histRes','qqPlotRes'};



return

function PhysProperties = addVPCPhysProperties(PopRunSet,flag)

switch flag
    
    case {'parallelComparison','ratioComparison'}
        
        PhysProperties.name = '';
        PhysProperties.reportName = '';
        
        % list demographic properties which should be plotted
        % pathID, reortname, display unit,categorial text
        PhysProperties.yList = {'Organism|Age', 'Age', 'year(s)',{};...
            'Organism|Weight', 'Body weight', 'kg',{};...
            'Organism|Height', 'Height', 'cm',{};...
            'Organism|BMI', 'BMI', 'kg/m²',{};...
            'Organism|BSA', 'Body surface area', 'm²',{};...
            'Gender', 'Gender', '',{'Male','Female'}};
        
        % list propety which shall be displayed on the x-axes, if empty a
        % histogramm will be created
        PhysProperties.xList = {};
        
        % list indices of popSet for that the VPC should be done
        [~,ix] = unique({PopRunSet.popcsv},'stable');
        PhysProperties.ixOfPopRunSets =  ix;
        
        % how to handle PopRunSets:
        % merged populations are merged,
        % serial, for each PopRunSet an extra plot is generated
        PhysProperties.popRunSetMerger =  'serial';
        
        % name of Reference Population
        PhysProperties.ixPopRunSetRef = [];
        
    case 'pediatric'

        PhysProperties.name = 'pediatricPopulation';
        PhysProperties.reportName = 'virtual pediatric population';

        
        % list demographic properties which should be plotted
        % pathID, reortname, display unit
        PhysProperties.yList = {...
            'Organism|Weight', 'Body weight', 'kg',{};...
            'Organism|Height', 'Height', 'cm',{};...
            'Organism|BMI', 'BMI', 'kg/m²',{};...
            'Organism|BSA', 'Body surface area', 'm²',{};...
           'Gender', '', '',{'Male','Female'};...
            '<addOntogenyFactor>','','',{}}; % if the code finds the key word <addOntogenyFactor> 
        % in this cell array then all ontogeny factors included in the xml file are added to this structure
        % automatically during initialization of the workflow.
        
        % list propety which shall be displayed on the x-axes, if empty a
        % histogramm will be created
        PhysProperties.xList = {'Organism|Age', 'age', 'year(s)'};
        
        % list popset for that the VPC should be done
        jj_ref = [PopRunSet.isReference];
        [~,ix_unique] = unique({PopRunSet.popcsv},'stable');
        ix = intersect(find(~jj_ref),ix_unique);
        PhysProperties.ixOfPopRunSets = ix ;
        
        % how to handle PopRunSets:
        % merged populations are merged,
        % serial, for each PopRunSet an extra plot is generated
        PhysProperties.popRunSetMerger =  'merged';
        
        % name of Reference Population
        ix = intersect(find(jj_ref),ix_unique);
        if ~isempty(ix)
            PhysProperties.ixPopRunSetRef = ix(1);
        else
            PhysProperties.ixPopRunSetRef = [];
        end
    otherwise
        error('unknown flag')
end



return



function PKParameter = addVPCPKParameter(PopRunSet,flag)

% set Properties to plot
switch flag
    
    case {'parallelComparison','ratioComparison'}
        
        PKParameter.name = '';
        PKParameter.reportName = '';
                
        % list propety which shall be displayed on the x-axes, if empty a
        % histogramm will be created
        PKParameter.xList = {};
        
    case 'pediatric'

        PKParameter.name = 'pediatricPopulation';
        PKParameter.reportName = 'virtual pediatric population';

        % list propety which shall be displayed on the x-axes, if empty a
        % histogramm will be created
        PKParameter.xList = {'Organism|Age', 'age', 'year(s)';...
            'Organism|Weight', 'weigth', 'kg'};
    otherwise
                error('unknown flag')

end


% set indices
switch flag
    
    case {'parallelComparison'}
        % list indices of popSet for that the VPC should be done
        [~,ix] = unique({PopRunSet.popcsv},'stable');
        PKParameter.ixOfPopRunSets =  ix;
                
        % name of Reference Population
        PKParameter.ixPopRunSetRef = [];
                
    case {'pediatric','ratioComparison'}

        
        % list popset for that the VPC should be done
        jj_ref = [PopRunSet.isReference];
        [~,ix_unique] = unique({PopRunSet.popcsv},'stable');
        ix = intersect(find(~jj_ref),ix_unique);
        PKParameter.ixOfPopRunSets = ix ;
                
        % name of Reference Population
        ix = intersect(find(jj_ref),ix_unique);
        PKParameter.ixPopRunSetRef = ix(1);
    otherwise
        error('unknown flag')
  
end



return
