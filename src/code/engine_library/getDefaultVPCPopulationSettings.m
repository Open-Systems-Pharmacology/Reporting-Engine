function VPC = getDefaultVPCPopulationSettings(WSettings,PopRunSet,flag) %#ok<INUSL>
%GETDEFAULTVPCPOPULATIONSETTINGS get WSettings for visual predicitve check
%
%  VPC = getDefaultVPCPopulationSettings(WSettings,PopRunSet,flag)
%
% Inputs:
%       WSettings (structure)    definition of properties used in all
%                   workflow functions see GETDEFAULTWORKFLOWSETTINGS
%       PopRunSet (structure)   list of population simulations see GENERATEWORKFLOWINPUTFORPOPULATIONSIMULATION
%       flag (string)  defines type of population, default = 'paralellComparison'
%                               other possibilities are 'pediatric';

% Open Systems Pharmacology Suite;  http://open-systems-pharmacology.org


% add function handle for figure and legend text
VPC.textFunctionHandle = @textVPCPopulation;

% add default structure for demographi plots
VPC.PhysProperties = addVPCPhysProperties(PopRunSet,flag);

% add default structure for timeprofile plots
VPC.Timeprofile = addVPCTimeprofile(PopRunSet,flag);

% add default structure for PK-Parameter plots
VPC.PKParameter = addVPCPKParameter(PopRunSet,flag);

return


function timeprofile = addVPCTimeprofile(PopRunSet,~)

% list indices of popSet for that the VPC should be done
jj_ref = [PopRunSet.isReference];
timeprofile.ixOfPopRunSets =  find(~jj_ref);

% name of Reference Population
timeprofile.ixPopRunSetRef = find(jj_ref,1);

% display Unit for tim
timeprofile.timeDisplayUnit = 'h';

% zoom int time axes
timeprofile.timelimit = [];

return

function PhysProperties = addVPCPhysProperties(PopRunSet,flag)

switch flag
    
    case 'parallelComparison'
        
        PhysProperties.name = '';
        PhysProperties.reportName = '';
        
        % list demographic properties which should be plotted
        % pathID, reortname, display unit,categorial text
        PhysProperties.yList = {'Organism|Age', 'Age', 'year(s)',{};...
            'Organism|Weight', 'Body weight', 'kg',{};...
            'Organism|Height', 'Height', 'cm',{};...
            'Organism|BMI', 'BMI', 'kg/m²',{};...
            'Gender', '', '',{'Male','Female'}};
        
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
            'Gender', '', '',{'Male','Female'};...
            '<addOntogenyFactor>','','',{}};
        
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
        PhysProperties.ixPopRunSetRef = ix(1);
        
end



return



function PKParameter = addVPCPKParameter(PopRunSet,flag)

switch flag
    
    case 'parallelComparison'
        
        PKParameter.name = '';
        PKParameter.reportName = '';
                
        % list propety which shall be displayed on the x-axes, if empty a
        % histogramm will be created
        PKParameter.xList = {};
        
        % list indices of popSet for that the VPC should be done
        [~,ix] = unique({PopRunSet.popcsv},'stable');
        PKParameter.ixOfPopRunSets =  ix;
                
        % name of Reference Population
        PKParameter.ixPopRunSetRef = [];
                
    case 'pediatric'

        PKParameter.name = 'pediatricPopulation';
        PKParameter.reportName = 'virtual pediatric population';

        % list propety which shall be displayed on the x-axes, if empty a
        % histogramm will be created
        PKParameter.xList = {'Organism|Age', 'age', 'year(s)';...
            'Organism|Weight', 'weigth', 'kg'};
        
        % list popset for that the VPC should be done
        jj_ref = [PopRunSet.isReference];
        [~,ix_unique] = unique({PopRunSet.popcsv},'stable');
        ix = intersect(find(~jj_ref),ix_unique);
        PKParameter.ixOfPopRunSets = ix ;
                
        % name of Reference Population
        ix = intersect(find(jj_ref),ix_unique);
        PKParameter.ixPopRunSetRef = ix(1);
        
end



return
