function VPC = getDefaultVPCSettings(WSettings,SimulationSet) 
%GETDEFAULTVPCSETTINGS get Settings for visual predicitve check
%
%  VPC = getDefaultVPCSettings(WSettings,SimulationSet) 
%
% Inputs:
%       WSettings (structure)    definition of properties used in all
%                   workflow functions see GETDEFAULTWORKFLOWSETTINGS
%       SimulationSet (structure)   list of  simulations 

% Open Systems Pharmacology Suite;  http://open-systems-pharmacology.org

% add function handle for figure and legend text
VPC.textFunctionHandle = @textVPC;


if strcmp(WSettings.workflowType,'popModel')
    
    % add default structure for demographi plots
    VPC.PhysProperties = addVPCPhysProperties(SimulationSet,WSettings.workflowMode);
end

% add default structure for timeprofile plots
VPC.Timeprofile = addVPCTimeprofile(SimulationSet,WSettings.workflowType,WSettings.workflowMode);

if strcmp(WSettings.workflowType,'popModel')
    % add default structure for PK-Parameter plots
    VPC.PKParameter = addVPCPKParameter(SimulationSet,WSettings.workflowMode);
end

% get List of optimized Parameters
VPC.optimizedParameters = {};

return


function Timeprofile = addVPCTimeprofile(SimulationSet,workflowType,workflowMode)

switch workflowType
    case 'popModel'
        % list indices of popSet for that the VPC should be done
        jj_ref = [SimulationSet.isReference];
        Timeprofile.ixOfRunSets =  find(~jj_ref);
        
        % name of Reference Population
        Timeprofile.ixRunSetRef = find(jj_ref,1);
        
        switch workflowMode
            case {'pediatric','ratioComparison'}
                Timeprofile.plotMeanModel = false;
            case 'parallelComparison'
                Timeprofile.plotMeanModel = true;
        end
        
    case 'meanModel'
        % list indices of MeanModelSet for that the VPC should be done
        Timeprofile.ixOfRunSets = 1:length(SimulationSet);
        
        % name of Reference Population
        Timeprofile.ixRunSetRef = [];
    otherwise 
        error('unknown workflowtype')
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
        % pathID, reportname (short), display unit,categorial text, report
        % name long (in Caption)
        PhysProperties.yList = {'Organism|Age', 'age', 'year(s)',{},'age';...
            'Organism|Weight', 'body weight', 'kg',{},'body weight';...
            'Organism|Height', 'body height', 'cm',{},'body height';...
            'Organism|BMI', 'BMI', 'kg/m²',{},'BMI';...
            'Organism|BSA', 'body surface area', 'm²',{},'body surface area';...
            'Gender', 'gender', '',{'Male','Female'},'gender'};
        
        % list propety which shall be displayed on the x-axes, if empty a
        % histogramm will be created
        PhysProperties.xList = {};
        
        % list indices of popSet for that the VPC should be done
        [~,ix] = unique(strcat({PopRunSet.popcsv},{PopRunSet.dataTpFilter}),'stable');
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
        % pathID, reportname (short), display unit,categorial text, report
        % name long (in Caption)
        PhysProperties.yList = {'Organism|Age', 'age', 'year(s)',{},'age';...
            'Organism|Weight', 'body weight', 'kg',{},'body weight';...
            'Organism|Height', 'body height', 'cm',{},'body height';...
            'Organism|BMI', 'BMI', 'kg/m²',{},'BMI';...
            'Organism|BSA', 'body surface area', 'm²',{},'body surface area';...
            'Gender', 'gender', '',{'Male','Female'},'gender';...
            '<addOntogenyFactor>','','',{},''}; % if the code finds the key word <addOntogenyFactor> 
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
            'Organism|Weight', 'body weight', 'kg'};
    otherwise
                error('unknown flag')

end


% set indices
switch flag
    
    case {'parallelComparison'}
        % list indices of popSet for that the VPC should be done
        PKParameter.ixOfPopRunSets =  1:length(PopRunSet);
                
        % name of Reference Population
        PKParameter.ixPopRunSetRef = [];
                
    case {'pediatric','ratioComparison'}

        
        % list popset for that the VPC should be done
        jj_ref = [PopRunSet.isReference];
        PKParameter.ixOfPopRunSets = find(~jj_ref) ;
                
        % name of Reference Population
        PKParameter.ixPopRunSetRef = find(jj_ref);
    otherwise
        error('unknown flag')
  
end


% do also boxwhsiker for pediatrics
if strcmp('pediatric',flag)
    PKParameter(2) = PKParameter(1);
    PKParameter(2).xList={};
    PKParameter(2).ixOfPopRunSets = 1:length(PopRunSet);
    PKParameter(2).ixPopRunSetRef = [];
end


return
