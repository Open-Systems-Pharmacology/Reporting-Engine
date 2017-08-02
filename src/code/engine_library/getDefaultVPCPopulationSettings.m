function VPC = getDefaultVPCPopulationSettings(Settings,PopRunSet,flag)
%GETDEFAULTVPCPOPULATIONSETTINGS get Settings for visual predicitve check
%
% Inputs:
%       Settings (structure)    definition of properties used in all
%                   workflow functions see GETDEFAULTWORKFLOWSETTINGS
%       PopRunSet (structure)   list of population simulations see GETDEFAULTPOPRUNSET
%       flag (string)  [optional] defines type of population, default = 'default'
%                               other possibilities are 'pediatric';

% Open Systems Pharmacology Suite;  http://forum.open-systems-pharmacology.org
% Date: 26-July-2017

% check inputs
if ~exist('flag','var')
    flag = 'default';
end

% demographic
VPC.demographic = addVPCDemographic(PopRunSet,flag);


return

function demographic = addVPCDemographic(PopRunSet,flag)

switch flag
    
    case 'default'
        
        demographic.name = 'demographic';
        demographic.reportName = '';
        
        % list demographic properties which should be plotted
        % pathID, reortname, display unit
        demographic.yList = {'Organism|Age', 'Age', 'year(s)';...
            'Organism|Weight', 'Body weight', 'kg';...
            'Organism|Height', 'Height', 'cm';...
            'Organism|BMI', 'BMI', 'kg/m²'};
        
        % list propety which shall be displayed on the x-axes, if empty a
        % histogramm will be created
        demographic.xList = {};
        
        % list indices of popSet for that the VPC should be done
        [~,ix] = unique({PopRunSet.popcsv},'stable');
        demographic.ixOfPopRunSets =  ix;
        
        % how to handle PopRunSets:
        % merged populations are merged,
        % serial, for each PopRunSet an extra plot is generated
        demographic.popRunSetMerger =  'serial';
        
        % name of Reference Population
        demographic.ixPopRunSetRef = [];
        
    case 'pediatric'

        demographic.name = 'pediatricPopulation';
        demographic.reportName = 'virtual pediatric population';

        
        % list demographic properties which should be plotted
        % pathID, reortname, display unit
        demographic.yList = {'Organism|Weight', 'Body weight', 'kg';...
            'Organism|Height', 'Height', 'cm';...
            'Organism|BMI', 'BMI', 'kg/m²'};
        
        % list propety which shall be displayed on the x-axes, if empty a
        % histogramm will be created
        demographic.xList = {'Organism|Age', 'age', 'year(s)'};
        
        % list popset for that the VPC should be done
        jj_ref = [PopRunSet.isReference];
        [~,ix_unique] = unique({PopRunSet.popcsv},'stable');
        ix = intersect(find(~jj_ref),ix_unique);
        demographic.ixOfPopRunSets = ix ;
        
        % how to handle PopRunSets:
        % merged populations are merged,
        % serial, for each PopRunSet an extra plot is generated
        demographic.popRunSetMerger =  'merged';
        
        % name of Reference Population
        ix = intersect(find(jj_ref),ix_unique);
        demographic.ixPopRunSetRef = ix(1);
        
end



return
