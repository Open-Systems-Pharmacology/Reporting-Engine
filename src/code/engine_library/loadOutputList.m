function OutputList = loadOutputList(RunSetName)
% LOADOUTPUTLIST loads the temporary Outputlist
%
% Inputs: 
%       RunSetName (structure)   name of correspondedn runset 
% 
%   Outputs:
%       OutputList (structure)  lists of all defined outputs

% Open Systems Pharmacology Suite;  http://open-systems-pharmacology.org

tmp = load(fullfile('tmp',RunSetName,'outputList.mat'));

OutputList = tmp.OutputList;

return
