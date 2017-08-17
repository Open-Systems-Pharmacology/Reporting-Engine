function MBS = getDefaultMassbalanceSettings(MeanModelSet) 
%GETDEFAULTVPCMEANMODELSETTINGS get WSettings for visual predicitve check
%
%  MBS = getDefaultMassbalanceSettings(MeanModelSet) 
%
% Inputs:
%       MeanModelSet (structure)   list of population simulations see GENERATEWORKFLOWINPUTFORMEANMODELSIMULATION

% Open Systems Pharmacology Suite;  http://open-systems-pharmacology.org


%display unit of time
MBS.displayUnit = 'h';

% defines time for pieplots (-1 is end of simulation time)
MBS.piePlotTime = -1;

% set list of compounds which should be no part of massbalance
MBS.excludedCompounds = {};

