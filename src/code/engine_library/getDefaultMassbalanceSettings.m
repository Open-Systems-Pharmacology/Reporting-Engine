function MBS = getDefaultMassbalanceSettings
%GETDEFAULTMASSBALANCESETTINGS get WSettings for massbalance plots
%
% MBS = getDefaultMassbalanceSettings
%
% Inputs:
%
% Outputs
%   MBS.displayUnit (display unit of time)
%   MBS.piePlotTime (defines time for pieplots (-1 is end of simulation
%   time))
%   MBS.excludedCompounds = {};list of compounds which should be no part of massbalance

% Open Systems Pharmacology Suite;  http://open-systems-pharmacology.org


%display unit of time
MBS.displayUnit = 'h';

% defines time for pieplots (-1 is end of simulation time)
MBS.piePlotTime = -1;

% set list of compounds which should be no part of massbalance
MBS.excludedCompounds = {};

