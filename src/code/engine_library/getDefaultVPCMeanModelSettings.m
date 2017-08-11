function VPC = getDefaultVPCMeanModelSettings(MeanModelSet) 
%GETDEFAULTVPCMEANMODELSETTINGS get WSettings for visual predicitve check
%
%  VPC = getDefaultVPCMeanModelSettings(WSettings,MeanModelSet,flag)
%
% Inputs:
%       MeanModelSet (structure)   list of population simulations see GENERATEWORKFLOWINPUTFORMEANMODELSIMULATION

% Open Systems Pharmacology Suite;  http://open-systems-pharmacology.org


% add function handle for figure and legend text
VPC.textFunctionHandle = @textVPCMeanModel;

% add default structure for timeprofile plots
VPC.Timeprofile = addVPCTimeprofile(MeanModelSet);

return


function Timeprofile = addVPCTimeprofile(MeanModelSet)

% list indices of MeanModelSet for that the VPC should be done
Timeprofile.ixOfRunSets = 1:length(MeanModelSet);

% display Unit for tim
Timeprofile.timeDisplayUnit = 'h';

% zoom int time axes
Timeprofile.timelimit = [];

% set scale of y axes
Timeprofile.yScale = {'lin','log'};

% get Types of Plot
Timeprofile.plotTypes = {'yVsTime','resVsTime','resVsY','predVsObs','histRes','qqPlotRes'};

return

