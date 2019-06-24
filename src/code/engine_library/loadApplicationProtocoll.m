function [ApplicationProtocol,isValid,PKParameterTemplate] = loadApplicationProtocoll(RunSetName)

% LOADOUTPUTLIST loads the temporary Outputlist
%
% Inputs: 
%       RunSetName (string)  name of Runset
% 
%   Outputs:
%   ApplicationProtocol structure with following fields
%       name                (string)    name of application
%       startTime           (double)    [min] starting time of application
%       dose                (double)    [kg]  applicated drug amount
%       dosePerBodyWeight   (double)    [kg/kg] applicated drug amount per body weight
%       drugMass            (double)    [µmol] applicated drug amount
%       isDosePerBodyweight (boolean)   if true, dose is adjusted to bodyweight
%                                       false dose is absolute  
%       infusionTime        (double)    [min] time of infusion, 
%                                       zero if the application is no infusion
%       compound        (string)        name of applicated compound
%       isValid         (boolean)  if true,  all necessary field were available in the xml file, 
%                           if false, some fields could not be found, may be the xml is a 
%                               MoBi file with a user defined application. Special attention is needed 
%      PKParameterTemplate (structure) defines PK Parameter defined for
%                                   this  Application protocoll

% Open Systems Pharmacology Suite;  http://open-systems-pharmacology.org

tmp = load(fullfile('tmp',RunSetName,'applicationProtocol.mat'));

ApplicationProtocol = tmp.ApplicationProtocol;
isValid = tmp.isValid;
PKParameterTemplate = tmp.PKParameterTemplate;

return
