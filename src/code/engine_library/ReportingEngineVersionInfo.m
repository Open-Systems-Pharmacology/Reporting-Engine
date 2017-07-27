function Info = ReportingEngineVersionInfo
% configuration function, here properties of the reporting engine are listed
%
% Outputs:
%   Info (structure) with following fields
%       ReportingEngineVersion (string)  current version of the reporting engine package
%       Settings (structure)  see GETDEFAULTWORKFLOWSETTINGS
%
%

% Open Systems Pharmacology Suite;  http://forum.open-systems-pharmacology.org
% Date: 19-July-2017

% version is still under construction
Info.ReportingEngineVersion = 0.1; 

% add here the names of the computers, for which the code was validated
Info.ListOfValidatedComputers = {}; %{'BY0LBF'}; 
