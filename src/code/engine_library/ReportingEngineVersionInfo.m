function Info = ReportingEngineVersionInfo
% configuration function, here properties of the reporting engine are listed
%
% Info = ReportingEngineVersionInfo
%
% Outputs:
%   Info (structure) with following fields
%       ReportingEngineVersion (string)  current version of the reporting engine package
%       ListOfValidatedComputers (cellarray of strings) list of computer
%               names which are validated
%

% Open Systems Pharmacology Suite;  http://open-systems-pharmacology.org


% version is still under construction
Info.ReportingEngineVersion = 3.0; 

% add here the names of the computers, for which the code was validated
Info.ListOfValidatedComputers = {'BY-SPMREPRD', 'BY0S3H'}; 