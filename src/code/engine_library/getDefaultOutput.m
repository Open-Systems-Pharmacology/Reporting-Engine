function Output = getDefaultOutput(pathID,reportName,displayUnit)
%GETDEFAULTOUTPUT defines the properties of an Output
%
% Inputs:  
%   - pathID    (string)    path to identify simulation parameter
%   - reportName   (string)    used for graphics (please start with upper case)
%   - displayUnit (string)  unit used for graphical display (must be a unit, known in PK-Sim) 
%           
% Outputs:
%   - Output (structure) with following fields
%       - pathID    (string)    path to identify simulation parameter
%       - reportName   (string)    used for graphics (please start with upper case)
%       - displayUnit (string)  unit used for graphical display (must be a unit, known in PK-Sim) 
%       - unitFactor (double) factor to convert internal unit to target
%           unit is calculated during initialisation
%       - dataTpFilter (string) this is an optional input, it may be empty
%          if timeprofile data are available, and a filter per Output is necessary 
%          please enter here the filter conditions. 
%          use nonmem column headers as variable name.
%          e.g. 'CMT==1 & EVID=1'
%       - pKParameterList (cellarray)
%                   first line identifer of PK Parameter
%                   second line display Unit
%                   third line unit factor to transfer from base unit to
%                   display unit


% Open Systems Pharmacology Suite;  http://open-systems-pharmacology.org




% list of all outputs defined as structure

Output.pathID = pathID;
% reportName used for graphics (please start with upper case)
Output.reportName = reportName;
% displayUnit unit used for graphical display (must be a unit, known in PK-Sim) 
Output.displayUnit = displayUnit;

% timeprofile filter for nonmem data file
Output.dataTpFilter = '';

% factor to convert interanl unit to target unit
Output.unitFactor = nan;

% list of PK parameter
Output.pKParameterList = {};

return
