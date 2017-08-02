function Output = getDefaultOutput(pathID,nameInReport,displayUnit)

%GETDEFAULTOUTPUT defines the properties of an Output
%
% Inputs:  
%   - pathID    (string)    path to identify simulation parameter
%   - nameInReport   (string)    used for graphics (please start with upper case)
%   - displayUnit (string)  unit used for graphical display (must be a unit, known in PK-Sim) 
%           
% Outputs:
%   - Output (structure) with following fields
%       - pathID    (string)    path to identify simulation parameter
%       - nameInReport   (string)    used for graphics (please start with upper case)
%       - displayUnit (string)  unit used for graphical display (must be a unit, known in PK-Sim) 
%       - unitFactor (double) factor to convert internal unit to target
%           unit is calculated during initialisation
%       - dataTpFilter (string) this is an optional input, it may be empty
%          if timeprofile data are available, and a filter per Output is necessary 
%          please enter here the filter conditions. 
%          use nonmem column headers as variable name.
%          e.g. 'CMT==1 & EVID=1'


% Open Systems Pharmacology Suite;  http://forum.open-systems-pharmacology.org
% Date: 21-July-2017



% list of all outputs defined as structure

Output.pathID = pathID;
% nameInReport used for graphics (please start with upper case)
Output.nameInReport = nameInReport;
% displayUnit unit used for graphical display (must be a unit, known in PK-Sim) 
Output.displayUnit = displayUnit;

% timeprofile filter for nonmem data file
Output.dataTpFilter = '';

% factor to convert interanl unit to target unit
Output.unitFactor = nan;

return
