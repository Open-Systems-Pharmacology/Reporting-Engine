function success = checkTestResult

% read logfile
fid = fopen('logfile.txt');
tmp=textscan(fid,'%s','delimiter','\n');
C=tmp{1};
fclose(fid);

% check if error message for non existing output occured
ij = strfind(C{end-4},'ERROR: Outputpath "Organism|PeripheralVenousBlood|Hydroxy_Itraconazole|Plasma (Peripheral Venous Blood) with Typo" could not be found in model');
success = ~isempty(ij);

% check if error message for non consitent display unit occured
ij = strfind(C{end-3},'ERROR: For unit "µmol/l", there is no common dimension with display unit "cm"');
success = success & ~isempty(ij);

% check if error message for non existent display unit occured
ij = strfind(C{end-2},'ERROR: unit "typo" for seems to be no default OSPSuite unit');
success = success & ~isempty(ij);

% check if error message for non existent display unit occured
ij = strfind(C{end-1},'ERROR: unit "typo2" for seems to be no default OSPSuite unit');
success = success & ~isempty(ij);

% check if error message for non existent display unit occured
ij = strfind(C{end},'ERROR: For unit "µmol/l", there is no common dimension with display unit "h"');
success = success & ~isempty(ij);

return