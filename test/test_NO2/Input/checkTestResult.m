function success = checkTestResult

% read logfile
fid = fopen('errorLogfile.txt');
tmp=textscan(fid,'%s','delimiter','\n');
C=tmp{1};
fclose(fid);

% check if error message for non existing output occured
jj = cellfun(@(x) contains(x,'ERROR: Outputpath "Organism|PeripheralVenousBlood|Hydroxy_Itraconazole|Plasma (Peripheral Venous Blood) with Typo" could not be found in model'),C);
success = any(jj);

% check if error message for non consitent display unit occured
jj = cellfun(@(x) contains(x,'ERROR: For unit "µmol/l", there is no common dimension with display unit "cm"'),C);
success = success & any(jj);

% check if error message for non existent display unit occured
jj = cellfun(@(x) contains(x,'ERROR: unit "typo" seems to be no default OSPSuite unit'),C);
success = success & any(jj);

% check if error message for non existent display unit occured
jj = cellfun(@(x) contains(x,'ERROR: unit "typo2" seems to be no default OSPSuite unit'),C);
success = success & any(jj);

% check if error message for non existent display unit occured
jj = cellfun(@(x) contains(x,'ERROR: For unit "µmol/l", there is no common dimension with display unit "h"'),C);
success = success & any(jj);

return
