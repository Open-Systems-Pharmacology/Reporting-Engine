function success = checkTestResult

% read logfile
fid = fopen('logfile.txt');
tmp=textscan(fid,'%s','delimiter','\n');
C=tmp{1};
fclose(fid);

% check if error message for non existing xml occurd
ij = strfind(C{end-2},'ERROR: "Children_OralSingle_IV_Multi_withTypo.xml" does not exist');
success = ~isempty(ij);

% check if error message for non existing csv
ij = strfind(C{end-1},'ERROR: "Children_OralSingle_IV_Multi_withTypo.csv" does not exist');
success = success & ~isempty(ij);


% check if error message for non existing csv
ij = strfind(C{end},'ERROR: "Output_withTypo.csv" does not exist');
success = success & ~isempty(ij);


return