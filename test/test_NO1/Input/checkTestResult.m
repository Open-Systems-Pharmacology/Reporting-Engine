function success = checkTestResult

% read logfile
fid = fopen('logfile.txt');
tmp=textscan(fid,'%s','delimiter','\n');
C=tmp{1};
fclose(fid);

% check if error message for non existing xml occurd
jj = cellfun(@(x) contains(x,'ERROR: "Children_OralSingle_IV_Multi_withTypo.xml" does not exist'),C);
success = any(jj);

% check if error message for non existing csv
jj = cellfun(@(x) contains(x,'ERROR: "Children_OralSingle_IV_Multi_withTypo.csv" does not exist'),C);
ij = strfind(C{end-2},'ERROR: "Children_OralSingle_IV_Multi_withTypo.csv" does not exist');
success = success & ~isempty(ij);



return
