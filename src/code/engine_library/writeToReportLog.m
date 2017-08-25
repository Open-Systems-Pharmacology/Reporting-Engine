function writeToReportLog(type,logText,isNewFile)
%WRITETOREPORTLOG Support function: Writes text to logfile. Each entry starts with a time stamp
%
%   writeToReportLog(type,logText,isNewFile)
%       type (string) type of message, possible are INFO,WARNRING, ERROR    
%       logText (string): text which is added to logFile
%       isNewFile (Boolean): 
%           false: (default)  text is appended to the existing logfile
%           true:   a new file is created 
%

% Open Systems Pharmacology Suite;  http://open-systems-pharmacology.org


disp(sprintf('%s: %s',type,logText));

% check inputs
if ~exist('isNewFile','var')
    isNewFile=false;
end


%get filename
switch type
    case  'INFO'
        logFile = 'logfile.txt';
    case {'WARNING', 'ERROR'}    
        logFile = 'errorLogfile.txt';
        logText = sprintf('%s: %s',type,logText);
    otherwise
        error('unknown filetype');
        
end


if isNewFile
    fid = fopen(logFile,'w');
else
     fid = fopen(logFile,'a');
end
fprintf(fid,'%s',[datestr(now)  ' ' logText]);
fprintf(fid,'\r\n');

if ismember(type,{'WARNRING', 'ERROR'})
    [st] = dbstack();
    for iSt = length(dbstack):-1:1
        fprintf(fid,'   %s line %d',st(iSt).file,st(iSt).line);
        fprintf(fid,'\r\n');
    end        
end

fclose(fid);

return
