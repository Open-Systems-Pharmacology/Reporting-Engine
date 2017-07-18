function outputList = readOutputSelection(outputList_csv)
%READOUTPUTSELECTION reads list of selected outputs from csv file
%
% Inputs: 
%   outputList_csv (string) name of the csv file containing the outputs
% Outputs
%   outputs (cellarray of string) List of selected outputs


fid = fopen(outputList_csv);
C=textscan(fid,'%s',2,'delimiter',';');
outputList=C{1};
fclose(fid);

return
