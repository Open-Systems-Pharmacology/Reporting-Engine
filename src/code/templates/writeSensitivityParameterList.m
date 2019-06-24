function writeSensitivityParameterList(fid,sensParameterList)
% WRITEDEFINITIONSENSITIVITYANALYSIS write sensitivity parameter to workflow.m
%
% writeDefinitionSensitivityAnalysis(fid,sensParameterList)
%
% Inputs:
%   fid (double) id of file
%   sensParameterList (cellarray) list of parameter for sensitivity analysis

% Open Systems Pharmacology Suite;  http://open-systems-pharmacology.org


fprintf(fid,'%% List of Parameters for sensitivity analysis:');
fprintf(fid,'\r\n');
fprintf(fid,'%% set to {}, if not needed');
fprintf(fid,'\r\n');
fprintf(fid,'%%  columns: 1. path, 2. number of steps, 3. variation range, 4. minValue 5. maxValue');
fprintf(fid,'\r\n');

if isempty(sensParameterList)
    fprintf(fid,'sensParameterList = {};');
    fprintf(fid,'\r\n');
else

    for iF = 1:size(sensParameterList,1)
        fprintf(fid,'sensParameterList(%d,:) = {''%s'',%g,%g,''%s''};',iF,...
            sensParameterList{iF,1},sensParameterList{iF,2},sensParameterList{iF,3},...
            sensParameterList{iF,4});
        fprintf(fid,'\r\n');
    end
end
fprintf(fid,'\r\n');

return
