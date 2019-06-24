function writeTabCellArray(array,fname)
% WRITETABCELLARRAY write a cell array of strings as table

% Open Systems Pharmacology Suite;  http://open-systems-pharmacology.org

t = cell2table(array);
writetable(t,fname,'Delimiter',';','WriteVariableNames',false,'WriteRowNames',false);

return
