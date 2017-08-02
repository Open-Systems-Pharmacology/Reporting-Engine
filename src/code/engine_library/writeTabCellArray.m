function writeTabCellArray(x,fname)
% WRITETABCELLARRAY write a cell array of strings as table

% Open Systems Pharmacology Suite;  http://forum.open-systems-pharmacology.org
% Date: 27-July-2017

for iCol = 1:size(x,2)
    data{1,iCol} = sprintf('column %d',iCol);
    data{2,iCol} = 'string';
    jj = cellfun(@isnumeric,x(:,iCol));
    if any(jj)
        x(jj,iCol) = cellfun(@num2str,x(jj,iCol),'uniformoutput',false);
    end
    data{3,iCol} = x(:,iCol);
end

writetab(fname, data, ';', 0, 0, 0,0);


return