function writeTabCellArray(array,fname)
% WRITETABCELLARRAY write a cell array of strings as table

% Open Systems Pharmacology Suite;  http://open-systems-pharmacology.org



data = cell(3,size(array,2));

for iCol = 1:size(array,2)
    data{1,iCol} = sprintf('column %d',iCol); 
    data{2,iCol} = 'string';
    jj = cellfun(@isnumeric,array(:,iCol));
    if any(jj)
        array(jj,iCol) = cellfun(@num2str,array(jj,iCol),'uniformoutput',false);
    end
    data{3,iCol} = array(:,iCol);
end

writetab(fname, data, ';', 0, 0, 0,0);


return