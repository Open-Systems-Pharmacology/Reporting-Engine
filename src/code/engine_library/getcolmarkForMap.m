function [col,mark] = getcolmarkForMap(map,nData)
%GETCOLMARKFORMAP get colormatrix and marker-vector for a specified colormap
% 
% [col,mark] = getcolmarkForMap(map,nData)
% 	
% Input:
%   map (matrix)   colormap
%   nData (double) number of different colors
% 	
% Output:
% 	col:	(nData x 3) - Matrix with colors
% 	mark:	Vector of  Markers

% Open Systems Pharmacology Suite;  http://open-systems-pharmacology.org

% get index for colormap
if nData >1
    step = (size(map,1)-1)/(nData-1);
    ix = round(1:step:size(map,1));
else
    ix = round(size(map,1)/2);
end

% make sure last point is not lost because of numeric
if length(ix) < nData
    ix(end+1) = size(map,1);
end

% get colormatrix
col = map(ix,:);

% markers
markers = 'ovs^d<p>h';
mark = [];
for i=1:ceil(nData/9)
	mark = strcat(mark,markers);
end
mark = mark(1:nData);

return
