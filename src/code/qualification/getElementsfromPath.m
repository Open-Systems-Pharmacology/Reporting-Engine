function C = getElementsfromPath(Path, p)
% GETELEMENTFROMPATH: Get all the elements from a path
%
% Path (string) input path to be read
%
% C (cells) array of cell with elements with path delemited by delemiter
% p (string) optional, delimiter of the path


if ~exist('p')
    p=object_path_delimiter;
end

delims = strfind(Path, p);

if isempty(delims)
    C={Path};
else
    delims = [0 delims length(Path)+1];
    for i=1:length(delims)-1
        C{i}=Path(delims(i)+1:delims(i+1)-1);
    end
end