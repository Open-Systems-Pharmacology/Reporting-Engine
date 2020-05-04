function addLegendWithLineBreaks(currentFig, curvesHandle, curvesLegend)

set(currentFig, 'unit', 'inches');
h_legend = legend(curvesHandle, cellfun(@(x) sprintf(x),curvesLegend,'UniformOutput',false));
set(h_legend, 'unit', 'inches');
nCharPerLine = max(cellfun(@(x) length(x),curvesLegend));
while 1.2*h_legend.Position(3) > currentFig.Position(3)
    nCharPerLine = nCharPerLine-1;
    for ll = 1:size(curvesLegend,2)
        curvesLegend_lineBreaks{ll} = insertLineBreaks(curvesLegend{ll},nCharPerLine);
    end
    h_legend = legend(curvesHandle, cellfun(@(x) sprintf(x),curvesLegend_lineBreaks,'UniformOutput',false));
end

function res = insertLineBreaks(str, NoCharPerLine)

if size(str,2) > NoCharPerLine
    II1=strfind(str,' ');
    II2=strfind(str,'-');
    II = sort([II1 II2]);
    for ii = 0:(length(II)-1)
        if size(str(1:II(end-ii)),2) <= NoCharPerLine
            str1 = str(1:II(end-ii));
            str2 = str((II(end-ii)+1):end);
            res = sprintf('%s%s%s',str1,'\n',str2);
            break
        end
    end
else
    res = str;
end
