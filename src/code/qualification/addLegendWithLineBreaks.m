function addLegendWithLineBreaks(currentFig, curvesHandle, curvesLegend)

set(currentFig, 'unit', 'inches');
h_legend = legend(curvesHandle, cellfun(@(x) sprintf(x),curvesLegend,'UniformOutput',false));
set(h_legend, 'unit', 'inches');
nCharPerLine = max(cellfun(@(x) length(x),curvesLegend));
nCharPerLineMax = nCharPerLine;
maxNoLines = 2;
while 1.2*h_legend.Position(3) > currentFig.Position(3)
    nCharPerLine = nCharPerLine-1;
    if nCharPerLine >= (nCharPerLineMax/maxNoLines)
        for ll = 1:size(curvesLegend,2)
            curvesLegend_lineBreaks{ll} = insertLineBreaks(curvesLegend{ll},nCharPerLine,maxNoLines);
        end
        h_legend = legend(curvesHandle, cellfun(@(x) sprintf(x),curvesLegend_lineBreaks,'UniformOutput',false));
    else
        nCharPerLine = nCharPerLineMax;
        maxNoLines = maxNoLines + 1;
    end
end

function res = insertLineBreaks(str, NoCharPerLine, maxNoLines)
residStr = str;
res = '';
for ll = 1:maxNoLines
    if length(residStr) > NoCharPerLine && ll < maxNoLines
        II1=strfind(residStr,' ');
        II2=strfind(residStr,'-');
        II = sort([II1 II2]);
        if ~isempty(II)
            for ii = 0:(length(II)-1)
                if size(residStr(1:II(end-ii)),2) <= NoCharPerLine
                    str1 = residStr(1:II(end-ii));
                    residStr = residStr((II(end-ii)+1):end);
                    res = sprintf('%s%s%s',res,str1,'\n');
                    break
                end
            end
        else
            str1 = residStr(1:floor(NoCharPerLine));
            residStr = residStr((ceil(NoCharPerLine)+1):end);
            res = sprintf('%s%s%s',res,str1,'\n');
        end
    else
        res = sprintf('%s%s',res,residStr);
        break
    end
end
