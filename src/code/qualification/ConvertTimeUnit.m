function [Yout, YoutUnit] = ConvertTimeUnit(Yin, YinUnit, YoutUnit)

% List of Time units
TimeUnitList={'s', 'min', 'h', 'd', 'y'};
ScaleFactor=cumprod([1, 60, 60, 24, 365]);

kin=[]; kout=[];
for k=1:length(TimeUnitList)
    if ~isempty(strfind(YoutUnit, TimeUnitList{k}))
        kout=k;
    end
    if ~isempty(strfind(YinUnit, TimeUnitList{k}))
        kin=k;
    end
end

% If the units were found, then convert    
if ~isempty(kin) && ~isempty(kout)
    Yout = Yin*ScaleFactor(kin)/ScaleFactor(kout);
else
    Yout=Yin;
end
