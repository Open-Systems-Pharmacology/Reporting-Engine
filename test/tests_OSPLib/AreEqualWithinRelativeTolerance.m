function success = AreEqualWithinRelativeTolerance(y1, y2, relTol)
    success = 0;
    
    if length(y1) ~= length(y2)
        return;
    end
    
    compare1 = abs(y1-y2) <= abs(y1)*relTol;
    compare2 = abs(y1-y2) <= abs(y2)*relTol;
    
    compare = compare1 & compare2;
    
    success = all(compare);