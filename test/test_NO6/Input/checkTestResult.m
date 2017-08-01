function success = checkTestResult

% initialize success
success = true;
% compare simulation results
% New generated file
data = readtab(fullfile('simulations','SingleIvBolus-Results.csv'),';',0,0,1,0);
% old comparison file
dataComp  = readtab(fullfile('simulationExportPKSIm','Europeans_IVSingle_BSA-Results.csv'),';',0,0,1,0);

if ~all(size(data) == size(dataComp))
    success = false;
    return
end

% check if all values are incraed by a factor 10
for iD = 3:size(data,2)
    
    xNew = data{3,iD};
    xOld = dataComp{3,iD}.*10;
    
    jjZero = xOld ==0;
    
    % compare zeros
    if any(abs(xNew(jjZero)) >= 1e-20)
        success = false;
    end
    
    % compare relative difference of non Zeros
    
    diff = abs(xNew(~jjZero) - xOld(~jjZero))./xNew(~jjZero);
    
    if any(diff >= 1e-4)
        success = false;
        return
    end
    
end


return