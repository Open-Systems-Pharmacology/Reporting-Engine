function success = checkTestResult

% initialize success
success = true;
% compare simulation results
% New generated file
data = readtab(fullfile('simulations','OralSingle_IV_Multi-Results.csv'),';',0,0,1,0);
% old comparison file
dataComp  = readtab(fullfile('simulationExportPKSIm','Children_OralSingle_IV_Multi-Results.csv'),';',0,0,1,0);

if ~all(size(data) == size(dataComp))
    success = false;
    return
end

for iD = 1:size(data,2)
    
    xNew = data{3,iD};
    xOld = dataComp{3,iD};
    
    jjZero = xOld ==0;
    
    % compare zeros
    if any(abs(xNew(jjZero)) >= 1e-20)
        success = false;
    end
    
    % compare relative difference of non Zeros
    
    diff = abs(xNew(~jjZero) - xOld(~jjZero))./xNew(~jjZero);
    
    if any(diff >= 1e-5)
        success = false;
        return
    end
    
end

% compare PK Parameter
% New generated file
data = readtab(fullfile('simulations','OralSingle_IV_Multi-PK-Analyses.csv'),';',0,0,1,0);
% old comparison file
dataComp  = readtab(fullfile('simulationExportPKSIm','Children_OralSingle_IV_Multi-PK-Analyses.csv'),';',0,0,1,0);
% get rid of strange letter befor µ
dataComp{3,5} = strrep(dataComp{3,5},'Â','');

if ~all(size(data) == size(dataComp))
    success = false;
    return
end

for iD = 1:size(data,2)
    
    switch data{2,iD}
        
        case 'double'
            xNew = data{3,iD};
            xOld = dataComp{3,iD};
            
            jjZero = xOld ==0;
            
            % compare zeros
            if any(abs(xNew(jjZero)) >= 1e-20)
                success = false;
                return
            end
            
            % compare relative difference of non Zeros
            
            diff = abs(xNew(~jjZero) - xOld(~jjZero))./xNew(~jjZero);
            
            if any(diff >= 1e-5)
                success = false;
                return
            end
            
        case 'string'
            if any(~strcmp(data{3,iD},dataComp{3,iD}))
                success = false;
                return
            end
    end
    
end

return