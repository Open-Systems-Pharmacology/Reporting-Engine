function success = checkResult
% automatic testscript for test case 8

success = true;

try
   
    
    % check if xml file was loaded and interpreted correctly
    if ~exist(fullfile('tmp','SingleIvBolus','applicationProtocol.mat'),'file')
        success = false;
        disp('failed in loading xml simulation file')
        return
    end
    
    load(fullfile('tmp','SingleIvBolus','applicationProtocol.mat'));
    if isempty(ApplicationProtocol)
        success = false;
        disp('failed in loading xml simulation file')
        return
    end
    
    % cehck if population file was loaded correctly
    if ~exist(fullfile('tmp','SingleIvBolus','pop.mat'),'file')
        success = false;
        disp('failed in loading pop csv')
        return
    end
    
    load(fullfile('tmp','SingleIvBolus','pop.mat'));
    if length(parPaths)~=104
        success = false;
        disp('failed in loading pop csv, length parPaths not ok')
        return
    end

    if abs(mean(mean(parValues(:,5:end)))-48.0348)>0.01
         success = false;
        disp('failed in loading pop csv, mean of parValues not ok')
        return
    end    
    
catch exception
    
    disp(exception.message)
    success = false;
end

disp('automatic test was successful');
