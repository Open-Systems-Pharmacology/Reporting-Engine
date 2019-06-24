function success = checkResult
% automatic testscript for population test 9

success = true;

try
       
    % check number of files
    
    figureSubdirectories = { 'physiology' ,'timeprofile', 'pKParameter' };
    extension = {'csv','fig','emf'};
    nFile = [28 24 24; 8 8 8; 7 12 12];
        
    for iDir = 1:length(figureSubdirectories)
        
        for iExt = 1:length(extension)
            tmp = dir(fullfile('figures',figureSubdirectories{iDir},['*.',extension{iExt}]));
            
            if length(tmp) ~=  nFile(iDir,iExt)
                
                disp(sprintf('there must be %d %s files in directory %s, there are %d',nFile(iDir,iExt),extension{iExt},figureSubdirectories{iDir},length(tmp)));
                success = false;
                return
            end
        end
    end
    
       
    
catch exception
    
    disp(exception.message)
    success = false;
end

disp('automatic test was successful');

return