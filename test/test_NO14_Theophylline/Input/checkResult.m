function success = checkResult
% automatic testscript population testcase theophylline
success = true;

try
   
    % check number of files
    figureSubdirectories = { 'physiology' ,'timeprofile', 'pKParameter' };
    extension = {'csv','fig','emf'};
    nFile = [7 6 6; 3 7 7; 8 14 14];
    
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
    
    % check goodnss of Fir
    goodness = readtable(fullfile('figures','timeprofile','goodnessOfFit_PO320mg.csv'));
    
    % test if residuals are distributed around zero, and there is no bias
    % vs time and concentration
    if ~all(goodness.success(2:3))
        disp('Goodness of fit failes, check table goodnessOfFit_PO320mg');
        success = false;
        return
    end
    
catch exception
    
    disp(exception.message)
    success = false;
end

disp('automatic test was successful');

return