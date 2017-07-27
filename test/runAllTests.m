function runAllTests

testList = dir('test_NO*');

maindir = cd;

%new testdirectory
testDir = sprintf('test_%s',datestr(now,'yyyy_mm_dd_HHMM'));
mkdir(testDir);

logfile = fullfile(maindir,testDir,'logfile.txt');
        

for iTest = 1:length(testList)

    try

        writeToLog(sprintf('Start Test %d',iTest),logfile,true,false);
        
        logfile = fullfile(maindir,testDir,'logfile.txt');

        sourceDir = fullfile(maindir,testList(iTest).name,'Input');
        targetDir = fullfile(maindir,testDir,testList(iTest).name);
        
        copyfile(sourceDir,targetDir);
        
        cd(targetDir);
        workflow;
        if exist('checkTestResult.m','file')
            checkTestResult;
        else
            writeToLog(sprintf('No automatic test available for test %d',iTest),logfile,true,false);
        end
        
        
    catch exception
        writeToLog(sprintf('Test %d failed: %s',iTest,exception.message),logfile,true,false);
    end
    
    cd(maindir);
end
    
    
