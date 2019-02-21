function [NumberOfSuccessfulTests, SuccessfulTestsInfo, NumberOfFailedTests, FailedTestsInfo] = runAllTests

NumberOfSuccessfulTests = 0;
SuccessfulTestsInfo = [];
NumberOfFailedTests = 0;
FailedTestsInfo = [];

testList = dir('test_NO*');

maindir = cd;

%new test directory
testDir = sprintf('test_%s',datestr(now,'yyyy_mm_dd_HHMM'));
mkdir(testDir);

logfile = fullfile(maindir,testDir,'logfile.txt');
        
% loop on tests
for iTest = 1:length(testList)

    testInfo.TestName = testList(iTest).name;
    testInfo.ErrorMessage = '';

    try
        writeToLog(sprintf('Start %s',testList(iTest).name),logfile,true,false);
        
        logfile = fullfile(maindir,testDir,'logfile.txt');

        sourceDir = fullfile(maindir,testList(iTest).name,'Input');
        targetDir = fullfile(maindir,testDir,testList(iTest).name);
        copyfile(sourceDir,targetDir);
        
        cd(targetDir);
        workflow;
        
        if exist('exception.mat','file')
            NumberOfFailedTests = NumberOfFailedTests + 1;
            testInfo.ErrorMessage = 'S. exception.mat for details';
            writeToLog(sprintf('%s exception occured',testList(iTest).name),logfile,true,false);
        elseif exist('checkResult.m','file')
            success = checkResult;
            if success
                NumberOfSuccessfulTests = NumberOfSuccessfulTests + 1;
                writeToLog(sprintf('%s was successful',testList(iTest).name),logfile,true,false);
            else
                NumberOfFailedTests = NumberOfFailedTests + 1;
                testInfo.ErrorMessage = 'checkResult failed';
                writeToLog(sprintf('%s automatic test failed',testList(iTest).name),logfile,true,false);
            end
        else
            NumberOfSuccessfulTests = NumberOfSuccessfulTests + 1;
            writeToLog(sprintf('No automatic test available for %s',testList(iTest).name),logfile,true,false);
        end
        
    catch exception
        NumberOfFailedTests = NumberOfFailedTests + 1;
        testInfo.ErrorMessage = exception.message;
        writeToLog(sprintf('%s crashed: %s',testList(iTest).name,iTest,exception.message),logfile,true,false);
    end
    
    if NumberOfFailedTests > length(FailedTestsInfo)
        %test run failed - attach test info to the list of FAILED tests
        if isempty(FailedTestsInfo)
            FailedTestsInfo = testInfo;
        else
            FailedTestsInfo(end+1) = testInfo; %#ok
        end
    else
        %test run succeeded - attach test info to the list of SUCESSFUL tests
        if isempty(SuccessfulTestsInfo)
            SuccessfulTestsInfo = testInfo;
        else
            SuccessfulTestsInfo(end+1) = testInfo; %#ok
        end
    end
    
    close all;
    cd(maindir);
end
    
    
