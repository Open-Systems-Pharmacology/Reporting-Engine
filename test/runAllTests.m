function runAllTests

testList = dir('test_NO*');

maindir = cd;

%new testdirectory
testDir = sprintf('test_%s',datestr(now,'yyyy_mm_dd_HHMM'));
mkdir(testDir);

logfile = fullfile(maindir,testDir,'logfile.txt');
        
% loop on tests
for iTest = 9:length(testList)

    try

        writeToLog(sprintf('Start %s',testList(iTest).name),logfile,true,false);
        
        logfile = fullfile(maindir,testDir,'logfile.txt');

        sourceDir = fullfile(maindir,testList(iTest).name,'Input');
        targetDir = fullfile(maindir,testDir,testList(iTest).name);
        preparationDir = fullfile(maindir,testList(iTest).name,'prepareInput');
        
        cd(preparationDir)
        if exist('Workflow.xlsx','file')
            preparePopulationWorkflow('Workflow.xlsx');
        else
            prepareMeanModelWorkflow('WorkflowMean.xlsx');
        end
        copyfile(fullfile(preparationDir,'workflow.m'),fullfile(sourceDir,'workflow.m'));
        
        
        
        copyfile(sourceDir,targetDir);
        
        cd(targetDir);
        workflow;
        if exist('exception.mat','file')
            success = false;
                writeToLog(sprintf('%s exception occured',testList(iTest).name),logfile,true,false);
        elseif exist('checkTestResult.m','file')
            success = checkTestResult;
            if success
                writeToLog(sprintf('%s was successful',testList(iTest).name),logfile,true,false);
            else
                writeToLog(sprintf('%s automatic test failed',testList(iTest).name),logfile,true,false);
            end
        else
            writeToLog(sprintf('No automatic test available for %s',testList(iTest).name),logfile,true,false);
        end
        
        
    catch exception
        writeToLog(sprintf('%s crashed: %s',testList(iTest).name,iTest,exception.message),logfile,true,false);
    end
    
    cd(maindir);
end
    
    
