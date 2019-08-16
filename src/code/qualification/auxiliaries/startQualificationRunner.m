function startQualificationRunner(qualificationRunnerPath, qualificationPlan, outputFolder, varargin)
% starts the qualification runner and creates inputs for the reporting engine
%
% Arguments:
%  qualificationRunnerPath: folder where QualificationRunner.exe is located
%
%  qualificationPlan:       full path of the input qualification plan
%
%  outputFolder:            folder where outputs of the qualification runner (= inputs for RE) will be created
%                           If <outputFolder> is not empty and the "force" option (s. below) is not set: 
%                           qualificationRunner will exit with an error
%
%  optional arguments:
%    '-f' or '--force':    If set to true, the contents of the output folder will be deleted, even if it not empty. Default is false
%
%    '-n <name>':          Name of the report qualification plan to be generated. Default is 'report-configuration-plan.json'
%
%    '-p <PKSim Folder>':  Path of PK-Sim installation folder. If not specified, installation path will be read
%                          from the registry (available only in case of full (non-portable) installation).
%                          This option is MANDATORY for the portable version of PK-Sim
%
%    '-l <logFile>':       Full path of log file where log output will be written. A log file will not be created if this value is not provided.
%
%    '-a' or '--append':   true to append data to the file; false to overwrite the file (default). 
%                          If the specified file does not exist, this parameter has no effect, and a new file is created.
%
%    '--logLevel <Level>': Log verbosity (Debug, Information, Warning, Error). Default is Information.
%
%    '--version':          Display version information.

    optionalArguments = '';
    for i=1:length(varargin)
        optionalArguments = [optionalArguments ' ' varargin{i}]; %#ok<AGROW>
    end
    
    qualificationRunner = ['"' qualificationRunnerPath filesep 'QualificationRunner.exe"'];
    
    arguments = [' -i "' qualificationPlan '" -o "' outputFolder '"' optionalArguments];
    
    status = system([qualificationRunner arguments]);
        
    if status~=0
        error('QualificationRunner failed');
    end
    