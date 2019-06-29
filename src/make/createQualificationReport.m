function createQualificationReport(qualificationWorkflow)

% the following lines are interpreted by the matlab compiler mcc, at compile time
% WINDOWS: DCIMatlabR2017b6_1Copy.dll must be present in matlab path
%#function DCIMatlabR2017b6_1Copy
% LINUX: DCIMatlabR2017b6_1.mexglx must be present in matlab path
%#function DCIMatlabR2017b6_1

%    slCharacterEncoding();
%    which('MoBiSettings')
    
    libPath = fullfile(fileparts(which('MoBiSettings')),'../../lib');
    setenv('path', [libPath ';' getenv('path')]);
%    getenv('path')
    simModelSchemaPath = fullfile(libPath, 'OSPSuite.SimModel.xsd');
    simModelCompConfigPath = fullfile(libPath, 'OSPSuite_SimModelComp.xml');
    MoBiSettings(simModelSchemaPath, simModelCompConfigPath);

    invalidFileIdentifier = -1;
    fid = invalidFileIdentifier;

    try
        [fid, message] = fopen(qualificationWorkflow);

        if fid==invalidFileIdentifier
            disp(message);
            return;
        end

        tline = fgetl(fid);
        
        while ischar(tline)
            tline = strtrim(tline);
            
            if isempty(tline) || strcmp(tline(1),'%') || strcmp(tline, 'clear all')
                tline = fgetl(fid);
                continue
            end
            
%            disp(tline);
            eval(tline);
            tline = fgetl(fid);
        end
        
    catch ex
        if fid ~= invalidFileIdentifier
            fclose(fid);
        end
        rethrow(ex);
    end
    