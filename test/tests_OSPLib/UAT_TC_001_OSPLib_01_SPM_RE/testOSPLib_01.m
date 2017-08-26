function success = testOSPLib_01(pathToReportingEngine, simulationFile)

success = 0;

%setup environment
addpath(genpath(pathToReportingEngine));
libPath = [pathToReportingEngine filesep 'lib'];
setenv('path', [libPath ';' getenv('path')]);

%Create the component
comp=DCIMatlabR2013b6_0('LoadComponent', [libPath filesep 'OSPSuite_SimModelComp.xml']);

%Get parameter table
tab=DCIMatlabR2013b6_0('GetParameterTable',comp,1);

%Simulation Schema
tab.Variables(1).Values={[libPath filesep 'OSPSuite.SimModel.xsd']};

%Simulation File
tab.Variables(2).Values={simulationFile};

%Set parameter table into the component and calculate
DCIMatlabR2013b6_0('SetParameterTable',comp,1,tab);
DCIMatlabR2013b6_0('Configure',comp);
DCIMatlabR2013b6_0('ProcessMetaData',comp);
DCIMatlabR2013b6_0('ProcessData',comp);

%Get 1st Output table (Simulation Times)
outTabTimes=DCIMatlabR2013b6_0('GetOutputTable',comp,1);

%Get 2nd Output table (Simulation values)
outTabValues=DCIMatlabR2013b6_0('GetOutputTable',comp,2);

%set output time vector
t = outTabTimes.Variables(1).Values;

%set all system outputs (ode variables and observers)
for i=1:length(outTabValues.Variables)
    values(i).Name = outTabValues.Variables(i).Name;
    values(i).Values = outTabValues.Variables(i).Values;
end

%================= verify model outputs =========================

%relative tolerance for values comparison
relTol = 1e-5;

% output time vector must be equal to [0,1,2,3,4,5,6,7,8,9,10] within relative tolerance
if ~AreEqualWithinRelativeTolerance(t, [0:10]', relTol)
    error('output time vector does not match');
end

% Variable values must have length 4
if length(values) ~= 4
    error('number of model outputs does not match');
end

% values(1).Name equals to 'y1'
if ~strcmp(values(1).Name, 'y1')
    error('Name of the 1st output does not match y1');
end

% values(1).Values equals to exp(t)+exp(-t) within relative tolerance
if ~AreEqualWithinRelativeTolerance(values(1).Values, exp(t)+exp(-t), relTol)
    error('output time vector for y1 does not match');
end

% values(2).Name equals to 'y2'
if ~strcmp(values(2).Name, 'y2')
    error('Name of the 2nd output does not match y2');
end

% values(2).Values equals to exp(t)-exp(-t) within relative tolerance
if ~AreEqualWithinRelativeTolerance(values(2).Values, exp(t)-exp(-t), relTol)
    error('output time vector for y2 does not match');
end

% values(3).Name equals to 'y3'
if ~strcmp(values(3).Name, 'y3')
    error('Name of the 3rd output does not match y3');
end

% values(3).Values is constant value 2 within relative tolerance
if ~AreEqualWithinRelativeTolerance(values(3).Values, 2, relTol)
    error('output time vector for y3 does not match');
end

% values(4).Name equals to 'Obs1'
if ~strcmp(values(4).Name, 'Obs1')
    error('Name of the 4th output does not match Obs1');
end

% values(4).Values equals to 2* values(1).Values within relative tolerance
if ~AreEqualWithinRelativeTolerance(values(4).Values, 2* values(1).Values, relTol)
    error('output time vector does not match');
end

%all tests successfull
success = 1;