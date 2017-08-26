function success = testOSPLib_02(pathToReportingEngine, simulationFile)

%ODE system used for test is taken from the CVODES example "A serial dense example: cvsRoberts_FSA_dns"
%s. https://computation.llnl.gov/sites/default/files/public/cvs_examples.pdf

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

%define column indices for input parameters tables (both all and variable)
ParamIDIdx=1; ParamPathIdx=2; ParamValueIdx=3; ParamUnitIdx=4; ParamIsFormulaIdx=5; 
ParamFormulaIdx=6; ParamDescriptionIdx=7; ParamIsVariableIdx=8; ParamCalculateSensitivityIdx=9;

%get variable parameters table
inTab2 = DCIMatlabR2013b6_0('GetInputTable',comp,2);

%calculate sensitivity for parameters p1, p2, p3
inTab2.Variables(ParamCalculateSensitivityIdx).Values(1:3)=[1 1 1];

%Save input table into the component and set sensitivity parameters
DCIMatlabR2013b6_0('SetInputTable',comp,2,inTab2);

DCIMatlabR2013b6_0('ProcessMetaData',comp);
DCIMatlabR2013b6_0('ProcessData',comp);

%Get 1st Output table (Simulation Times)
outTabTimes=DCIMatlabR2013b6_0('GetOutputTable',comp,1);

%Get 2nd Output table (Simulation values)
outTabValues=DCIMatlabR2013b6_0('GetOutputTable',comp,2);

%Get 3rd Output table (parameter sensitivities)
outTabSensitivities=DCIMatlabR2013b6_0('GetOutputTable',comp,3);

%set output time vector
t = outTabTimes.Variables(1).Values;

%number of sensitivity vectors should be 12: d(y1,y2,y3,Obs1)/d(p1,p2,p3)
if length(outTabSensitivities.Variables)~=12
    error('not all expected sensitivities are calculated');
end

%get all sensitivities
for i=1:length(outTabSensitivities.Variables)
    sensitivities(i).OutputName = outTabSensitivities.Variables(i).Name;
    sensitivities(i).Values = outTabSensitivities.Variables(i).Values;
    sensitivities(i).SensitivityParameterId = outTabSensitivities.Variables(i).Attributes(5).Value;
end

%================= verify model outputs =========================

%relative tolerance for values comparison
relTol = 1e-3;

% output time vector must be equal to [0 0.4 4] within relative tolerance
if ~AreEqualWithinRelativeTolerance(t, [0 0.4 4]', relTol)
    error('output time vector does not match');
end

%now check sensitivity values. 
%Expected values are taken from the corresponding CVODES example cvsRoberts_FSA_dns (first 2 time steps )
%s. https://computation.llnl.gov/sites/default/files/public/cvs_examples.pdf

%dy1/dp1
checkSensitivityValues(sensitivities(1), 'y1', '111', [0 -3.5595e-01 -1.8761e+00]', relTol);

%dy1/dp2
checkSensitivityValues(sensitivities(2), 'y1', '112', [0 9.5431e-08 2.9614e-06]', relTol);

%dy1/dp3
checkSensitivityValues(sensitivities(3), 'y1', '113', [0 -1.5833e-11 -4.9334e-10]', relTol);

%dy2/dp1
checkSensitivityValues(sensitivities(4), 'y2', '111', [0 3.9025e-04 1.7922e-04]', relTol);

%dy2/dp2
checkSensitivityValues(sensitivities(5), 'y2', '112', [0 -2.1309e-10 -5.8305e-10]', relTol);

%dy2/dp3
checkSensitivityValues(sensitivities(6), 'y2', '113', [0 -5.2900e-13 -2.7626e-13]', relTol);

%dy3/dp1
checkSensitivityValues(sensitivities(7), 'y3', '111', [0 3.5556e-01 1.8759e+00]', relTol);

%dy3/dp2
checkSensitivityValues(sensitivities(8), 'y3', '112', [0 -9.5218e-08 -2.9608e-06]', relTol);

%dy3/dp3
checkSensitivityValues(sensitivities(9), 'y3', '113', [0 1.6362e-11 4.9362e-10]', relTol);

%observer Obs 1 is defined as y1+2*y2+3*y3
%so the same formula must apply for its sensitivities

dy1_dp1=sensitivities(1).Values; dy1_dp2=sensitivities(2).Values; dy1_dp3=sensitivities(3).Values;
dy2_dp1=sensitivities(4).Values; dy2_dp2=sensitivities(5).Values; dy2_dp3=sensitivities(6).Values;
dy3_dp1=sensitivities(7).Values; dy3_dp2=sensitivities(8).Values; dy3_dp3=sensitivities(9).Values;

%dObs1/dp1
checkSensitivityValues(sensitivities(10), 'Obs1', '111', [0 dy1_dp1(2)+2*dy2_dp1(2)+3*dy3_dp1(2) dy1_dp1(3)+2*dy2_dp1(3)+3*dy3_dp1(3)]', relTol);

%dObs1/dp2
checkSensitivityValues(sensitivities(11), 'Obs1', '112', [0 dy1_dp2(2)+2*dy2_dp2(2)+3*dy3_dp2(2) dy1_dp2(3)+2*dy2_dp2(3)+3*dy3_dp2(3)]', relTol);

%dObs1/dp3
checkSensitivityValues(sensitivities(12), 'Obs1', '113', [0 dy1_dp3(2)+2*dy2_dp3(2)+3*dy3_dp3(2) dy1_dp3(3)+2*dy2_dp3(3)+3*dy3_dp3(3)]', relTol);


%all tests successfull
success = 1;

function checkSensitivityValues(sensitivity, expectedName, expectedSensitivityParameterId, expectedSensitivityValues, relTol)

    if ~strcmp(sensitivity.OutputName, expectedName)
        error(['Name of the output does not match ' expectedName]);
    end
    
    if ~strcmp(sensitivity.SensitivityParameterId, expectedSensitivityParameterId)
        error(['Sensitivity parameter id does not match ' expectedSensitivityParameterId]);
    end
    
    if ~AreEqualWithinRelativeTolerance(sensitivity.Values, expectedSensitivityValues, relTol)
        error(['Sensitivity values of ' expectedName ' do not match']);
    end
    