function testString = removeForbiddenLetters(testString)
% REMOVEFORBIDDENLETTERS replaces letters not suited for filenames with '_'
%
% Inputs
%   testSTring (string) string which may contain letters not suitable
% 
% Outputs 
%   testString (string) with replade letters

% Open Systems Pharmacology Suite;  http://forum.open-systems-pharmacology.org
% Date: 14-July-2017

jj = ismember(testString,'./ ?§$%&()[]{}+~*#');
testString(jj)='_';

return
