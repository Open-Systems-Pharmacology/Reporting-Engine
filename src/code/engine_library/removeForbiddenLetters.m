function testString = removeForbiddenLetters(testString)
% REMOVEFORBIDDENLETTERS replaces letters not suited for filenames with '_'
%
% Inputs
%   testSTring (string) string which may contain letters not suitable
% 
% Outputs 
%   testString (string) with replade letters

% Open Systems Pharmacology Suite;  http://open-systems-pharmacology.org


jj = ismember(testString,'./ ?§$%&()[]{}+~*#:');
testString(jj)='_';

return
