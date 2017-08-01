function FP = plotLoopPopulationVPCPKparameter(WSettings,textFunctionHandle,Def,PopRunSet,FP) 
% PLOTLOOPPOPULATIONVPCPKPARAMETER generates pK Parameter plots for VPC
%
% FP = plotLoopPopulationVPCPKparameter(WSettings,Def,PopRunSet,FP) 
%
% Inputs: 
%       Settings (structure)    definition of properties used in all
%                   workflow functions see GETDEFAULTWORKFLOWSETTINGS
%       textFunctionHandle (function handle) function handle to set figure text
%       Def (structure) contains information which plots should be
%               generated  see GETDEFAULTVPCPOPULATIONSETTINGS
%               (timeprofile)
%       PopRunSet (structure)   list of population simulations see GETDEFAULTPOPRUNSET
%       FP (ReportFigurPrint)  objects with manages print of figures and
%                           tables
%   Outputs:
%       FP (ReportFigurPrint)  objects with manages print of figures and
%       tables

% Open Systems Pharmacology Suite;  http://forum.open-systems-pharmacology.org
% Date: 31-July-2017


% collect all simulated PK Parameter of interst
x = collectPK(PopRunSet(Def.ixOfPopRunSets));


if isempty(Def.xList)
    
    FP = doBWhiskers(WSettings,Def,textFunctionHandle,x,PopRunSet,FP);

else
    
    FP = doShadedArea(WSettings,Def,textFunctionHandle,x,PopRunSet,FP);
end

return

function FP = doBWhiskers(WSettings,Def,textFunctionHandle,x,PopRunSet,FP)


reportNames = {PopRunSet(Def.ixOfPopRunSets).reportName};
popReportNames = {PopRunSet(Def.ixOfPopRunSets).popReportName};

yscale = {'lin','log'};

for iO = 1:length(x)
    
    %initialize new sheet
    header = sprintf('PK parameter of %s',x(iO).reportName);
    sheet = removeForbiddenLetters(header);
    FP = FP.iniCaptiontextFigtextArray(header,sheet);
    
    for iPK = 1:size(x(iO).pKParameterList,2)
         
        y = squeeze(x(iO).values(:,:,iPK));
        yLabel = x(iO).pKParameterList{4,iPK};
        yUnit = x(iO).pKParameterList{2,iPK};
        
        % loop on scale
        for iScale = 1:length(yscale)
            
            % get name and figure description
            figureName = sprintf('bw%d_%s',iPK,removeForbiddenLetters(x(iO).pKParameterList{1,iPK}));
            [figtxt,figtxtTable,legendEntries] = feval(textFunctionHandle,WSettings,'pkBW',...
                {yLabel,x(iO).reportName,reportNames,popReportNames,yscale{iScale}});
            
            % do figure
            [csv] = plotReportBoxwhisker(WSettings,FP.figureHandle,y,yLabel,yUnit,yscale{iScale},legendEntries);
            
            % save figure
            FP = FP.printFigure(figureName,figtxt,csv,figtxtTable);
        end
    end

    % save caption array for this sheet
    FP.saveCaptiontextArray;

end

return

    
function FP = doShadedArea(WSettings,Def,PopRunSet,FP)
    
    
    % get Reference simulation
if ~isempty(Def.ixPopRunSetRef)
    PKParameterRef = loadPKParameter({PopRunSet(Def.ixPopRunSetRef).name});
    OutputListRef = fullfile(load('tmp',PopRunSet(Def.ixPopRunSetRef).name,'outputList.mat'),'OutputList');
else
    PKParameterRef = [];
end

return


function x = collectPK(PopRunSet)

x = [];
for iSet = 1:length(PopRunSet)
    
    load(fullfile('tmp',PopRunSet(iSet).name,'outputList.mat'),'OutputList');
    PKPList = loadPKPList(PopRunSet(iSet).name);
    
    load(fullfile('tmp',PopRunSet(iSet).name,'outputList.mat'),'OutputList');

    load(fullfile('tmp',PopRunSet(iSet).name,'applicationProtocol.mat'),'PKParameterTemplate')

    if isempty(x)
        
        for iO = 1:length(PKPList); 
            
            x(iO).reportName = OutputList(iO).reportName;
            
            x(iO).pKParameterList = OutputList(iO).pKParameterList; %#ok<AGROW>
            
            % add diplayname
            [~,ix] = ismember(x(iO).pKParameterList(1,:),{PKParameterTemplate.name});
            x(iO).pKParameterList(4,:) = {PKParameterTemplate(ix).reportName}; %#ok<AGROW>
            
            for iPK = 1:length(PKPList{iO})
                x(iO).values(:,iSet,iPK) = PKPList{iO}(iPK).value.*OutputList(iO).pKParameterList{3,iPK};%#ok<AGROW>
            end
                
        end
    else
        for iO = 1:length(PKPList); 
            % get indices
            [jj,ixTarget] = ismember(OutputList(iO).pKParameterList(1,:),x(iO).pKParameterListOutputList(iO).pKParameterList(1,:));
            ixSource = [find(jj) find(~jj)];
            ixTarget = ixTarget(jj);
            
            % add new values to structure
            for iNew = find(~jj)
                ixTarget(end+1) =  size(x(iO).pKParameterList,2)+1;%#ok<AGROW>
                x(iO).pKParameterList(1:3,ixTarget(end)) = OutputList(iO).pKParameterList(:,iNew);%#ok<AGROW>
                
                % add displayname
                [~,ix] = ismember(x(iO).pKParameterList(1,iNew),{PKParameterTemplate.name});
                x(iO).pKParameterList{4,ixTarget(end)} = {PKParameterTemplate(ix).reportName}; %#ok<AGROW>

            end
            
            % add values
            for k = 1:length(ixTarget)
                x(iO).values(ixTarget(k),iSet,:) = percentiles(PKPList{iO}(ixSource(k)).value,WSettings.displayPercentiles)...
                    .*OutputList(iO).pKParameterList{3,ixTarget(k)};%#ok<AGROW>
            end
        end
    end
        
end

return

function PKPList = loadPKPList(simulationName)

fname = fullfile('tmp',simulationName,'pKPList.mat');

if exist(fname,'file')
    load(fname);
else
    PKPList = readPKParameter(simulationName);
end

return
