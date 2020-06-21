function [fig_handle, DDIRatioQuali, GMFE, names] = plotSubunit(Subunit,DDIRatioGroups,Result,Observations,ObservedDataSets,PKParameter,WSettings,PlotSettings,AxesOptions)

DDISubunit ={};
for i=1:length(DDIRatioGroups)
    DDIRatios=DDIRatioGroups(i).DDIRatios;
    for j=1:length(DDIRatios)
        % Get the observed data
        for k=1:length(ObservedDataSets)
            if strcmp(DDIRatios(j).ObservedData, ObservedDataSets(k).Id)
                ObservedData = ObservedDataSets(k).y;
                break
            end
        end
        DDISubunit(size(DDISubunit,1)+1,1) = ObservedData{ObservedData.ID==DDIRatios(j).ObservedDataRecordId,{Subunit}};
        DDISubunit(size(DDISubunit,1),2) = num2cell(i);
        DDISubunit(size(DDISubunit,1),3) = num2cell(j);
        
    end
end
listSubunit = unique(DDISubunit(:,1));

for idxSub = 1:length(listSubunit)
    curSubunit = matlab.lang.makeValidName(listSubunit{idxSub});
    names.(curSubunit) = strrep(listSubunit{idxSub},'_',' ');
end

% Set the axis window
for k=1:length(PKParameter)
    for idxSub = 1:length(listSubunit)
        curSubunit = matlab.lang.makeValidName(listSubunit{idxSub});
        idSubunit = strcmp(DDISubunit(:,1),listSubunit{idxSub});
        subsetSubunit = DDISubunit(idSubunit,:);
        axisObsVsPred.(curSubunit)(k).min=NaN;
        axisObsVsPred.(curSubunit)(k).max=NaN;
        axisResVsObs.(curSubunit)(k).min=NaN;
        axisResVsObs.(curSubunit)(k).max=NaN;
        for jj = 1:size(subsetSubunit,1)
            i = cell2mat(subsetSubunit(jj,2));
            j = cell2mat(subsetSubunit(jj,3));
            axisObsVsPred.(curSubunit)(k).min=nanmin(nanmin(Result(i).RatioPK(k,j), Observations(i).RatioPK(k,j)),axisObsVsPred.(curSubunit)(k).min);
            axisObsVsPred.(curSubunit)(k).max=nanmax(nanmax(Result(i).RatioPK(k,j), Observations(i).RatioPK(k,j)),axisObsVsPred.(curSubunit)(k).max);
            axisResVsObs.(curSubunit)(k).min=nanmin(Result(i).RatioPK(k,j)./Observations(i).RatioPK(k,j),axisResVsObs.(curSubunit)(k).min);
            axisResVsObs.(curSubunit)(k).max=nanmax(Result(i).RatioPK(k,j)./Observations(i).RatioPK(k,j),axisResVsObs.(curSubunit)(k).max);
        end
    end
end

originalTitle = PlotSettings.title;
for idxSub = 1:length(listSubunit)
    curSubunit = matlab.lang.makeValidName(listSubunit{idxSub});
    idSubunit = strcmp(DDISubunit(:,1),listSubunit{idxSub});
    subsetSubunit = DDISubunit(idSubunit,:);
    PlotSettings.title = [originalTitle ' ' names.(curSubunit)];
    for ParameterIndex=1:length(PKParameter)
        % Initialize the Qualification Measures
        QualiMeasure.(curSubunit)(ParameterIndex).PointsTotal = 0;
        QualiMeasure.(curSubunit)(ParameterIndex).PointsGuest= 0;
        QualiMeasure.(curSubunit)(ParameterIndex).Points2fold= 0;
        
        % Initialize legend labels
        leg_labels={};
        
        % Get default axis limits for DDI Ratio plot and Guest equation
        if axisObsVsPred.(curSubunit)(ParameterIndex).min >= 1
            minAxisObsVsPred.(curSubunit)(ParameterIndex) = min([0.8*axisObsVsPred.(curSubunit)(ParameterIndex).min, 0.8]);
        else
            minAxisObsVsPred.(curSubunit)(ParameterIndex) = min([0.8*axisObsVsPred.(curSubunit)(ParameterIndex).min, 0.25]);
        end
        if axisObsVsPred.(curSubunit)(ParameterIndex).max >= 1
            maxAxisObsVsPred.(curSubunit)(ParameterIndex) = max([1.25*axisObsVsPred.(curSubunit)(ParameterIndex).max, 4]);
        else
            maxAxisObsVsPred.(curSubunit)(ParameterIndex) = max([1.25*axisObsVsPred.(curSubunit)(ParameterIndex).max, 1.25]);
        end
        
        if  axisResVsObs.(curSubunit)(ParameterIndex).min >= 1
            minAxisResVsObs.(curSubunit)(ParameterIndex) = 0.25;
        else
            minAxisResVsObs.(curSubunit)(ParameterIndex) = min([0.8*axisResVsObs.(curSubunit)(ParameterIndex).min, 0.25]);
        end
        if axisResVsObs.(curSubunit)(ParameterIndex).max >= 1
            maxAxisResVsObs.(curSubunit)(ParameterIndex) = max([1.25*axisResVsObs.(curSubunit)(ParameterIndex).max, 4]);
        else
            maxAxisResVsObs.(curSubunit)(ParameterIndex) = 4;
        end
        
        % GuestRatio vectors for plotting Guest et al. equation
        GuestRatio(ParameterIndex).x=10.^(log10(minAxisObsVsPred.(curSubunit)):0.001:log10(maxAxisObsVsPred.(curSubunit)));
        
        GuestRatio(ParameterIndex).x = 10.^(log10(minAxisObsVsPred.(curSubunit)(ParameterIndex)):0.001:log10(maxAxisObsVsPred.(curSubunit)(ParameterIndex)));
        [GuestRatio(ParameterIndex).yup, GuestRatio(ParameterIndex).ylo] = DDIRatioGuestEquation(GuestRatio(ParameterIndex).x);
        
        % Create figure for Obs vs Pred with lines from Guest equation
        [ax, fig_handle.(curSubunit)(ParameterIndex).predictedVsObserved] = getReportFigureQP(WSettings,1,1,[],PlotSettings);
        setFigureOptions(AxesOptions.DDIRatioPlotsPredictedVsObserved);
        plot(GuestRatio(ParameterIndex).x, GuestRatio(ParameterIndex).x, '-k', 'Linewidth', 1.5, 'HandleVisibility','off');
        plot(GuestRatio(ParameterIndex).x, GuestRatio(ParameterIndex).x/2, ':k', 'Linewidth', 1, 'HandleVisibility','off');
        plot(GuestRatio(ParameterIndex).x, GuestRatio(ParameterIndex).x*2, ':k', 'Linewidth', 1, 'HandleVisibility','off');
        plot(GuestRatio(ParameterIndex).x, GuestRatio(ParameterIndex).yup, '-k', 'Linewidth', 1, 'HandleVisibility','off');
        plot(GuestRatio(ParameterIndex).x, GuestRatio(ParameterIndex).ylo, '-k', 'Linewidth', 1, 'HandleVisibility','off');
        
        xlabel(sprintf('Observed %s Ratio', PKParameter{ParameterIndex})); ylabel(sprintf('Predicted %s Ratio', PKParameter{ParameterIndex}));
        axis([minAxisObsVsPred.(curSubunit)(ParameterIndex) maxAxisObsVsPred.(curSubunit)(ParameterIndex) minAxisObsVsPred.(curSubunit)(ParameterIndex) maxAxisObsVsPred.(curSubunit)(ParameterIndex)]);
        
        % Create figure for Residuals vs Obs with lines from Guest equation
        [ax, fig_handle.(curSubunit)(ParameterIndex).residualsVsObserved] = getReportFigureQP(WSettings,1,1,[],PlotSettings);
        setFigureOptions(AxesOptions.DDIRatioPlotsResidualsVsObserved);
        plot(GuestRatio(ParameterIndex).x, ones(size(GuestRatio(ParameterIndex).x)), '-k', 'Linewidth', 1.5, 'HandleVisibility','off');
        plot(GuestRatio(ParameterIndex).x, ones(size(GuestRatio(ParameterIndex).x))/2, ':k', 'Linewidth', 1, 'HandleVisibility','off');
        plot(GuestRatio(ParameterIndex).x, ones(size(GuestRatio(ParameterIndex).x))*2, ':k', 'Linewidth', 1, 'HandleVisibility','off');
        plot(GuestRatio(ParameterIndex).x, GuestRatio(ParameterIndex).yup./GuestRatio(ParameterIndex).x, '-k', 'Linewidth', 1, 'HandleVisibility','off');
        plot(GuestRatio(ParameterIndex).x, GuestRatio(ParameterIndex).ylo./GuestRatio(ParameterIndex).x, '-k', 'Linewidth', 1, 'HandleVisibility','off');
        
        xlabel(sprintf('Observed %s Ratio', PKParameter{ParameterIndex})); ylabel(sprintf('Predicted %s Ratio / Observed %s Ratio', PKParameter{ParameterIndex}, PKParameter{ParameterIndex}));
        axis([minAxisObsVsPred.(curSubunit)(ParameterIndex) maxAxisObsVsPred.(curSubunit)(ParameterIndex) minAxisResVsObs.(curSubunit)(ParameterIndex) maxAxisResVsObs.(curSubunit)(ParameterIndex)]);
        
        listGroups = unique([subsetSubunit{:,2}]);
        for ii = 1:length(listGroups)
            i = listGroups(ii);
            jj = [subsetSubunit{[subsetSubunit{:,2}]==i,3}];
            
            % Update the number of points for each group
            conditionPoints = ~isnan(Observations(i).RatioPK(ParameterIndex,jj));
            QualiMeasure.(curSubunit)(ParameterIndex).PointsTotal = QualiMeasure.(curSubunit)(ParameterIndex).PointsTotal + length(Observations(i).RatioPK(ParameterIndex,conditionPoints));
            [UpperBound, LowerBound] = DDIRatioGuestEquation(Observations(i).RatioPK(ParameterIndex,jj));
            conditionGuest = (Result(i).RatioPK(ParameterIndex,jj) >= LowerBound & Result(i).RatioPK(ParameterIndex,jj) <= UpperBound);
            QualiMeasure.(curSubunit)(ParameterIndex).PointsGuest= QualiMeasure.(curSubunit)(ParameterIndex).PointsGuest + length(Result(i).RatioPK(ParameterIndex,conditionGuest));
            condition2fold = (Result(i).RatioPK(ParameterIndex,jj) >= Observations(i).RatioPK(ParameterIndex,jj)/2 & Result(i).RatioPK(ParameterIndex,jj) <= Observations(i).RatioPK(ParameterIndex,jj)*2);
            QualiMeasure.(curSubunit)(ParameterIndex).Points2fold= QualiMeasure.(curSubunit)(ParameterIndex).Points2fold+ length(Result(i).RatioPK(ParameterIndex,condition2fold));
            
            % Plot part
            set(0, 'CurrentFigure', fig_handle.(curSubunit)(ParameterIndex).predictedVsObserved);
            pp=plot(Observations(i).RatioPK(ParameterIndex,jj), Result(i).RatioPK(ParameterIndex,jj), 'o', 'Linewidth',1);
            setCurveOptions(pp, DDIRatioGroups(i));
            
            set(0, 'CurrentFigure', fig_handle.(curSubunit)(ParameterIndex).residualsVsObserved);
            pp=plot(Observations(i).RatioPK(ParameterIndex,jj), Result(i).RatioPK(ParameterIndex,:)./Observations(i).RatioPK(ParameterIndex,jj), 'o', 'Linewidth',1);
            setCurveOptions(pp, DDIRatioGroups(i));
            
            if isfield(DDIRatioGroups(i), 'Caption')
                leg_labels=[leg_labels DDIRatioGroups(i).Caption];
            end
  
        end
        set(0, 'CurrentFigure', fig_handle.(curSubunit)(ParameterIndex).predictedVsObserved);
        %legend(leg_labels, 'Location', 'northoutside');
        addLegendWithConstantPlotArea(gcf,gca,leg_labels);
        
        set(0, 'CurrentFigure', fig_handle.(curSubunit)(ParameterIndex).residualsVsObserved);
        %legend(leg_labels, 'Location', 'northoutside');
        addLegendWithConstantPlotArea(gcf,gca,leg_labels);
    end
end

% Calculation of GMFE
for idxSub = 1:length(listSubunit)
    curSubunit = matlab.lang.makeValidName(listSubunit{idxSub});
    idSubunit = strcmp(DDISubunit(:,1),listSubunit{idxSub});
    subsetSubunit = DDISubunit(idSubunit,:);
    RatioPK.(curSubunit)=[];
    listGroups = unique([subsetSubunit{:,2}]);
    for ii = 1:length(listGroups)
        i = listGroups(ii);
        jj = [subsetSubunit{[subsetSubunit{:,2}]==i,3}];
        RatioPK.(curSubunit) = [RatioPK.(curSubunit) Result(i).RatioPK(:,jj)./Observations(i).RatioPK(:,jj)];
        
    end
    RatioPK.(curSubunit) = reshape(RatioPK.(curSubunit)', [], length(PKParameter));
    for k=1:length(PKParameter)
        GMFE.(curSubunit)(k) = 10.^(sum(abs(log10(RatioPK.(curSubunit)(~isnan(RatioPK.(curSubunit)(:,k)),k))))./length(RatioPK.(curSubunit)(~isnan(RatioPK.(curSubunit)(:,k)),k)));
        fprintf('%s\t%s: GMFE = %f \n', curSubunit,PKParameter{k}, GMFE.(curSubunit)(k));
    end
end

% Get the DDI Ratio Qualification
for idxSub = 1:length(listSubunit)
    curSubunit = matlab.lang.makeValidName(listSubunit{idxSub});
    for k=1:length(PKParameter)
        DDIRatioQualiHeader = {PKParameter{k}, 'Number', 'Ratio [%]'};
        DDIRatioQuali_1st_Column = {'Points total'; 'Points within Guest et al.'; 'Points within 2-fold'};
        DDIRatioQuali_Other_Columns = num2cell([QualiMeasure.(curSubunit)(k).PointsTotal NaN; ...
            QualiMeasure.(curSubunit)(k).PointsGuest 100.*QualiMeasure.(curSubunit)(k).PointsGuest./QualiMeasure.(curSubunit)(k).PointsTotal; ...
            QualiMeasure.(curSubunit)(k).Points2fold 100.*QualiMeasure.(curSubunit)(k).Points2fold./QualiMeasure.(curSubunit)(k).PointsTotal]);
        DDIRatioQuali_Other_Columns{1,2}='-';
        
        DDIRatioQuali.(curSubunit)(k).Output = [DDIRatioQualiHeader  ; ...
            DDIRatioQuali_1st_Column DDIRatioQuali_Other_Columns];
        
        fprintf('%s: Qualification Measures for PK parameter %s \n', curSubunit, PKParameter{k});
        disp(DDIRatioQuali.(curSubunit)(k).Output);
        
    end
end

% % Get the DDI Ratio Qualification
% DDIRatioQualiHeader = {Subunit,'PK Parameter',...
%             'Points total (n)', 'Points within Guest et al. (n)','Points within Guest et al. (%)',...
%             'Points within 2-fold (n)','Points within 2-fold (%)'};
% DDIRatioQuali_1st_Columns = {};
% DDIRatioQuali_Other_Columns = {};
% for idxSub = 1:length(listSubunit)
%     curSubunit = listSubunit{idxSub};
%     for k=1:length(PKParameter)
%         
%         DDIRatioQuali_1st_Columns = [DDIRatioQuali_1st_Columns;{curSubunit,PKParameter{k}}];
%         DDIRatioQuali_Other_Columns = [DDIRatioQuali_Other_Columns;num2cell([QualiMeasure.(curSubunit)(k).PointsTotal, ...
%             QualiMeasure.(curSubunit)(k).PointsGuest 100.*QualiMeasure.(curSubunit)(k).PointsGuest./QualiMeasure.(curSubunit)(k).PointsTotal, ...
%             QualiMeasure.(curSubunit)(k).Points2fold 100.*QualiMeasure.(curSubunit)(k).Points2fold./QualiMeasure.(curSubunit)(k).PointsTotal])];
%         %DDIRatioQuali_Other_Columns{1,2}='-';
% 
%         %fprintf('Qualification Measures for PK parameter %s \n', PKParameter{k});
%         %disp(DDIRatioQuali(k).Output);
%         
%     end
% end
% DDIRatioQuali = [DDIRatioQualiHeader  ; ...
%             DDIRatioQuali_1st_Columns DDIRatioQuali_Other_Columns];


function [UpperLimit, LowerLimit] = DDIRatioGuestEquation(Robs, delta)
% Upper and Lower limits for DDI Ratio plots as proposed by Guest et al.


% Delta might be different from 1 to use different limits
if ~exist('delta')
    delta=1;
end

symRobs=Robs;

symRobs(Robs<1)=1./symRobs(symRobs<1);

Limit = (delta + 2.*(symRobs-1))./symRobs;

symUpperLimit = symRobs.*Limit;
symLowerLimit = symRobs./Limit;

UpperLimit(isnan(Robs))=NaN;
LowerLimit(isnan(Robs))=NaN;

UpperLimit(Robs>=1)=symUpperLimit(Robs>=1);
LowerLimit(Robs>=1)=symLowerLimit(Robs>=1);
UpperLimit(Robs<1)=1./symLowerLimit(Robs<1);
LowerLimit(Robs<1)=1./symUpperLimit(Robs<1);