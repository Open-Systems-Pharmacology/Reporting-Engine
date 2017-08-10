classdef ReportFigurePrint
%REPORTFIGUREPRINT class which manages the print of figures

% Open Systems Pharmacology Suite;  http://open-systems-pharmacology.org

    
    
    properties
        % formats
        EMF = false;
        FIG = false;
        CSV = true;
        PNG = false;
        
        % Ps Format
        PsOrientation='portrait';
        PaperPosition=[0 0 20 27];

        % filing directory for figures
        figureDir='';
        
        % figure handle
        figureHandle=nan;
        
        % Position of figure
        Position = [];
        
        % Structure to manage captiontext xls, used for word macro
        CaptiontextFigtextArray = {'header'}
        CaptiontextSheet = 'sheet';
        CaptiontextSheetNo = 0;
    end
    
    methods
        function obj=ReportFigurePrint(resultDir,printFormatList,format)
            
            % set figur directory and create it
            obj.figureDir=resultDir;
            if ~exist(obj.figureDir,'dir')
                mkdir(obj.figureDir);
            else
                delete(fullfile(obj.figureDir,'*'));
            end
            
            % set print format
            for iP = 1:length(printFormatList)
                obj.(printFormatList{iP}) = true;
            end
            

            if exist('format','var')
                switch upper(format)
                    case 'PORTRAIT'
                        obj.PsOrientation='Portrait';
                        obj.PaperPosition=[2 2 20 27];
                    case 'SQUARE'
                        obj.PsOrientation='Portrait';
                        obj.PaperPosition=[2 2 20 20];
                    case 'LANDSCAPE'
                        obj.PsOrientation='Landscape';
                        obj.PaperPosition=[0 0 27 20];
                    otherwise
                        error('unknonw format')
                end
            end

            
            obj.figureHandle=figure;
            
        end
        
        
        function obj=printFigure(obj,figureName,figtxt,csv,figtxt_csv)
            % prints figure in selected format
            
            figureName = sprintf('S%d_%s',obj.CaptiontextSheetNo,figureName);
            
            if isa(obj.figureHandle,'matlab.ui.Figure')
                iFig = ['-f' num2str(obj.figureHandle.Number)];
            else
                iFig = ['-f' num2str(obj.figureHandle)];
            end
            
                        
            if ~isempty(obj.Position)
                set(obj.figureHandle,'position',obj.Position);
                pause(0.1);
            end

            % change figure format
            set(obj.figureHandle,'PaperType','a4');
            set(obj.figureHandle, 'PaperPositionMode', 'manual');
            set(obj.figureHandle, 'PaperOrientation',obj.PsOrientation);
            old_PaperPosition=get(obj.figureHandle, 'PaperPosition');
            set(obj.figureHandle, 'PaperUnits', 'centimeters');
            set(obj.figureHandle, 'PaperPosition', obj.PaperPosition);
           
            
            % emf
            if  obj.EMF
                print('-dmeta', iFig, fullfile(obj.figureDir, [figureName '.emf']));
            end
            
            %fig
            if  obj.FIG
                saveas(obj.figureHandle,fullfile(obj.figureDir, [figureName '.fig']))
            end
            
            if  obj.PNG
                saveas(obj.figureHandle,fullfile(obj.figureDir, figureName),'png');
            end
            
            % reset figure format
            set(obj.figureHandle, 'PaperPosition', old_PaperPosition);
 
            % add to captiontext    
            obj.CaptiontextFigtextArray{end+1,1}=figureName;
            obj.CaptiontextFigtextArray{end,2}=figtxt;
 
            % write figure as table
            if exist('csv','var') && ~isempty(csv)
                
                % add header to table
                csv = [repmat({''},2,size(csv,2));csv];
                csv{1,1} = strrep(figtxt_csv,';',',');
                writeTabCellArray(csv,fullfile(obj.figureDir, [figureName '.csv']));
            end
        end
        
        function obj = iniCaptiontextFigtextArray(obj,outline,sheet)
            % Initialize Structure to manage captiontext xls, used for word macro
            
            obj.CaptiontextFigtextArray = {outline};            
            obj.CaptiontextSheetNo = obj.CaptiontextSheetNo+1;
            obj.CaptiontextSheet = sprintf('S%d__%s.csv',obj.CaptiontextSheetNo,sheet);
            
        end
        
        
        
        function saveCaptiontextArray(obj)
            
            % Save the Caption text  used for word macro
            writeTabCellArray(obj.CaptiontextFigtextArray,fullfile(obj.figureDir,obj.CaptiontextSheet));

        end
    end
    
end

