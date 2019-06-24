function generateCaptionTextxls
% GENERATECAPTIONTEXTXLS  generate caption text for report generation with word macro

% Open Systems Pharmacology Suite;  http://open-systems-pharmacology.org

figureMainDir = fullfile('figures');

figureDirs = dir(figureMainDir);
jj = cellfun(@(x) length(x)>2,{figureDirs.name});
figureDirs = figureDirs(jj);

for iDir = 1:length(figureDirs)
    
    iSection = 1;
    csvfile = dir(fullfile(figureMainDir,figureDirs(iDir).name,'S1__*.csv'));
    
    while ~isempty(csvfile)
        % read file
        fid = fopen(fullfile(figureMainDir,figureDirs(iDir).name,csvfile.name));
        tmp=textscan(fid,'%s','delimiter','\n');
        lines=tmp{1};
        fclose(fid);

        export = {};
        for iL = 1:length(lines)
            tmp = regexp(lines(iL),';','split');
            export(iL,1:length(tmp{1})) = tmp{1};
            
            if strcmp('pKParameter',figureDirs(iDir).name) 
                csvName = fullfile(figureMainDir,figureDirs(iDir).name,[export{iL,1} '.csv']);
                if exist(csvName,'file')
                    % read file
                    fid = fopen(csvName);
                    tmp=textscan(fid,'%s','delimiter','\n');
                    linesXls=tmp{1};
                    fclose(fid);
                    
                    tmp = regexp(linesXls(1),';','split');
                    export{iL,3} = tmp{1}{1};
                    
                    exportxls = {figureDirs(iDir).name};
                    for iLXls = 3:length(linesXls)
                        tmp = regexp(linesXls(iLXls),';','split');
                        tmp{1}{1} = [char(39) tmp{1}{1}];
                        tmp{1} = strrep(tmp{1},'\tau',char(964));
                        exportxls(iLXls-1,1:length(tmp{1})) = tmp{1};
                    end
                    xlswrite(fullfile(figureMainDir,figureDirs(iDir).name,[export{iL,1} '.xls']),exportxls,1);
                    
                end
            end
            
                
        end
        
        % repalce tau
        export(:,2) = strrep(export(:,2),'\tau',char(964));
        export(:,3) = strrep(export(:,3),'\tau',char(964));
        
        xlswrite(fullfile(figureMainDir,figureDirs(iDir).name,'Captiontext.xls'),export,sprintf('section_%d',iSection));
        

        iSection = 1 + iSection;
        csvfile = dir(fullfile(figureMainDir,figureDirs(iDir).name,sprintf('S%d__*.csv',iSection)));
        
        
    end
    
end
    

