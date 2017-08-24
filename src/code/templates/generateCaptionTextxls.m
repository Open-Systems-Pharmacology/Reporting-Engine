function generateCaptionTextxls
% GENERATECAPTIONTEXTXLS  generate caption text for report generation with word macro

% Open Systems Pharmacology Suite;  http://open-systems-pharmacology.org

figureDirs = dir('figures');
jj = cellfun(@(x) length(x)>2,{figureDirs.name});
figureDirs = figureDirs(jj);

for iDir = 1:length(figureDirs)
    
    iSection = 1;
    csvfile = dir(fullfile('figures',figureDirs(iDir).name,'S1__*.csv'));
    
    while ~isempty(csvfile)
        % read file
        fid = fopen(fullfile('figures',figureDirs(iDir).name,csvfile.name));
        tmp=textscan(fid,'%s','delimiter','\n');
        lines=tmp{1};
        fclose(fid);

        export = {};
        for iL = 1:length(lines)
            tmp = regexp(lines(iL),';','split');
            export(iL,:) = tmp{1};
        end
        
        xlswrite(fullfile('figures',figureDirs(iDir).name,'Captiontext.xls'),export,sprintf('section_%d',iSection));
        

        iSection = 1 + iSection;
        csvfile = dir(fullfile('figures',figureDirs(iDir).name,sprintf('S%d__*.csv',iSection)));
        
        
    end
    
end
    

