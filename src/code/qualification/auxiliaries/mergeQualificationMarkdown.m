function mergeQualificationMarkdown(inputFolder, outputMarkdownFile)
    % Function to merge all markdown files within a folder(including all subfolders) into one markdown file
    % inputFolder - Complete path to the folder where MD files are saved:
    % outputMarkdownFile - Complete path to the folder and file name of the new output markdown.

    % Generate the list of all the md files with the folder structure:
    [MD_FileNames] = func_folder_N_MD(inputFolder);

    % Remove the summary md file
    MD_FileNames_wo_summary = MD_FileNames(~strcmp(MD_FileNames,'summary.md'));

    %% Report specific word and syntax:
    ReportSpecific_word    = 'final input parameters';
    ReportSpecific_hashtag = '#';

    FontSizeSyntaxOpen = '<h1>';
    FontSizeSyntaxClose = '</h1>';
    PageBreak = '<div style="page-break-after: always;"></div>';

    %% Merge the content of the various md into cell array
    MD_FinalContent = {};

    % Read the content of '_intro.md'
    Fullfilename_Src = fullfile(inputFolder, MD_FileNames_wo_summary{1});
    fid = fopen(Fullfilename_Src,'r','n','UTF-8');
    tmp = textscan(fid,'%s','delimiter',newline,'whitespace','');
    fclose(fid);

    % Remove the hash from the title, change the font size of the title and add table of content:
    modtmp = tmp{1};
    flag_Hash = contains(modtmp,ReportSpecific_hashtag);
    tmp_str = modtmp{flag_Hash};
    Index_Hash = strfind(tmp_str,ReportSpecific_hashtag);
    tmp_str2 =  strcat(FontSizeSyntaxOpen, tmp_str((Index_Hash+1):end), FontSizeSyntaxClose);
    modtmp(flag_Hash) = {tmp_str2};

    %
    MD_FinalContent = vertcat(MD_FinalContent,modtmp);

    % Add a table of content
    MD_FinalContent = vertcat(MD_FinalContent,PageBreak);
    MD_FinalContent = vertcat(MD_FinalContent,'[toc]');

    %% Read other md files and merge into the cell array
    for i = 2:length(MD_FileNames_wo_summary)
        Fullfilename_Src = fullfile(inputFolder, MD_FileNames_wo_summary{i});

        fid = fopen(Fullfilename_Src,'r','n','UTF-8');
        tmp = textscan(fid,'%s','delimiter',newline,'whitespace','');
        fclose(fid);
        %%
        modtmp = tmp{1};
        if contains(modtmp{1},ReportSpecific_word,'IgnoreCase',true)
            % Check 
            NrHashInFirstLine = length(strfind(modtmp{1},ReportSpecific_hashtag));

            flag_LinesWithHash = contains(modtmp,ReportSpecific_hashtag);
            flag_LinesWithHash_NotFirstLine = flag_LinesWithHash;
            flag_LinesWithHash_NotFirstLine(1) = false;

            LinesWithHash_NotFirstLine = modtmp(flag_LinesWithHash_NotFirstLine); 

            modtmp(flag_LinesWithHash_NotFirstLine) = cellfun(@(X)(strcat(repmat(ReportSpecific_hashtag,1,NrHashInFirstLine),X)),...
                                                                   LinesWithHash_NotFirstLine,...
                                                                   'UniformOutput',false);
        end
        %%
        MD_FinalContent = vertcat(MD_FinalContent,modtmp);
    end

    %% Write the final content into new MD file
    fid = fopen(outputMarkdownFile,'w','n','UTF-8');
    for j = 1:length(MD_FinalContent)
        fprintf(fid,'%s',MD_FinalContent{j}); % Experiment ID
        fprintf(fid,'%s',newline);
    end

    fclose(fid);
end

function [MD_FileNames,FolderNames] = func_folder_N_MD(CompleteFolderName) 
    %% Find the list of all MD files recursively within a folder:
    
    %% Folder names
    AllFolderNFiles = dir(CompleteFolderName);
    
    %% Select only folders:
    flagDot = strcmp({AllFolderNFiles.name},'.');
    flagDot2 = strcmp({AllFolderNFiles.name},'..');
    flagMd = contains({AllFolderNFiles.name},'.md');
    flagPng = contains({AllFolderNFiles.name},'.png');
    
    % MD file
    fileMD = AllFolderNFiles(flagMd);
    folderNameOnly = AllFolderNFiles((~flagDot) & (~flagDot2) & ~(flagMd) & ~(flagPng));
    
    %% MD file in this folder 
    MD_FileNames = {fileMD.name};
    
    %% Folders within this folder:
    FolderNames  = {folderNameOnly.name};
    
    %%
    if isempty(FolderNames)
        return % If there is no folder, finish the function
    else
        % If there is a folder(s), run this function again and update the filenames list
        for i = 1:length(FolderNames)
            [MD_FileNames_new] = func_folder_N_MD([CompleteFolderName filesep FolderNames{i}]);
            MD_FileNames_new_withFolderName = fullfile(FolderNames{i},MD_FileNames_new);
            MD_FileNames = horzcat(MD_FileNames,MD_FileNames_new_withFolderName);
        end
    end
end
