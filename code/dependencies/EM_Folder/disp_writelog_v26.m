function [] = disp_writelog_v26(pathOut, nameFile,textTowrite)
%% writelog_v26 
% write textTowrite in a file - file should exist before calling this
% function
%% E Mougel CIBM 02.2026

disp(textTowrite)
fileIDLog = fopen(fullfile(pathOut,nameFile),'a');
fprintf(fileIDLog,'%s\n',textTowrite);
fclose(fileIDLog);

end