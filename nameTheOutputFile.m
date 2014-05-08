function outputFileName = nameTheOutputFile()
%name the output file 
%   Example:
%   
%   Written by LiC    
%   05.08, 2014
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
currentTime = fix(clock);
currentTimeNameElement = arrayfun(@num2str, currentTime, 'UniformOutput', false);
currentTimeName = strcat(currentTimeNameElement(1),currentTimeNameElement(2),currentTimeNameElement(3),currentTimeNameElement(4),currentTimeNameElement(5),currentTimeNameElement(6));
prompt = {'Enter file name:','Enter condition:','Enter FM/P2'};
dlg_title = '';
num_lines = 1;
def = {currentTimeName{1,1},'300','P2'};
outputFileNameElement = inputdlg(prompt,dlg_title,num_lines,def);
outputFileName = strcat(outputFileNameElement{1},'_',outputFileNameElement{2},'_',outputFileNameElement{3},'.xlsx');