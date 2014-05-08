
function A = readFile()
%read .dat file or files into A
%   startPoint : the start point to be import (ms)
%   endPoint : the end point to be import (ms)
%   Example:
%   
%   Written by LiC    
%   05.08, 2014

[fileName,inputFilePath]=uigetfile('*.dat','Select the LabView Data','MultiSelect','on');

b = whos('fileName');
if strcmp(b.class, 'cell')
    numberFile = length(fileName);
    A = cell(1,numberFile);
     for k = 1:numberFile
         f = fullfile(inputFilePath,fileName{k});
         A{1,k} = importdata(f);
     end
else
     f = fullfile(inputFilePath,fileName);
     A = importdata(f);
end
