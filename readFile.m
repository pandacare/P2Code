
function [A,fileName,numberFile] = readFile()
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
    numberFile = 1;
    fileNameTemp = fileName;
    f = fullfile(inputFilePath,fileNameTemp);
    B= importdata(f);
    A = cell(1,1);
    A{1,1} = B;
    fileName = cell(1,1); 
    fileName{1,1} = fileNameTemp;
     
end
