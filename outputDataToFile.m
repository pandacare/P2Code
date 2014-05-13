function b = outputDataToFile(outputFileName,trialData)
%built .exel file for store result 
%   Example:
%   file = 'myExample.xlsx';
%   header = {'a','a','b'};
%   outputDataToFile(outputFileName);
%   Written by LiC    
%   05.08, 2014

file = outputFileName;
[a0,a1,a2] = xlsread(file);
nRows = (size(a2,1));
nRows = nRows +1;

% convert number to string
b = num2str(nRows);

% if you want to add data to the collum A you make concat strings
c = strcat('C', b);

% right to file the data t on the sheet Folha1 begining in the row c (e.g. A20)
xlswrite(file,trialData,1,c);