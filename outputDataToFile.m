%function outputDataToFile(outputFileName,header)
%built .exel file for store result 
%   Example:
%   
%   Written by LiC    
%   05.08, 2014

%file = outputFileName;
file = 'myExample.xlsx';
header = {'a','a','b'};
xlswrite(file,header);
a = xlsread(file);
x = size(a,1);
if ~size(a,1)
    %we build header into the file
    xlswrite(file, header);
end

[a0,a1,a2] = xlsread(file);
nRows = (size(a2,1));

nRows = nRows +1;

t = zeros(1,10);

% convert number to string
b = num2str(nRows);

% if you want to add data to the collum A you make concat strings
c = strcat('A', b);

% right to file the data t on the sheet Folha1 begining in the row c (e.g. A20)
xlswrite(file,t,1,c);