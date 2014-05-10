%values = {-1,2,2};
headers = {'First','Second','Third'};
file = 'myExample.xlsx';
xlswrite(file, headers)

% store data of file in variable a
[a,a1,a3] = xlsread(file);

% new data to write in file
t = zeros(1,10);

% last row with data in the file

nRows = (size(a3,1));

% plus 1 to write in the next line (if you have an header it should be + 2)
nRows = nRows +1;

% convert number to string
b = num2str(nRows);

% if you want to add data to the collum A you make concat strings
c = strcat('A', b);

% right to file the data t on the sheet Folha1 begining in the row c (e.g. A20)
xlswrite(file,t,1,c);