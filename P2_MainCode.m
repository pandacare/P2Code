%function
clc;
clear all;
close all;


%Sampling frequency 1000 Hz
F=1000; 
Fs=F/2;

%cuttoff frequency for grip force
w1=20;
WnPass=w1/Fs;

%channel option
channelForceXaxisThumb = 1; %Thumb shear force
channelForceXaxisIndex = 7; %Index shear force
channelForceZaxisThumb = channelForceXaxisThumb + 2; %Thumb grip
channelForceZaxisIndex = channelForceXaxisIndex + 2; %Index grip
channelEMGThumb = 13; % APB EMG
channelEMGIndex = 14; % FDI EMG

%Initialization
startPoint =1;                    %signal start reference ms
endPoint =45001;                  %signal stop reference ms
lengthTime=endPoint-startPoint;
dat=zeros(lengthTime,14);         %Matrix for holding raw data
datMod=zeros(lengthTime,14);      %Matrix for holding filtered data
meanStart = 11;
meanEnd   = 110; %for offset cancellation

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Read the raw data file %%%%%%%%%%%%%%%%%%%%

%read .dat file
[A,fileName,numberFile] = readFile();

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Output file setup %%%%%%%%%%%%%%%%%%%%%%%%%%

%Name the output .xsls file
[outputFileName, loadCondition] = nameTheOutputFile();

%Define parameter Header as the first row of the .xsls file
 global parameterArray;
 parameterArray        = {'InputFileName',...
                          'LoadCondition',...
                          'ThumbPerturbation',...
                          'ThumbGripForce',...
                          'ThumbEMGMethod0',...
                          'ThumbEMGMethod1',...
                          'ThumbEMGMethod2',...
                          'ThumbMaxGripForceValue',...
                          'ThumbMaxGripForceValueTimePoint',...
                          'ThumbMaxGripForceRate',...
                          'ThumbMaxGripForceRateTimePoint',...
                          'IndexPerturbation',...
                          'IndexGripForce',...
                          'IndexEMGMethod0',...
                          'IndexEMGMethod1',...
                          'IndexEMGMethod2',...
                          'IndexMaxGripForceValue',...
                          'IndexMaxGripForceValueTimePoint',...
                          'IndexMaxGripForceRate',...
                          'IndexMaxGripForceRateTimePoint'};
                      

xlswrite(outputFileName,parameterArray,1,'A1');
                      
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Data analysis %%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% Example of conducting data analysis with 'eventOnset' function:
%
%   %step 1: initial variables
%
%   inputData = zeros(1,100);
%
%   event = 'ThumbPerturbation'; % note: need to match with the element in
%                                        parameterArray
%
%   baseline = [0,0,0];          % note: baseline[3]:  baseLineStart, baseLineLength, nTimesSD
%
%   method = [0,0,0,0,0];        % note: method[5]:searchRangeStart, searchRangeLength,
%                                                  methodFunctionName:(singleCompare  0, meanCompare   1, meanWindowCompare 2)
%                                                  duration1(e.g.length for mean calculation in meanWindowMethod),
%                                                  duration2(e.g.length for comparison windown lengh in meanWindowMethod)
%   %step 2: call 'eventOnset' function
%
%   [eventValue, eventTime, eventID] = eventOnset(inputData,event,baseline,method);
%
%   %step 3: store result to result array
%
%   resultArray(eventID) = eventTime; % need to adjudt to eighter output time or value

for k = 1:numberFile
   resultArray           = zeros(1,length(parameterArray));
   inputFileName         = fileName{k};
   rawData               = A{1,k};
   %magnetRelease         = rawData(1,1);
   magnetRelease         = 20000;
   defaultBaselineLength = 100;
   searchRangeLength     = 1000;
   
   b1 = magnetRelease;
   b2 = defaultBaselineLength;
   m1 = b1 + 1;
   m2 = searchRangeLength;
   
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%Thumb%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   
   %ThumbPerturbation
   inputData = rawData((startPoint+1):endPoint,channelForceXaxisThumb);
   event = 'ThumbPerturbation'; 
   baseline = [b1-b2+1,b2,2]; %2 times SD for threshold calculation
   method = [m1,m2,0,5,0]; %singleCompare, duration 5 ms
   [eventValue, eventTime, eventID] = eventOnset(inputData,event,baseline,method);
   resultArray(eventID) = eventTime; 
        
   
    inputData = rawData((startPoint+1):endPoint,channelForceXaxisThumb);
   event = 'ThumbEMGMethod5'; 
   baseline = [b1-b2+1,b2,2]; %2 times SD for threshold calculation
   method = [m1,m2,0,5,0]; %singleCompare, duration 5 ms
   [eventValue, eventTime, eventID] = eventOnset(inputData,event,baseline,method);
   resultArray(eventID) = eventTime; 
   %ThumbGripForce
   %ThumbEMGMethod0
   %ThumbEMGMethod1                      
   %ThumbEMGMethod2
   %ThumbMaxGripForceValue
   %ThumbMaxGripForceValueTimePoint
   %ThumbMaxGripForceRate
   %ThumbMaxGripForceRateTimePoint
                           
                          
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%% Output data to file %%%%%%%%%%%%
  
   %Write result data to file
   currentRow = outputDataToFile(outputFileName,resultArray);
   %Write input file name to file
   xlswrite(outputFileName,{inputFileName},1,strcat('A', currentRow));
   %write load condition of this trial
   xlswrite(outputFileName,{loadCondition},1,strcat('B', currentRow));
   
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Display if only one file as input%%%%%%%%%%%%  




