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
endPoint =44001;                  %signal stop reference ms
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
                          'ThumbMaxGripForceValue',...
                          'ThumbMaxGripForceValueTimePoint',...
                          'ThumbGripForceRate',...
                          'ThumbMaxGripForceRate',...
                          'ThumbMaxGripForceRateTimePoint',...
                          'ThumbEMGMethod0',...
                          'ThumbEMGMethod1',...
                          'ThumbEMGMethod2',...
                          'IndexPerturbation',...
                          'IndexGripForce',...
                          'IndexMaxGripForceValue',...
                          'IndexMaxGripForceValueTimePoint',...
                          'IndexGripForceRate',...
                          'IndexMaxGripForceRate',...
                          'IndexMaxGripForceRateTimePoint',...
                          'IndexEMGMethod0',...
                          'IndexEMGMethod1',...
                          'IndexEMGMethod2'};
                      

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
   bufferData            = zeros((endPoint-startPoint),1);
   magnetRelease         = round(rawData(1,1));
   displayWindowLength   = 2000;
   defaultBaselineLength = 100;
   searchRangeLength     = 1000;
   
   displayTimeStart = magnetRelease-500;
   displayTimeEnd   = displayTimeStart + 1500;
   
   
   b1 = magnetRelease;
   b2 = defaultBaselineLength;
   m1 = b1;
   m2 = searchRangeLength;
   
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%Thumb%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   
   %ThumbPerturbation
   %Raw data 
   inputData = rawData((startPoint+1):endPoint,channelForceXaxisThumb);
   event = 'ThumbPerturbation'; 
   baseline = [b1-b2+1,b2,2]; %2 times SD for threshold calculation
   method = [m1,m2,0,5,0]; %singleCompare, duration 5 ms
   [eventValue, eventTime, eventID] = eventOnset(inputData,event,baseline,method);
   resultArray(eventID) = eventTime; 
  
   %ThumbGripForce
   %Filtered data 
   [b,a] = butter(4,WnPass,'low');
   inputData =abs(filter(b,a,rawData((startPoint+1):endPoint,channelForceZaxisThumb)));
   event = 'ThumbGripForce'; 
   baseline = [b1-b2+1,b2,3]; %3 times SD for threshold calculation
   method = [m1,m2,0,100,0]; %singleCompare, duration 100 ms
   [eventValue, eventTime, eventID] = eventOnset(inputData,event,baseline,method);
   resultArray(eventID) = eventTime; 
   
   %ThumbMaxGripForceValue
   inputData =abs(filter(b,a,rawData((startPoint+1):endPoint,channelForceZaxisThumb)));
   event = 'ThumbMaxGripForceValue';
   baseline = [b1-b2+1,b2,3]; %3 times SD for threshold calculation
   method = [m1,m2,3,100,0];
   [eventValue, eventTime, eventID] = eventOnset(inputData,event,baseline,method);
   resultArray(eventID) = eventValue; 
   
   %ThumbMaxGripForceValueTimePoint
    inputData =abs(filter(b,a,rawData((startPoint+1):endPoint,channelForceZaxisThumb)));
   event = 'ThumbMaxGripForceValueTimePoint';
   baseline = [b1-b2+1,b2,3]; %3 times SD for threshold calculation
   method = [m1,m2,3,100,0];
   [eventValue, eventTime, eventID] = eventOnset(inputData,event,baseline,method);
   resultArray(eventID) = eventTime; 
   
   %ThumbGripForceRate
   h = waitbar(0, 'waiting ...');
   inputData =abs(filter(b,a,rawData((startPoint+1):endPoint,channelForceZaxisThumb)));
   tempData = inputData;
   GripForceDuration = 30;
   GripForceRateData = zeros(1,(endPoint-startPoint));
    for j = (displayTimeStart+1):displayTimeEnd
        a1 = 0;
        a2 = 0;
        a3 = 0;
        a4 = 0;
        for i = j:(j+GripForceDuration-1)
            a1 = tempData(i)*((i-j+1)/F)+a1;
            a2 = tempData(i)+a2;
            a3 = ((i-j+1)/F)+a3;
            a4 = ((i-j+1)/F)^2+a4;
        end
        GripForceRateData(j) = (GripForceDuration*a1-a2*a3)/(GripForceDuration*a4-a3*a3);
        waitbar((j-(displayTimeStart+1))/(displayTimeEnd-(displayTimeStart+1)),h);
    end
   delete(h);
   inputData = abs(GripForceRateData);
   event = 'ThumbGripForceRate'; 
   baseline = [b1-b2+1,b2,3]; %2 times SD for threshold calculation
   method = [m1,m2,0,50,0]; %singleCompare, duration 50 ms
   [eventValue, eventTime, eventID] = eventOnset(inputData,event,baseline,method);
   resultArray(eventID) = eventTime;
    
   %ThumbMaxGripForceRate
   inputData = abs(GripForceRateData);
   event = 'ThumbMaxGripForceRate'; 
   baseline = [b1-b2+1,b2,3]; %2 times SD for threshold calculation
   method = [m1,m2,3,0,0]; %max
   [eventValue, eventTime, eventID] = eventOnset(inputData,event,baseline,method);
   resultArray(eventID) = eventValue; 
   
   %ThumbMaxGripForceRateTimePoint
   inputData = abs(GripForceRateData);
   event = 'ThumbMaxGripForceRateTimePoint'; 
   baseline = [b1-b2+1,b2,3]; %2 times SD for threshold calculation
   method = [m1,m2,3,0,0]; %max: method 3
   [eventValue, eventTime, eventID] = eventOnset(inputData,event,baseline,method);
   resultArray(eventID) = eventTime; 
   
   %%%%%ThumbEMG%%%%%%%
   [b50Hz,a50Hz]=butter(4,50/Fs,'low');
   [b5Hz,a5Hz]=butter(4,5/Fs,'low');
   averageMeanThumb = mean(rawData(meanStart:meanEnd,channelEMGThumb));
   tempData = rawData(:,channelEMGThumb)- averageMeanThumb;
   filtered50HzData = filtfilt(b50Hz,a50Hz,tempData);
   filtered5HzData = filtfilt(b5Hz,a5Hz,tempData);
   inputData = abs(filtered50HzData);%used as inputData for following
   
   %ThumbEMGMethod0
   event = 'ThumbEMGMethod0'; 
   baseline = [b1-b2+1,b2,3]; %3 times SD for threshold calculation
   method = [m1,m2,0,50,0]; %single compare
   [eventValue, eventTime, eventID] = eventOnset(inputData,event,baseline,method);
   resultArray(eventID) = eventTime; 
   
   %ThumbEMGMethod1 
   event = 'ThumbEMGMethod1'; 
   baseline = [b1-b2+1,b2,3]; %2 times SD for threshold calculation
   method = [m1,m2,1,50,0]; %mean compare
   [eventValue, eventTime, eventID] = eventOnset(inputData,event,baseline,method);
   resultArray(eventID) = eventTime; 
   
   %ThumbEMGMethod2
   event = 'ThumbEMGMethod2'; 
   baseline = [b1-b2+1,b2,3]; %2 times SD for threshold calculation
   method = [m1,m2,2,0,0]; %mean window compare
   [eventValue, eventTime, eventID] = eventOnset(inputData,event,baseline,method);
   resultArray(eventID) = eventTime; 
   
   %%%%%%%%%%%%%%%%%%%%%%% Index %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %IndexPerturbation
   %Raw data 
   inputData = rawData((startPoint+1):endPoint,channelForceXaxisIndex);
   event = 'IndexPerturbation'; 
   baseline = [b1-b2+1,b2,2]; %2 times SD for threshold calculation
   method = [m1,m2,0,5,0]; %singleCompare, duration 5 ms
   [eventValue, eventTime, eventID] = eventOnset(inputData,event,baseline,method);
   resultArray(eventID) = eventTime; 
  
   %IndexGripForce
   %Filtered data 
   [b,a] = butter(4,WnPass,'low');
   inputData =abs(filter(b,a,rawData((startPoint+1):endPoint,channelForceZaxisIndex)));
   event = 'IndexGripForce'; 
   baseline = [b1-b2+1,b2,3]; %3 times SD for threshold calculation
   method = [m1,m2,0,100,0]; %singleCompare, duration 100 ms
   [eventValue, eventTime, eventID] = eventOnset(inputData,event,baseline,method);
   resultArray(eventID) = eventTime; 
   
   %IndexMaxGripForceValue
   inputData =abs(filter(b,a,rawData((startPoint+1):endPoint,channelForceZaxisIndex)));
   event = 'IndexMaxGripForceValue';
   baseline = [b1-b2+1,b2,3]; %3 times SD for threshold calculation
   method = [m1,m2,3,100,0];
   [eventValue, eventTime, eventID] = eventOnset(inputData,event,baseline,method);
   resultArray(eventID) = eventValue; 
   
   %IndexMaxGripForceValueTimePoint
    inputData =abs(filter(b,a,rawData((startPoint+1):endPoint,channelForceZaxisIndex)));
   event = 'IndexMaxGripForceValueTimePoint';
   baseline = [b1-b2+1,b2,3]; %3 times SD for threshold calculation
   method = [m1,m2,3,100,0];
   [eventValue, eventTime, eventID] = eventOnset(inputData,event,baseline,method);
   resultArray(eventID) = eventTime; 
   
   %IndexGripForceRate
   h = waitbar(0, 'waiting ...');
   inputData =abs(filter(b,a,rawData((startPoint+1):endPoint,channelForceZaxisIndex)));
   tempData = inputData;
   GripForceDuration = 30;
   GripForceRateData = zeros(1,(endPoint-startPoint));
    for j = (displayTimeStart+1):displayTimeEnd
        a1 = 0;
        a2 = 0;
        a3 = 0;
        a4 = 0;
        for i = j:(j+GripForceDuration-1)
            a1 = tempData(i)*((i-j+1)/F)+a1;
            a2 = tempData(i)+a2;
            a3 = ((i-j+1)/F)+a3;
            a4 = ((i-j+1)/F)^2+a4;
        end
        GripForceRateData(j) = (GripForceDuration*a1-a2*a3)/(GripForceDuration*a4-a3*a3);
        waitbar((j-(displayTimeStart+1))/(displayTimeEnd-(displayTimeStart+1)),h);
    end
   delete(h);
   inputData = abs(GripForceRateData);
   event = 'IndexGripForceRate'; 
   baseline = [b1-b2+1,b2,3]; %2 times SD for threshold calculation
   method = [m1,m2,0,50,0]; %singleCompare, duration 50 ms
   [eventValue, eventTime, eventID] = eventOnset(inputData,event,baseline,method);
   resultArray(eventID) = eventTime;
    
   %IndexMaxGripForceRate
   inputData = abs(GripForceRateData);
   event = 'IndexMaxGripForceRate'; 
   baseline = [b1-b2+1,b2,3]; %2 times SD for threshold calculation
   method = [m1,m2,3,0,0]; %max
   [eventValue, eventTime, eventID] = eventOnset(inputData,event,baseline,method);
   resultArray(eventID) = eventValue; 
   
   %IndexMaxGripForceRateTimePoint
   inputData = abs(GripForceRateData);
   event = 'IndexMaxGripForceRateTimePoint'; 
   baseline = [b1-b2+1,b2,3]; %2 times SD for threshold calculation
   method = [m1,m2,3,0,0]; %max: method 3
   [eventValue, eventTime, eventID] = eventOnset(inputData,event,baseline,method);
   resultArray(eventID) = eventTime; 
   
   %%%%%IndexEMG%%%%%%%
   [b50Hz,a50Hz]=butter(4,50/Fs,'low');
   [b5Hz,a5Hz]=butter(4,5/Fs,'low');
   averageMeanIndex = mean(rawData(meanStart:meanEnd,channelEMGIndex));
   tempData = rawData(:,channelEMGIndex)- averageMeanIndex;
   filtered50HzData = filtfilt(b50Hz,a50Hz,tempData);
   filtered5HzData = filtfilt(b5Hz,a5Hz,tempData);
   inputData = abs(filtered50HzData);%used as inputData for following
   
   %IndexEMGMethod0
   event = 'IndexEMGMethod0'; 
   baseline = [b1-b2+1,b2,3]; %3 times SD for threshold calculation
   method = [m1,m2,0,50,0]; %single compare
   [eventValue, eventTime, eventID] = eventOnset(inputData,event,baseline,method);
   resultArray(eventID) = eventTime; 
   
   %IndexEMGMethod1 
   event = 'IndexEMGMethod1'; 
   baseline = [b1-b2+1,b2,3]; %2 times SD for threshold calculation
   method = [m1,m2,1,50,0]; %mean compare
   [eventValue, eventTime, eventID] = eventOnset(inputData,event,baseline,method);
   resultArray(eventID) = eventTime; 
   
   %IndexEMGMethod2
   event = 'IndexEMGMethod2'; 
   baseline = [b1-b2+1,b2,3]; %2 times SD for threshold calculation
   method = [m1,m2,2,0,0]; %mean window compare
   [eventValue, eventTime, eventID] = eventOnset(inputData,event,baseline,method);
   resultArray(eventID) = eventTime; 
                           
                          
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%% Output data to file %%%%%%%%%%%%
  
   %Write result data to file
   currentRow = outputDataToFile(outputFileName,resultArray);
   %Write input file name to file
   xlswrite(outputFileName,{inputFileName},1,strcat('A', currentRow));
   %write load condition of this trial
   xlswrite(outputFileName,{loadCondition},1,strcat('B', currentRow));
   
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Display if only one file as input%%%%%%%%%%%%  




