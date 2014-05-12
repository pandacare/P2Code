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
startPoint =0;                    %signal start reference ms
endPoint =44000;                  %signal stop reference ms
lengthTime=endPoint-startPoint;
dat=zeros(lengthTime,14);         %Matrix for holding raw data
datMod=zeros(lengthTime,14);      %Matrix for holding filtered data
lengthData = lengthTime;
GripForceData = zeros(lengthData,1);
GripForceDataIndexFinger = zeros(lengthData,1);
meanStart = 11;
meanEnd   = 110; %for offset cancellation


%Define parameter Header
parameterArray = ['InputFileName',...
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
                  'IndexMaxGripForceRateTimePoint'];
              
%Define parameter is nonData amount "inputFileName" and "loadCondition" in
%this case
nonDataparameters = 2; 

%Calculate parameterArray length
resultArray = zeros(1, (length(parameterArray)-nonDataparameters));


% Define Output file name



%%%%%%%%%%%%%%%%%%%%%%%%% Data analysis %%%%%%%%%%%%%%%%%%%%%%%%%%%%

%call inputfile
%read data to A
%

     %variable definitions
     baseLineStart = magnetRelease-n ;
     
     perturbationOnsetThumb = 1;
     GripForceResponseOnsetThumb = 1;
     EMGResponseThumbMethod0 = 1;
     EMGResponseThumbMethod1 = 1;
     EMGResponseThumbMethod2 = 1;
     
     perturbationOnsetIndex = 1;
     GripForceResponseOnsetIndex = 1;
     EMGResponseIndexMethod0 = 1;
     EMGResponseIndexMethod1 = 1;
     EMGResponseIndexMethod2 = 1;
     gripForceRateMaxThumb   = 1;
     gripForceRateMaxIndex   = 1;
     
     
     finalIndexForceXT       = perturbationOnsetThumb;
     finalIndexForceZT       = GripForceResponseOnsetThumb;
     finalIndexEMGT          = EMGResponseThumbMethod0;
     finalIndexEMGTAverage2  = EMGResponseThumbMethod1;
     finalIndexEMGTAverage   = EMGResponseThumbMethod2;
     finalIndexForceXI       = perturbationOnsetIndex;
     finalIndexForceZI       = GripForceResponseOnsetIndex;
     finalIndexEMGI          = EMGResponseIndexMethod0;
     finalIndexEMGIAverage2  = EMGResponseIndexMethod1;
     finalIndexEMGIAverage   = EMGResponseIndexMethod2;
     finalIndexForceZTRate   = gripForceRateMaxThumb;
     finalIndexForceZIRate   = gripForceRateMaxIndex;
     
     
     

%---------------------------Thumb----------------------------%
             %%%%%%%%%%%%%%Thumb Shear%%%%%%%%%%%%%%%%




        %%%%%%%%%%%%%%Thumb Grip rate %%%%%%%%%%%%%%%% 
h = waitbar(0, 'waiting ...');
for j = theTime1*1000:theTime2*1000
    a1 = 0;
    a2 = 0;
    a3 = 0;
    a4 = 0;
    for i = j:(j+GripForceDuration-1)
        a1 = rec_zT(i)*((i-j+1)/F)+a1;
        a2 = rec_zT(i)+a2;
        a3 = ((i-j+1)/F)+a3;
        a4 = ((i-j+1)/F)^2+a4;
    end
    GripForceData(j) = (GripForceDuration*a1-a2*a3)/(GripForceDuration*a4-a3*a3);
    waitbar((j-theTime1*1000)/(theTime2*1000-theTime1*1000),h);
end
delete(h);

         %%%%%%%%%%%%%%Thumb EMG%%%%%%%%%%%%%%%%

% Method 0 _Continous Duration points

% Method average 1_the mean of continous Duration points comparing to
% threshold


% Method average 2_average 1+continous Continuation seconds

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%Repeat for Index %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


resultArray(1)  =  absloteIndexPerturbationThumb;
resultArray(2)  =  absloteIndexGripForceThumb;
resultArray(3)  =  absloteIndexEMGT;
resultArray(4)  =  absloteIndexEMGTAverage ;
resultArray(5)  =  absloteIndexEMGTAverage2;
resultArray(6)  =  absloteIndexPerturbationIndex ;
resultArray(7)  =  absloteIndexGripForceIndex ;
resultArray(8)  =  absloteIndexEMGI;
resultArray(9)  =  absloteIndexEMGIAverage;
resultArray(10) =  absloteIndexEMGIAverage2;
resultArray(11) =  maxGripForceThumb;
resultArray(12) =  maxGripForceIndexThumb;
resultArray(13) =  maxGripForceRateThumb;
resultArray(14) =  maxGripForceRateIndexThumb ;
resultArray(15) =  maxGripForceIndex ;
resultArray(16) =  maxGripForceIndexIndex;
resultArray(17) =  maxGripForceRateIndex;
resultArray (18)=  maxGripForceRateIndexIndex;

%get file name
xlswrite(outputFileName,parameterHeader,1,'C1');
xlswrite(outputFileName,resultArray,1,'C2');
xlswrite(outputFileName,{'FileName'},1,'A1');
myString = {FileName};
newString ={ myString{1}(1:14)};
xlswrite(outputFileName,newString,1,'A2');
xlswrite(outputFileName,{'condition'},1,'B1');
xlswrite(outputFileName,{outputFileNameElement{2}},1,'B2');



