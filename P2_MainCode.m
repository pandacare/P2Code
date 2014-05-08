clc;
clear all;
close all;

F=1000; %Sampling frequency 1000 Hz
Fs=F/2;

%focus range (equals to 1 second plot on the figure)

%No.of points used to calculate the baseline left to the magnetRelease point  
n=100;

%lasting length for continious comparation with threshold event 
DURATION =30; %  window width used for EMG signal to determine onset
Continuation=50;% Continious poits beyound the threshold   ( EMG )
DURATION_F=100; %  window width used for grip signal to determine onset Force (Grip force )
GripForceDuration = 50; % window width used for grip force rate (Grip force rate)
delayTime =0;
%SD used to determine the threshold
num_std =3;

%channel option
channelForceXT = 1; %Thumb shear force
channelForceXI = 7; %Index shear force
channelForceZT = channelForceXT + 2; %Thumb grip
channelForceZI = channelForceXI + 2; %Index grip
channelEMGT = 13; % APB EMG
channelEMGI = 14; % FDI EMG

%cuttoff frequency for grip force
w1=20;
WnPass=w1/Fs;

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


%Open the data to be processed
%file name length should be 18, max file number as input is 10

[FileName,inputFilePath]=uigetfile('*.dat','Select the LabView Data','MultiSelect','on');
numberFile = length(FileName);


%Define parameter Header
parameterHeader =  {'ThumbPerturbation', 'ThumbGripForce','ThumbEMGMethod0','ThumbEMGMethod1','ThumbEMGMethod2','IndexPerturbation','IndexGripForce','IndexEMGMethod0','IndexEMGMethod1','IndexEMGMethod2','MaxThumbGrip','MaxThumbGripIndex','MaxThumbRate','MaxThumbRateIndex','MaxIndexGrip','MaxIndexGripIndex','MaxIndexRate','MaxIndexRateIndex'};
parameterNumber = length(parameterHeader);

%Nniitilize output parameter array 
resultArray = zeros(1, parameterNumber);


% Output file name
currentTime = fix(clock);
currentTimeNameElement = arrayfun(@num2str, currentTime, 'UniformOutput', false);
currentTimeName = strcat(currentTimeNameElement(1),currentTimeNameElement(2),currentTimeNameElement(3),currentTimeNameElement(4),currentTimeNameElement(5),currentTimeNameElement(6));
prompt = {'Enter file name:','Enter condition:','Enter FM/P2'};
dlg_title = '';
num_lines = 1;
def = {currentTimeName{1,1},'300','P2'};
outputFileNameElement = inputdlg(prompt,dlg_title,num_lines,def);
outputFileName = strcat(outputFileNameElement{1},'_',outputFileNameElement{2},'_',outputFileNameElement{3},'.xlsx');


%%%%%%%%%%%%%%%%%%%%%%%%% Data analysis %%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%% Single file  %%%%%%%%%%%%%%%%%%%%%%%%%%%%

 if (numberFile == 18 || numberFile == 1)
 
     f = fullfile(inputFilePath,FileName);
     A = importdata(f);
     
     %Get weight drop release point
     MagnetStart =round( A(1,1));
     
     %Set data analysis window length
     magnetRelease = MagnetStart;
     theTime1 = (magnetRelease -500)/1000;
     theTime2 = theTime1+1.5;
     
     %Transfer raw data into analysis array "dat"
     for i= 1:14
     dat(:,i)=A(startPoint+1:endPoint,i);
     end
     
     %Build filter for filtering raw data
     [b1,a1] = butter(4,WnPass,'low');
     
     % filter raw Force data and store in "datMod" 
     datMod(:,1)=filter(b1,a1,dat(:,1));
     datMod(:,2)=filter(b1,a1,dat(:,2));
     datMod(:,3)=filter(b1,a1,dat(:,3));
     datMod(:,4)=filter(b1,a1,dat(:,4));
     datMod(:,5)=filter(b1,a1,dat(:,5));
     datMod(:,6)=filter(b1,a1,dat(:,6));
     datMod(:,7)=filter(b1,a1,dat(:,7));
     datMod(:,8)=filter(b1,a1,dat(:,8));
     datMod(:,9)=filter(b1,a1,dat(:,9));
     datMod(:,10)=filter(b1,a1,dat(:,10));
     datMod(:,11)=filter(b1,a1,dat(:,11));
     datMod(:,12)=filter(b1,a1,dat(:,12));

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
figure % Shear force plot and the onset of the perturbation (Figure 1)

plot(dat((theTime1*1000+1):theTime2*1000,channelForceXT),'k');
title('Thumb Shear Force');
ylabel('Force (N)') ;% label for y axis
xlabel('Time (ms)'); % label for x axis
hold on

rec_xT = dat(:,channelForceXT); %shear force was analyzed with raw data
PxT = dat((baseLineStart+1):(n+baseLineStart),channelForceXT);
thresholdForceXT = mean(PxT) + std(PxT)*2;

% 5 indicats the 5ms window
% assumpation that onset only happens in a 500 ms range after the magnetic
% relase
for j = (magnetRelease +1):(theTime2*1000) 
    w = 0;
for i = j:(j+5-1)
    if rec_xT(i) > thresholdForceXT
      w = w + 1;  
    end   
end
    if w == 5
        finalIndexForceXT = j - theTime1*1000-1;
        break;
    end
    
end
absloteIndexPerturbationThumb = finalIndexForceXT+theTime1*1000;
plot(ones(21)*finalIndexForceXT,-10:1:10,'r');
hold off


       %%%%%%%%%%%%%%Thumb Grip %%%%%%%%%%%%%%%% 
figure % Figure 2

plot(datMod((theTime1*1000+1):theTime2*1000,channelForceZT),'k');
hold on
title('Thumb Grip force');
ylabel('Force (N)');% label for y axis
xlabel('Time (ms)');% label for x axis

plot(ones(21)*finalIndexForceXT,-10:1:10,'r');
hold on

rec_zT = datMod(:,channelForceZT);
PzT = datMod((baseLineStart+1):(n+baseLineStart),channelForceZT);
thresholdForceZT = mean(PzT)-std(PzT)*3;

for j = (magnetRelease + 1):(theTime2*1000)
    w = 0;
for i = j:(j+DURATION_F-1)
    if abs(rec_zT(i))> abs(thresholdForceZT)
      w = w + 1;  
    end   
end
    if w == DURATION_F
        finalIndexForceZT = j-theTime1*1000-1;
        break;
    end
end
[maxGripForceThumb, maxGripForceIndexThumb] = max(abs(rec_zT((magnetRelease + 1):(theTime2*1000))));
absloteIndexGripForceThumb = finalIndexForceZT + theTime1*1000;
plot(ones(21)*finalIndexForceZT,-10:1:10,'b');
hold on
plot(ones(21)*(maxGripForceIndexThumb+magnetRelease-theTime1*1000 -1 ),-10:1:10,'y');
hold on
y = rec_zT((finalIndexForceXT+theTime1*1000 + 1):(finalIndexForceXT+theTime1*1000 + 300));
[ymax,imax,ymin,imin] = extrema(y);
plot(imax+finalIndexForceXT,ymax,'r*',imin+finalIndexForceXT,ymin,'g*');
hold off
maxGripForceIndexThumb = maxGripForceIndexThumb + magnetRelease;
%plot(imax(ymaxindex1)+finalIndexForceXT,ymax1,'r*');



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

PzTRate = GripForceData((baseLineStart+1):(n+baseLineStart));
thresholdForceZTRate = mean(PzTRate)+ std(PzTRate)*3;
for j = (magnetRelease + 1):(theTime2*1000)
    w = 0;
for i = j:(j+DURATION_F)
    if abs(GripForceData(i))> abs(thresholdForceZT)
      w = w + 1;  
    end   
end
    if w == DURATION_F
        finalIndexForceZTRate = j-theTime1*1000 - 1;
        break;
    end
    
end
absloteIndexGripForceThumbRate = finalIndexForceZTRate+theTime1*1000;
[maxGripForceRateThumb, maxGripForceRateIndexThumb] = max(abs(GripForceData)) ;

figure %Figure 3 
plot(GripForceData((theTime1*1000+1):theTime2*1000),'k');
hold on
plot(ones(21)*finalIndexForceXT,-10:1:10,'r');
hold on
plot(ones(21)*(maxGripForceRateIndexThumb-theTime1*1000 - 1),-10:1:10,'y');
hold on
plot(ones(21)*finalIndexForceZTRate,-10:1:10,'b');
title('Thumb Grip force rate');
hold on
ylabel('Force (N)') ;% label for y axis
xlabel('Time (s)'); % label for x axis
hold off

         %%%%%%%%%%%%%%Thumb EMG%%%%%%%%%%%%%%%%
figure % Figure 4
title('EMG Thumb with filter');
ylabel('EMG (V)') ;% label for y axis
xlabel('Time (ms)'); % label for x axis
hold on
[b,a]=butter(4,50/Fs,'low');
[b5Hz,a5Hz]=butter(4,5/Fs,'low');
averageMeanThumb = mean(dat(meanStart:meanEnd,channelEMGT));
rec_EMGT = abs(dat(:,channelEMGT)- averageMeanThumb);
filter_EMGT=filtfilt(b,a,rec_EMGT);
filter5hz_EMGT=filtfilt(b5Hz,a5Hz,rec_EMGT);

PEMGT = filter_EMGT((baseLineStart+1):(n+baseLineStart));
thresholdEMGT = mean(PEMGT) + std(PEMGT)*num_std;

% Method 0 _Continous Duration points
for j = (absloteIndexPerturbationThumb +1):(theTime2*1000)
    w = 0;
for i = j:(j+Continuation-1)
    if filter_EMGT(i)> thresholdEMGT
      w = w + 1;  
    end   
end
    if w == Continuation
        finalIndexEMGT = j-theTime1*1000-1;
        break;
    end
    
end
absloteIndexEMGT = finalIndexEMGT+theTime1*1000;   


% Method average 1_the mean of continous Duration points comparing to
% threshold
for j = (absloteIndexPerturbationThumb +1):(theTime2*1000)
    ww = 0;
    www = 0;
for i = j:(j+DURATION-1)
    ww =ww+filter_EMGT(i);  
end
    www = ww/DURATION;
    if www > thresholdEMGT
        finalIndexEMGTAverage = j-theTime1*1000; 
        break;
    end
    
end
absloteIndexEMGTAverage = finalIndexEMGTAverage+theTime1*1000;  


% Method average 2_average 1+continous Continuation seconds
for j = (absloteIndexPerturbationThumb +1+delayTime):(theTime2*1000)
    ww2 = 0;
    kk2 =0;
for i = j:(j+Continuation-1)
    ww2 = 0;
    for ii = i:(i+DURATION-1)
    ww2 =ww2+filter_EMGT(ii);  
    end
    if (ww2/DURATION) > thresholdEMGT
        kk2 = kk2 +1;
    end
end    
    if kk2 ==Continuation
        finalIndexEMGTAverage2 = j-theTime1*1000-1;
        break;
    end    
end

absloteIndexEMGTAverage2 = finalIndexEMGTAverage2+theTime1*1000;  


plot(filter_EMGT((theTime1*1000+1):theTime2*1000), 'k');
hold on
% plot(filter5hz_EMGT((theTime1*1000+1):theTime2*1000), 'k');
% hold on
plot(ones(21)*finalIndexEMGTAverage,0:0.005:0.1,'g');% Average 1 with green
hold on
plot(ones(21)*finalIndexEMGTAverage2,0:0.005:0.1,'y');%Average 2 with yellow
hold on
plot(ones(21)*finalIndexForceXT,0:0.005:0.1,'r');% Onset of the Perturbation
hold on
plot(ones(21)*finalIndexEMGT,0:0.005:0.1,'b'); % Method 0 with blue
hold off

figure % Figure 5


plot(rec_EMGT((theTime1*1000+1):theTime2*1000),'k');
hold on
title('EMG Thumb withoutfilter');
ylabel('EMG (V)') ;% label for y axis
xlabel('Time (ms)'); % label for x axis

plot(ones(21)*finalIndexForceXT,0:0.005:0.1,'r');
hold on
plot(ones(21)*finalIndexEMGT,0:0.005:0.1,'b');
hold on
plot(ones(21)*finalIndexEMGTAverage,0:0.005:0.1,'g');
hold on
plot(ones(21)*finalIndexEMGTAverage2,0:0.005:0.1,'y');
hold off


%---------------------------Index----------------------------%
            %%%%%%%%%%%%%%Index Shear%%%%%%%%%%%%%%%%        
figure % Figure 6
plot(dat((theTime1*1000+1):theTime2*1000,channelForceXI),'k');
title('Index Shear Force_X');
ylabel('Force (N)') ;% label for y axis
xlabel('Time (ms)'); % label for x axis
hold on

rec_xI = dat(:,channelForceXI);
PxI = dat((baseLineStart+1):(n+baseLineStart),channelForceXI);
thresholdForceXI = mean(PxI) + std(PxI)*2;
for j = (magnetRelease +1):(theTime2*1000)
    w = 0;
for i = j:(j+5-1)
    if rec_xI(i)> thresholdForceXI
      w = w + 1;  
    end   
end
    if w == 5
        finalIndexForceXI = j-theTime1*1000 - 1;
        break;
    end
    
end
absloteIndexPerturbationIndex = finalIndexForceXI+theTime1*1000;
plot(ones(21)*finalIndexForceXI,-10:1:10,'r');
hold off





        %%%%%%%%%%%%%%Index Grisp%%%%%%%%%%%%%%%%  
 %%%%%%%%%%%%%%Index Grisp%%%%%%%%%%%%%%%%  
figure % Figure 7
title('Index Grip force');
ylabel('Force (N)') ;% label for y axis
xlabel('Time (ms)'); % label for x axis
hold on
plot(datMod((theTime1*1000+1):theTime2*1000,channelForceZI),'k');
hold on
plot(ones(21)*finalIndexForceXI,-10:1:10,'r');
hold on
rec_zI = datMod(:,channelForceZI);
PzI = datMod((baseLineStart+1):(n+baseLineStart),channelForceZI);
thresholdForceZI = mean(PzI) -std(PzI)*3;
for j = (magnetRelease + 1):(theTime2*1000)
    w = 0;
for i = j:(j+DURATION_F-1)
    if abs(rec_zI(i))> abs(thresholdForceZI)
      w = w + 1;  
    end   
end
    if w == DURATION_F
        finalIndexForceZI = j-theTime1*1000 - 1;
        break;
    end
    
end
absloteIndexGripForceIndex = finalIndexForceZI+theTime1*1000;



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[maxGripForceIndex, maxGripForceIndexIndex] = max(abs(rec_zI((magnetRelease + 1):(theTime2*1000)))) ;

plot(ones(21)*finalIndexForceZI,-10:1:10,'b');
hold on
plot(ones(21)*(maxGripForceIndexIndex+ magnetRelease -theTime1*1000 - 1),-10:1:10,'y');
hold off
maxGripForceIndexIndex = maxGripForceIndexIndex + magnetRelease ; 
        %%%%%%%%%%%%%%Index grip rate %%%%%%%%%%%%%%%% 
h = waitbar(0, 'waiting ...');
for j = theTime1*1000:theTime2*1000
    a1 = 0;
    a2 = 0;
    a3 = 0;
    a4 = 0;
    for i = j:(j+GripForceDuration-1)
        a1 = rec_zI(i)*((i-j+1)/F)+a1;
        a2 = rec_zI(i)+a2;
        a3 = ((i-j+1)/F)+a3;
        a4 = ((i-j+1)/F)^2+a4;
    end
    GripForceDataIndexFinger(j) = (GripForceDuration*a1-a2*a3)/(GripForceDuration*a4-a3*a3);
    waitbar((j-theTime1*1000)/(theTime2*1000-theTime1*1000),h);
end
delete(h);

PzIRate = GripForceDataIndexFinger((baseLineStart+1):(n+baseLineStart));
thresholdForceZIRate = mean(PzIRate)+ std(PzIRate)*3;
for j = (magnetRelease + 1):(theTime2*1000)
    w = 0;
for i = j:(j+DURATION_F)
    if abs(GripForceDataIndexFinger(i))> abs(thresholdForceZT)
      w = w + 1;  
    end   
end
    if w == DURATION_F
        finalIndexForceZIRate = j-theTime1*1000 - 1;
        break;
    end
    
end
absloteIndexGripForceIndexRate = finalIndexForceZIRate+theTime1*1000;
[maxGripForceRateIndex, maxGripForceRateIndexIndex] = max(abs(GripForceDataIndexFinger)) ;

figure % Figure 8
plot(GripForceDataIndexFinger((theTime1*1000+1):theTime2*1000),'k');
hold on
plot(ones(21)*finalIndexForceXI,-10:1:10,'r');
hold on
plot(ones(21)*(maxGripForceRateIndexIndex-theTime1*1000 - 1),-10:1:10,'y');
hold on
plot(ones(21)*finalIndexForceZIRate,-10:1:10,'b');
title('Index Grip force rate');
ylabel('Force (N)') ;% label for y axis
xlabel('Time (s)'); % label for x axis
hold off




        %%%%%%%%%%%%%%Index EMG%%%%%%%%%%%%%%%%  
figure % Figure 9


[b,a]=butter(4,50/Fs,'low');

averageMeanIndex = mean(dat(meanStart:meanEnd,channelEMGI));
rec_EMGI = abs(dat(:,channelEMGI)- averageMeanIndex);
filter_EMGI=filtfilt(b,a,rec_EMGI);

% Method0_Continous Duration points
PEMGI = filter_EMGI((baseLineStart+1):(n+baseLineStart));
thresholdEMGI = mean(PEMGI) + std(PEMGI)*num_std;
for j = (absloteIndexPerturbationIndex +1):(theTime2*1000)
    w = 0;
for i = j:(j+Continuation-1)
    if filter_EMGI(i)> thresholdEMGI
      w = w + 1;  
    end   
end
    if w == Continuation
        finalIndexEMGI = j-theTime1*1000-1;
        break;
    end
    
end
absloteIndexEMGI = finalIndexEMGI+theTime1*1000;   

%Method average 1
for j = (absloteIndexPerturbationIndex +1+delayTime):(theTime2*1000)
    ww = 0;
    www = 0;
for i = j:(j+DURATION-1)
   
     ww =ww+filter_EMGI(i);  
     
end
    www = ww/DURATION;
    if www > thresholdEMGI
        finalIndexEMGIAverage = j-theTime1*1000-1;
        break;
    end
    
end
absloteIndexEMGIAverage = finalIndexEMGIAverage+theTime1*1000;  



%Method average 2
for j = (absloteIndexPerturbationIndex +1+delayTime):(theTime2*1000)
    ww2 = 0;
    kk2 =0;
for i = j:(j+Continuation-1)
    ww2 = 0;
    for ii = i:(i+DURATION-1)
     ww2 =ww2+filter_EMGI(ii);  
    end
    if (ww2/DURATION) > thresholdEMGI
        kk2 = kk2 +1;
    end
end
    
    if kk2 ==Continuation
        finalIndexEMGIAverage2 = j-theTime1*1000-1;
        break;
    end
    
end
absloteIndexEMGIAverage2 = finalIndexEMGIAverage2+theTime1*1000;  



plot(filter_EMGI((theTime1*1000+1):theTime2*1000), 'k');
title('EMG Index with filter');
ylabel('EMG (V)') ;% label for y axis
xlabel('Time (ms)'); % label for x axis
hold on
plot(ones(21)*finalIndexForceXI,0:0.005:0.1,'r');
hold on
plot(ones(21)*finalIndexEMGI,0:0.005:0.1,'b');
hold on
plot(ones(21)*finalIndexEMGIAverage,0:0.005:0.1,'g');
hold on
plot(ones(21)*finalIndexEMGIAverage2,0:0.005:0.1,'y');
hold off

figure %Fig 10
plot(rec_EMGI((theTime1*1000+1):theTime2*1000),'k');
title('EMG Index withoutfilter');
ylabel('EMG (V)') ;% label for y axis
xlabel('Time (ms)'); % label for x axis
hold on
plot(ones(21)*finalIndexForceXI,0:0.005:0.1,'r');
hold on
plot(ones(21)*finalIndexEMGI,0:0.005:0.1,'b');
hold on
plot(ones(21)*finalIndexEMGIAverage,0:0.005:0.1,'g');
hold on
plot(ones(21)*finalIndexEMGIAverage2,0:0.005:0.1,'y');
hold off


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



%multiple file input
     
 else
     fileNameString = {0;0;0;0;0;0;0;0;0;0;0;0};
     resultArray = zeros(numberFile, parameterNumber);
     myString = cell(1);
     for k = 1:numberFile
         f = fullfile(inputFilePath,FileName{k});
         A = importdata(f);
         %get weight drop release point
         MagnetStart =round( A(1,1));
         
         magnetRelease = MagnetStart;
         theTime1 = (magnetRelease -500)/1000;
         theTime2 = theTime1+1.5;
         
         
         for i= 1:14
            dat(:,i)=A(startPoint+1:endPoint,i);
         end
        [b1,a1] = butter(4,WnPass,'low');
        datMod(:,1)=filter(b1,a1,dat(:,1));
        datMod(:,2)=filter(b1,a1,dat(:,2));
        datMod(:,3)=filter(b1,a1,dat(:,3));
        datMod(:,4)=filter(b1,a1,dat(:,4));
        datMod(:,5)=filter(b1,a1,dat(:,5));
        datMod(:,6)=filter(b1,a1,dat(:,6));
        datMod(:,7)=filter(b1,a1,dat(:,7));
        datMod(:,8)=filter(b1,a1,dat(:,8));
        datMod(:,9)=filter(b1,a1,dat(:,9));
        datMod(:,10)=filter(b1,a1,dat(:,10));
        datMod(:,11)=filter(b1,a1,dat(:,11));
        datMod(:,12)=filter(b1,a1,dat(:,12));


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

rec_xT = dat(:,channelForceXT);
PxT = dat((baseLineStart+1):(n+baseLineStart),channelForceXT);
thresholdForceXT = mean(PxT) + std(PxT)*2;

% 5 indicats the 5ms window
% assumpation that onset only happens in a 500 ms range after the magnetic
% relase
for j = (magnetRelease +1):(theTime2*1000) 
    w = 0;
for i = j:(j+5-1)
    if rec_xT(i)> thresholdForceXT;
      w = w + 1;  
    end   
end
    if w == 5
        finalIndexForceXT = j-theTime1*1000 - 1;
        break;
    end
    
end
absloteIndexPerturbationThumb = finalIndexForceXT+theTime1*1000;


        %%%%%%%%%%%%%%Thumb Grip %%%%%%%%%%%%%%%% 
        
rec_zT = datMod(:,channelForceZT);
PzT = datMod((baseLineStart+1):(n+baseLineStart),channelForceZT);
thresholdForceZT = mean(PzT)-std(PzT)*3;
for j = (magnetRelease + 1):(theTime2*1000)
    w = 0;
for i = j:(j+DURATION_F-1)
    if abs(rec_zT(i))> abs(thresholdForceZT)
      w = w + 1;  
    end   
end
    if w == DURATION_F
        finalIndexForceZT = j-theTime1*1000 - 1;
        break;
    end
end
[maxGripForceThumb, maxGripForceIndexThumb] = max(abs(rec_zT((magnetRelease + 1):(theTime2*1000))));
maxGripForceIndexThumb = maxGripForceIndexThumb+magnetRelease;
absloteIndexGripForceThumb = finalIndexForceZT + theTime1*1000;

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

PzTRate = GripForceData((baseLineStart+1):(n+baseLineStart));
thresholdForceZTRate = mean(PzTRate)+ std(PzTRate)*3;
for j = (magnetRelease + 1):(theTime2*1000)
    w = 0;
for i = j:(j+DURATION_F)
    if abs(GripForceData(i))> abs(thresholdForceZT)
      w = w + 1;  
    end   
end
    if w == DURATION_F
        finalIndexForceZTRate = j-theTime1*1000 - 1;
        break;
    end
    
end
absloteIndexGripForceThumbRate = finalIndexForceZTRate+theTime1*1000;
[maxGripForceRateThumb, maxGripForceRateIndexThumb] = max(abs(GripForceData)) ;



         %%%%%%%%%%%%%%Thumb EMG%%%%%%%%%%%%%%%%
[b,a]=butter(4,50/Fs,'low');
averageMeanThumb = mean(dat(meanStart:meanEnd,channelEMGT));
rec_EMGT = abs(dat(:,channelEMGT)- averageMeanThumb);
filter_EMGT=filtfilt(b,a,rec_EMGT);

PEMGT = filter_EMGT((baseLineStart+1):(n+baseLineStart));
thresholdEMGT = mean(PEMGT) + std(PEMGT)*num_std;

% Method 0 _Continous Duration points
for j = (absloteIndexPerturbationThumb +1):(theTime2*1000)
    w = 0;
for i = j:(j+Continuation-1)
    if filter_EMGT(i)> thresholdEMGT
      w = w + 1;  
    end   
end
    if w == Continuation
        finalIndexEMGT = j-theTime1*1000-1;
        break;
    end
    
end
absloteIndexEMGT = finalIndexEMGT+theTime1*1000;   


% Method average 1_the mean of continous Duration points comparing to
% threshold
for j = (absloteIndexPerturbationThumb +1):(theTime2*1000)
    ww = 0;
    www = 0;
for i = j:(j+DURATION-1)
   
     ww =ww+filter_EMGT(i);  
     
end
    www = ww/DURATION;
    if www > thresholdEMGT
        finalIndexEMGTAverage = j-theTime1*1000-1; 
        break;
    end
    
end
absloteIndexEMGTAverage = finalIndexEMGTAverage+theTime1*1000;  


% Method average 2_average 1+continous Continuation seconds
for j = (absloteIndexPerturbationThumb +1+delayTime):(theTime2*1000)
    ww2 = 0;
    kk2 =0;
for i = j:(j+Continuation-1)
    ww2 = 0;
    for ii = i:(i+DURATION-1)
    ww2 =ww2+filter_EMGT(ii);  
    end
    if (ww2/DURATION) > thresholdEMGT
        kk2 = kk2 +1;
    end
end    
    if kk2 ==Continuation
        finalIndexEMGTAverage2 = j-theTime1*1000-1;
        break;
    end    
end

absloteIndexEMGTAverage2 = finalIndexEMGTAverage2+theTime1*1000;  


%---------------------------Index----------------------------%
            %%%%%%%%%%%%%%Index Shear%%%%%%%%%%%%%%%%        

rec_xI = dat(:,channelForceXI);
PxI = dat((baseLineStart+1):(n+baseLineStart),channelForceXI);
thresholdForceXI = mean(PxI) + std(PxI)*2;
for j = (magnetRelease +1):(theTime2*1000)
    w = 0;
for i = j:(j+5-1)
    if rec_xI(i)> thresholdForceXI
      w = w + 1;  
    end   
end
    if w == 5
        finalIndexForceXI = j-theTime1*1000 - 1;
        break;
    end
    
end
absloteIndexPerturbationIndex = finalIndexForceXI+theTime1*1000;







        %%%%%%%%%%%%%%Index Grisp%%%%%%%%%%%%%%%%  

rec_zI = datMod(:,channelForceZI);
PzI = datMod((baseLineStart+1):(n+baseLineStart),channelForceZI);
thresholdForceZI = mean(PzI) -std(PzI)*3;
for j = (magnetRelease + 1):(theTime2*1000)
    w = 0;
for i = j:(j+DURATION_F-1)
    if abs(rec_zI(i))> abs(thresholdForceZI)
      w = w + 1;  
    end   
end
    if w == DURATION_F
        finalIndexForceZI = j-theTime1*1000 - 1;
        break;
    end
    
end
absloteIndexGripForceIndex = finalIndexForceZI+theTime1*1000;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[maxGripForceIndex, maxGripForceIndexIndex] = max(abs(rec_zI((magnetRelease + 1):(theTime2*1000)))) ;
maxGripForceIndexIndex = maxGripForceIndexIndex + magnetRelease;

        %%%%%%%%%%%%%%Index grip rate %%%%%%%%%%%%%%%% 
h = waitbar(0, 'waiting ...');
for j = theTime1*1000:theTime2*1000
    a1 = 0;
    a2 = 0;
    a3 = 0;
    a4 = 0;
    for i = j:(j+GripForceDuration-1)
        a1 = rec_zI(i)*((i-j+1)/F)+a1;
        a2 = rec_zI(i)+a2;
        a3 = ((i-j+1)/F)+a3;
        a4 = ((i-j+1)/F)^2+a4;
    end
    GripForceDataIndexFinger(j) = (GripForceDuration*a1-a2*a3)/(GripForceDuration*a4-a3*a3);
    waitbar((j-theTime1*1000)/(theTime2*1000-theTime1*1000),h);
end
delete(h);

PzIRate = GripForceDataIndexFinger((baseLineStart+1):(n+baseLineStart));
thresholdForceZIRate = mean(PzIRate)+ std(PzIRate)*3;
for j = (magnetRelease + 1):(theTime2*1000)
    w = 0;
for i = j:(j+DURATION_F)
    if abs(GripForceDataIndexFinger(i))> abs(thresholdForceZT)
      w = w + 1;  
    end   
end
    if w == DURATION_F
        finalIndexForceZIRate = j-theTime1*1000 - 1;
        break;
    end
    
end
absloteIndexGripForceIndexRate = finalIndexForceZIRate+theTime1*1000;
[maxGripForceRateIndex, maxGripForceRateIndexIndex] = max(abs( GripForceDataIndexFinger)) ;







        %%%%%%%%%%%%%%Index EMG%%%%%%%%%%%%%%%%  

[b,a]=butter(4,50/Fs,'low');
averageMeanIndex = mean(dat(meanStart:meanEnd,channelEMGI));
rec_EMGI = abs(dat(:,channelEMGI)- averageMeanIndex);
filter_EMGI=filtfilt(b,a,rec_EMGI);

% Method0_Continous Duration points
PEMGI = filter_EMGI((baseLineStart+1):(n+baseLineStart));
thresholdEMGI = mean(PEMGI) + std(PEMGI)*num_std;
for j = (absloteIndexPerturbationIndex +1):(theTime2*1000-Continuation)
    w = 0;
for i = j:(j+Continuation-1)
    if filter_EMGI(i)> thresholdEMGI
      w = w + 1;  
    end   
end
    if w == Continuation
        finalIndexEMGI = j-theTime1*1000-1;
        break;
    end
    
end
absloteIndexEMGI = finalIndexEMGI+theTime1*1000;   


%Method average 1
for j = (absloteIndexPerturbationIndex +1+delayTime):(theTime2*1000)
    ww = 0;
    www = 0;
for i = j:(j+DURATION-1)
   
     ww =ww+filter_EMGI(i);  
     
end
    www = ww/DURATION;
    if www > thresholdEMGI
        finalIndexEMGIAverage = j-theTime1*1000-1;
        break;
    end
    
end
absloteIndexEMGIAverage = finalIndexEMGIAverage+theTime1*1000;  


%Method average 2
for j = (absloteIndexPerturbationIndex +delayTime+1):(theTime2*1000)
    ww2 = 0;
    kk2 =0;
for i = j:(j+Continuation-1)
    for ii = i:(i+DURATION-1)
     ww2 =ww2+filter_EMGI(ii);  
    end
    if (ww2/DURATION) > thresholdEMGI
        kk2 = kk2 +1;
    end
end
    
    if kk2 ==Continuation
        finalIndexEMGIAverage2 = j-theTime1*1000-1;
        break;
    end
    
end
absloteIndexEMGIAverage2 = finalIndexEMGIAverage2+theTime1*1000; 

fileNameString{k,1} = { FileName{1,k}(1:14)};
resultArray(k,1)    = absloteIndexPerturbationThumb;
resultArray(k,2)    = absloteIndexGripForceThumb;
resultArray(k,3)    = absloteIndexEMGT;
resultArray(k,4)    = absloteIndexEMGTAverage ;
resultArray(k,5)    = absloteIndexEMGTAverage2;
resultArray(k,6)    = absloteIndexPerturbationIndex;
resultArray(k,7)    = absloteIndexGripForceIndex ;
resultArray(k,8)    = absloteIndexEMGI;
resultArray(k,9)    = absloteIndexEMGIAverage;
resultArray(k,10)   = absloteIndexEMGIAverage2;
resultArray(k,11)   = maxGripForceThumb;
resultArray(k,12)   = maxGripForceIndexThumb;
resultArray(k,13)   = maxGripForceRateThumb;
resultArray(k,14)   = maxGripForceRateIndexThumb ;
resultArray(k,15)   = maxGripForceIndex ;
resultArray(k,16)   = maxGripForceIndexIndex;
resultArray(k,17)   = maxGripForceRateIndex;
resultArray(k,18)   = maxGripForceRateIndexIndex;

     end 
    xlswrite(outputFileName,fileNameString{1},1,'A2');
    xlswrite(outputFileName,fileNameString{2},1,'A3');
    xlswrite(outputFileName,fileNameString{3},1,'A4');
    xlswrite(outputFileName,fileNameString{4},1,'A5');
    xlswrite(outputFileName,fileNameString{5},1,'A6');
    xlswrite(outputFileName,fileNameString{6},1,'A7');
    xlswrite(outputFileName,fileNameString{7},1,'A8');
    xlswrite(outputFileName,fileNameString{8},1,'A9');
    xlswrite(outputFileName,fileNameString{9},1,'A10');
    xlswrite(outputFileName,fileNameString{10},1,'A11');
    xlswrite(outputFileName,fileNameString{11},1,'A12');
    
    xlswrite(outputFileName,parameterHeader,1,'C1');
    xlswrite(outputFileName,resultArray,1,'C2');
   xlswrite(outputFileName,{'condition'},1,'B1');
   xlswrite(outputFileName,{outputFileNameElement{2}},1,'B2');
   xlswrite(outputFileName,{'FileName'},1,'A1');
    
 end
% export result to excel

%get file name

     