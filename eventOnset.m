
function [eventValue, eventTime, eventID] = eventOnset(inputData,event,baseline,method)
% event:
% perturbation:  
% gripForce: 
% gripForceRate: 
% EMG: 
% baseline[3]:  baseLineStart, baseLineLength, nTimesSD
%%%%%%%%
% method[5]: searchRangeStart, searchRangeLength, methodFunctionName,duration1, duration2
% methodFunctionName
% singleCompare  0 
% meanCompare   1
% meanWindowCompare 2
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

global parameterArray;

%initial eventTime and eventValue
eventValue = -1;
eventTime  = -1;
eventID = -1;

%BaseLineData for baseLine window
baseLineData = inputData(baseline(1):(baseline(1) + baseline(2)-1));
        
%Calculation threshlod
threshhold = mean( baseLineData) + std( baseLineData)*baseline(3);

%Setup eventID
for i = 1 : length(parameterArray)
    if strcmp(parameterArray{i}, event)
        eventID = i;
        break;
    end
     %error('Wrong parameter number');
end

switch method(3)
    %'singleCompare'
    case 0
        for j = ((method(1)+1):(method(1) + method(2))) 
            w = 0;
            for i = j:(j+method(4)-1)
                if inputData(i) > threshhold
                w = w + 1;  
                end   
            end
            if w == method(4)
                eventTime  = j;
                eventValue = inputData(j);
                break;
            end
        end
    %'meanCompare'   
    case 1
        for j = ((method(1)+1):(method(1) + method(2)))
            w = 0;
            for i = j:(j+method(4)-1)
                w =w+inputData(i);  
            end
            www = w/method(4);
            if www > threshhold
                eventTime  = j; 
                eventValue = inputData(j);
                break;
            end
        end
    %'meanWindowCompare'    
    case 2
        for j = ((method(1)+1):(method(1) + method(2)))
            k =0;
            for i = j:(j+method(5)-1)
                w = 0;
                for ii = i:(i+method(4)-1)
                    w =w+inputData(ii);  
                end
                if (w/method(4)) > threshold
                    k = k +1;
                end
            end    
            if k == method(5)
               eventTime  = j;
               eventValue = inputData(j);
               break;
            end    
        end
    %'max'
    case 3
    [eventValue, eventTime] = max(inputData((method(1)+1):(method(1) + method(2)))) ;    
    otherwise 
        disp('wrong method setup');
end

