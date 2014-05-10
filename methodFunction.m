%function [eventValue, eventTime] = methodFunction(methodFunctionName)
% methodFunctionName
% singleCompare  
% meanCompare 
% meanWindowCompare 
switch methodFunctionName
    
    case 'singleCompare'
        for j = (method(1):(method(1) + method(2))) 
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
        
    case 'meanCompare'
        for j = (method(1):(method(1) + method(2)))
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
    case 'meanWindowCompare'
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
end