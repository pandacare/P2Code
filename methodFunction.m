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
        for j = (method(1):(method(1) + method(2)))
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
end