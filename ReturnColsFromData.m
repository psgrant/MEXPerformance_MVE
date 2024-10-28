function [ScaledCols,CMin,CMax] = ReturnColsFromData(Data,map,CMin,CMax)
% Returns the associated colour for a given data point and the data mins
% and maxes


% Check if the number of input arguments is less than 1
if nargin <= 2
    CMin = min(Data);
    CMax = max(Data);
end

% DataTemp = rescale([Data; CMin; CMax]);
if CMax == CMin
    y = 0.5+0*Data;
    CMax = mean(Data)*0.1 + CMax;
    CMin = CMin - mean(Data)*0.1;
    if CMax == 0
        CMin = -1e-8;
        CMax = 1e8;
    end
else
    y = (Data  - CMin )/(CMax - CMin);
end

y(y < 0) = 0;
y(y > 1) = 1;
Data = y;
ScaledCols = map(round((Data)*(length(map)-1))+1,:);
