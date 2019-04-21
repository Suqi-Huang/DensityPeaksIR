function [ conf ] = confidenceCal( ~, inputVector )
%CONFIDENCECAL Summary of this function goes here
%   Input: The GVR values of density peaks (candidate targets).
%   Output: The confidence of each density peak. 
%   (if a density peak's confidence is larger than thdQuatile, 
% then it can be recognized as a real target. )

num = length(inputVector);
[ inputVec, ~ ] = sort(inputVector, 'ascend');
q1 = inputVec( round(num/4) );
q3 = inputVec( round(3*num/4) );
iqr = q3 - q1;

% In special cases (extremely flat region),  iqr equals to 0;
% we give a simple solution here to avoid 0 in denominator.
if (iqr == 0) 
    iqr = 2;
end

k = (inputVector - q3) / iqr;
conf = k;
end

