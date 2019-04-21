function [ areaSize ] = calRefSize( obj, rows, cols )
%CALREFSIZE Compute the refrence size of small target based 
%   on image size.(0.15%)
%   Input: Length and width of an image.
%   Output: The refrence size of infrared small target.
totalArea = rows * cols;
areaSize = totalArea * obj.tarRate;

end

