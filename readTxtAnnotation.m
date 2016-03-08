function [ time, pitch ] = readTxtAnnotation( path )
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here
sizeTimePitchProb = [3 Inf];
s = strcat(path,'.txt')
fid = fopen(s);
timePitchProb = fscanf(fid,'%f %f %f',sizeTimePitchProb);
fclose(fid);
time = timePitchProb(1,:);
pitch = timePitchProb(2,:);
end

