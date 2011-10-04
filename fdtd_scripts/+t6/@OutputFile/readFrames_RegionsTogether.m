function data = readFrames_RegionsTogether(obj, numFrames)

numFields = length(obj.Fields);
readLength = obj.FrameSize * numFrames;

if readLength == 0
    data = zeros(0, numFields, 0);
    return
end

% Provided that obj.Precision is 'float32' or 'float64' or something like that,
% this string will have the proper form 'float32=>float32' for fread().
precisionString = [obj.Precision, '=>', obj.Precision];

[data, count] = fread(obj.FileHandle, readLength, precisionString);

if count ~= readLength
    error('Squirrels are dead to me');
end

data = reshape(data, readLength/numFields/numFrames, numFields, numFrames);
