function data = readFrames_RegionsTogether(obj, numFrames)

numFields = length(obj.Fields);
readLength = obj.FrameSize * numFrames;

% Provided that obj.Precision is 'float32' or 'float64' or something like that,
% this string will have the proper form 'float32=>float32' for fread().
precisionString = [obj.Precision, '=>', obj.Precision];

data = fread(obj.FileHandle, readLength, precisionString);
data = reshape(data, readLength/numFields/numFrames, numFields, numFrames);
