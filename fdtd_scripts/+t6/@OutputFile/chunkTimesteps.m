function [starts, ends, lengths] = chunkTimesteps(nFirst, nLast, frameSize, bytes)

bytesPerValue = 4;
frameBytes = sum(prod(frameSize,2))*bytesPerValue;

chunkTimesteps = max(1, floor(bytes / frameBytes));

starts = nFirst:chunkTimesteps:nLast;
ends = [starts(2:end)-1, nLast];
lengths = ends - starts + 1;



