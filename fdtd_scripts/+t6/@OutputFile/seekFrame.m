function seekFrame(obj, frame)
% outputFile.seekFrame(1) will go to the first frame in the file.

frameBytes = obj.FrameSize * obj.BytesPerValue;
fseek(obj.FileHandle, (frame-1)*frameBytes, 'bof');

