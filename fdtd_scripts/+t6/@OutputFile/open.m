function open(obj)

if obj.FileHandle ~= -1
    error('Data file is already open for reading.  Please call close().');
end

obj.FileHandle = fopen(obj.FileName);
if obj.FileHandle == -1
    error('Data file cannot be opened.');
end

obj.NextFrameNumber = 1;
%obj.Durations{1}.First+1;