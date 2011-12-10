function open(obj)

if obj.FileHandle ~= -1
    error('Data file is already open for reading.  Please call close().');
end

obj.FileHandle = fopen(obj.FileName, 'r');
if obj.FileHandle == -1
    error('Data file cannot be opened.');
end

