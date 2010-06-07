function n = numFramesAvailable(obj)

if ~exist(obj.FileName)
    n = 0;
else
    %determine size of file in bytes (this is gross; Matlab's fault!!! lame API)
    d = dir;
    for nn = 1:length(d)
        if strcmp(d(nn).name, obj.FileName)
            n = floor( d(nn).bytes / obj.BytesPerValue / obj.FrameSize );
        end
    end
end
