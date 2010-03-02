function trace = getOutputPoint(filename, point)
%   Figure it out
%
%   version 4.5
%   July 29, 2008


[fid, dim] = openOutputFile(filename);

done = 0;
count = 0;

while (done == 0)
    [data, frameCount] = getOutputFrame(fid, dim);
    
    if (frameCount ~= 0)
        data = squeeze(data);
        trace(count+1) = data(point(1), point(2));
        count = count + 1;
    else
        done = 1;
    end
end

fclose(fid);