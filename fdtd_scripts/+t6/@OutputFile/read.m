function data = read(obj, varargin)
% This is my function!

X.Regions = 'Separate';
X = parseargs(X, varargin{:});

% What does this do?
% It reads all the fields!  All of them!  A quick but slow approach
% would use readFrame.  I'll try that first.

obj.open
try

numFrames = obj.numFramesAvailable;

MAXREADBYTES = 1e8;   % Read up to 100 MB at a time, to reduce memory overhead.
bytesPerFrame = obj.FrameSize * obj.BytesPerValue;
framesPerChunk = min(numFrames, ceil(MAXREADBYTES / bytesPerFrame) );
numChunks = floor(numFrames / framesPerChunk);

% Reserve some space
% size: [x y z f t r]
if length(obj.Regions) > 1 && strcmp(X.Regions, 'Separate')
    data = cell(size(obj.Regions));
    for rr = 1:length(obj.Regions)
        data{rr} = zeros([obj.Regions{rr}.Size, length(obj.Fields), ...
            numFrames]);
    end
elseif length(obj.Regions) == 1 && strcmp(X.Regions, 'Separate')
    data = zeros([obj.Regions{1}.Size, length(obj.Fields), numFrames]);
elseif ~strcmp(X.Regions, 'Separate')
    data = zeros(obj.FrameSize/length(obj.Fields), length(obj.Fields), numFrames);
else
    error('No Regions available.  Read GridReports with readFrames(1).');
end

% Copy in the data by chunks as far as possible
for chunk = 1:numChunks
    frameNum = (chunk-1)*framesPerChunk + 1;
    frameRange = [frameNum frameNum+framesPerChunk-1]; % read on this iteration
%    someData = obj.readFrames(framesPerChunk);
    
    if iscell(data)
        someData = obj.readFrames('NumFrames', framesPerChunk);
        for rr = 1:length(someData)
            valuesPerRegion = obj.Regions{rr}.NumYeeCells*length(obj.Fields);
            i1 = (frameNum-1)*valuesPerRegion + 1;
            i2 = i1 + valuesPerRegion*framesPerChunk - 1;
            data{rr}(i1:i2) = someData{rr}(:);
        end
    else
        someData = obj.readFrames('NumFrames', framesPerChunk, 'Regions', ...
            'Together');
        i1 = (frameNum-1)*obj.FrameSize + 1;
        i2 = i1 + obj.FrameSize*framesPerChunk - 1;
        data(i1:i2) = someData(:);
    end
%{    
    if length(obj.Regions) > 1
        for rr = 1:length(obj.Regions)
            valsPerRegion = obj.Regions{rr}.NumYeeCells*length(obj.Fields);
            i1 = (frameNum-1)*valsPerRegion + 1;
            i2 = i1 + valsPerRegion*framesPerChunk - 1;
            data{rr}(i1:i2) = someData{rr}(:);
        end
    else
        valsPerRegion = obj.Regions{1}.NumYeeCells*length(obj.Fields);
        i1 = (frameNum-1)*valsPerRegion + 1;
        i2 = i1 + valsPerRegion*framesPerChunk - 1;
        data(i1:i2) = someData(:);
    end
%}
end

% If there is a remainder after copying by MAXREADBYTES-sized chunks, read it.
if frameRange(2) ~= numFrames
    frameNum = frameRange(2)+1;
    
    if iscell(data)
        someData = obj.readFrames('NumFrames', numFrames - frameRange(2));
        for rr = 1:length(obj.Regions)
            valuesPerRegion = obj.Regions{rr}.NumYeeCells*length(obj.Fields);
            i1 = (frameNum-1)*valuesPerRegion + 1;
            data{rr}(i1:end) = someData{rr}(:);
        end
    else
        someData = obj.readFrames('NumFrames', numFrames - frameRange(2),...
            'Regions', 'Together');
        i1 = (frameNum-1)*obj.FrameSize + 1;
        data(i1:end) = someData(:);
    end
    
    %{
    if iscell(someData)
        for rr = 1:length(obj.Regions)
            valuesPerRegion = obj.Regions{rr}.NumYeeCells*length(obj.Fields);
            i1 = (frameNum-1)*valuesPerRegion + 1;
            data{rr}(i1:end) = someData{rr}(:);
        end
    else
        i1 = (frameNum-1)*obj.FrameSize + 1;
        data(i1:end) = someData(:);
    end
    %}
end

catch exception
    obj.close
    rethrow(exception)
end

obj.close