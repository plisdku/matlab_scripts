function data = read(obj, varargin)
% data = read(obj, varargin)
%
% Reads every field in the file.
% 
% Named options:
%
%   'Separate'      Put separate regions in separate cells
%

X.Regions = 'Separate';
X = parseargs(X, varargin{:});

% What does this do?
% It reads all the fields!  All of them!  A simple but slow approach
% would use readFrames().  I'll try that first.

obj.open
try

numFrames = obj.numFramesAvailable;

MAXREADBYTES = 1e8;   % Read up to 100 MB at a time, to reduce memory overhead.
bytesPerFrame = obj.FrameSize * obj.BytesPerValue;
framesPerChunk = min(numFrames, ceil(MAXREADBYTES / bytesPerFrame) );
numChunks = floor(numFrames / framesPerChunk);

% Reserve some space
% size: [x y z f t r]
if obj.numRegions() > 1 && strcmp(X.Regions, 'Separate')
    data = cell(obj.numRegions(),1);
    for rr = 1:obj.numRegions()
        data{rr} = zeros([obj.Regions.Size(rr,:), length(obj.Fields), ...
            numFrames]);
    end
elseif obj.numRegions() == 1 && strcmp(X.Regions, 'Separate')
    data = zeros([obj.Regions.Size(1,:), length(obj.Fields), numFrames]);
elseif ~strcmp(X.Regions, 'Separate')
    data = zeros(obj.FrameSize/length(obj.Fields), length(obj.Fields), numFrames);
else
    error('No Regions available.  Read GridReports with readFrames(1).');
end

warning('I think I can unify the two paragraphs of chunk reading');
% Copy in the data by chunks as far as possible
for chunk = 1:numChunks
    frameNum = (chunk-1)*framesPerChunk + 1;
    frameRange = [frameNum frameNum+framesPerChunk-1]; % read on this iteration
%    someData = obj.readFrames(framesPerChunk);
    
    if iscell(data)
        someData = obj.readFrames('NumFrames', framesPerChunk);
        for rr = 1:length(someData)
            valuesPerRegion = obj.Regions.NumYeeCells(rr)*length(obj.Fields);
            i1 = (frameNum-1)*valuesPerRegion + 1; % start of region
            i2 = i1 + valuesPerRegion*framesPerChunk - 1; % end of region
            data{rr}(i1:i2) = someData{rr}(:);
        end
    else
        someData = obj.readFrames('NumFrames', framesPerChunk, 'Regions', ...
            'Together');
        i1 = (frameNum-1)*obj.FrameSize + 1; % start of frame
        i2 = i1 + obj.FrameSize*framesPerChunk - 1; % end of frame
        data(i1:i2) = someData(:);
    end
end

% If there is a remainder after copying by MAXREADBYTES-sized chunks, read it.
if frameRange(2) ~= numFrames
    frameNum = frameRange(2)+1;
    
    if iscell(data)
        someData = obj.readFrames('NumFrames', numFrames - frameRange(2));
        for rr = 1:obj.numRegions()
            valuesPerRegion = obj.Regions.NumYeeCells(rr,:)*length(obj.Fields);
            i1 = (frameNum-1)*valuesPerRegion + 1;
            data{rr}(i1:end) = someData{rr}(:);
        end
    else
        someData = obj.readFrames('NumFrames', numFrames - frameRange(2),...
            'Regions', 'Together');
        i1 = (frameNum-1)*obj.FrameSize + 1;
        data(i1:end) = someData(:);
    end
    
end

catch exception
    obj.close
    rethrow(exception)
end

obj.close