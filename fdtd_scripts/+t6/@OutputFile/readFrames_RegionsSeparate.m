function data = readFrames_RegionsSeparate(obj, X)
% X.NumFrames
% X.Regions
% X.Size
% X.Times

assert(isempty(X.Times));

numFields = obj.numFields();

if ~isempty(X.Positions)
    assert(iscell(X.Positions));
    assert(all(size(X.Positions) == [obj.numRegions(), 3]));
    
    % Here I need to cache the region & field positions as present in the
    % file (raw Yee positions, for instance).
    
    samplePositionsInFile = cell(obj.numRegions(), obj.numFields());
    for rr = 1:obj.numRegions
        for ff = 1:obj.numFields
            samplePositionsInFile{rr,ff} = obj.positions(...
                'Field', ff, 'InterpolateSpace', false);
        end
    end
    
    subReadData = @(obj, numFrames) readFromFile(obj, numFrames, ...
        samplePositionsInFile, X.Positions);
else
    subReadData = @(obj, numFrames) readFromFile(obj, numFrames);
end

% We read timesteps nFirst to nLast.
% The data available from these reads span n2t(nFirst) to n2t(nLast).
nFirst = obj.NextFrameNumber;
nLast = nFirst + X.NumFrames - 1;

frameSizes = obj.Regions.Size;
frameSizes(:,4) = numFields;

if ~isempty(X.Positions)
    for rr = 1:obj.numRegions()
    for xyz = 1:3
        frameSizes(rr,xyz) = numel(X.Positions{rr,xyz});
    end
    end
end

chunkBytes = 10;
[chunkStarts, chunkEnds] = t6.OutputFile.chunkTimesteps(nFirst, nLast, frameSizes, chunkBytes);
numChunks = length(chunkStarts);

% Multi-region
data = cell(obj.numRegions(), 1);
for rr = 1:obj.numRegions()
    data{rr} = zeros([frameSizes(rr,:), X.NumFrames]);
end
% Single-region
%data = zeros([frameSizes, numel(X.Times)]);

%% Read chunks one at a time and write into output array correctly

for cc = 1:numChunks
    n0 = chunkStarts(cc);
    n1 = chunkEnds(cc);
    
    % modify here to read multiple frames.
    % it looks like i need buffer.data to be a cell array.  but then i
    % have to copy the data region by region.  that's sloooowwww.
    % (Is it slow?  Measure first.)
    
    % Multi-region
    dataFromFile = subReadData(obj, n1-n0+1);
    for rr = 1:obj.numRegions()
        data{rr}(:,:,:,:,1+(n0:n1)-nFirst) = dataFromFile{rr};
        %fprintf('Writing %i to %i into %i to %i\n', n0, n1, ...
        %    1+n0-n0, 1+n1-n0);
    end
end

%% Keep the extra frame from the buffer
% Multi-region
for rr = 1:obj.numRegions()
    obj.SavedFrame{rr} = data{rr}(:,:,:,:,end);
end

obj.SavedFrameNumber = n1;

%% Single-region outputs don't return cell arrays.
if obj.numRegions() == 1
    data = data{1};
end




