function data = readFrames_RegionsSeparate_InterpTime(obj, X)
% X.NumFrames
% X.Regions
% X.Size
% X.Times
% X.Positions

numFields = obj.numFields();

if ~isempty(X.Positions)
    assert(iscell(X.Positions));
    assert(all(size(X.Positions) == [obj.numRegions(), 3]));
    
    samplePositionsInFile = cell(obj.numRegions(), obj.numFields());
    for rr = 1:obj.numRegions
        for ff = 1:obj.numFields
            samplePositionsInFile{rr,ff} = obj.positions('InterpolateSpace', ...
                false);
        end
    end
    
    subReadData = @(obj, numFrames) readFromFile(obj, numFrames, ...
        samplePositionsInFile, X.Positions);
else
    subReadData = @(obj, numFrames) readFromFile(obj, numFrames);
end

% The n2t function takes a frame number and provides a time.

if numFields == 1
    sampleTimes = [-1, obj.times()];
    n2t = @(n,f) sampleTimes(n+1);
else
    sampleTimes = obj.times();
    for ff = 1:numFields
        sampleTimes{ff} = [-1, sampleTimes{ff}];
    end
    n2t = @(n,f) sampleTimes{f}(n+1);
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
[chunkStarts, chunkEnds] = chunkTimesteps(nFirst, nLast, frameSizes, chunkBytes);
numChunks = length(chunkStarts);
chunkSize = chunkEnds(1)-chunkStarts(1)+1;
bufferLength = chunkSize+1;

% 12/16/11: note that only firstSampleInChunk and lastSampleInChunk
% actually seem to get used, out of the whole fields struct array.
fields = struct;
for ff = 1:numFields
    fields(ff).tChunkStart = n2t(chunkStarts, ff);
    fields(ff).tChunkEnd = n2t(chunkStarts, ff);
    [unused, fields(ff).whichChunks, fields(ff).whichTimes] = ...
        binSamples(X.Times, n2t([chunkStarts(1), chunkEnds], ff));
    
    fields(ff).firstSampleInChunk = fields(ff).whichTimes(1:end-1);
    fields(ff).lastSampleInChunk = fields(ff).whichTimes(2:end)-1;
    
    if any(fields(ff).whichChunks == 0)
        error('Not all times fall into measured duration');
    end
end

%fakeData = rand([frameSizes, nLast-nFirst+1]);

% Multi-region
data = cell(obj.numRegions(), 1);
for rr = 1:obj.numRegions()
    data{rr} = zeros([frameSizes(rr,:), numel(X.Times)]);
end
% Single-region
%data = zeros([frameSizes, numel(X.Times)]);

%%

% Multi-region
buffer.data = cell(obj.numRegions(), 1);
for rr = 1:obj.numRegions()
    buffer.data{rr} = zeros([frameSizes(rr,:), bufferLength]);
end
% Single-region
%assert(size(frameSizes, 1) == 1);
%buffer.data = zeros([frameSizes bufferLength]);

buffer.frameNumbers = 0*(1:bufferLength);

if numel(obj.SavedFrame) == numel(buffer.data)
    % Multi-region
    for rr = 1:obj.numRegions()
        buffer.data{rr}(:,:,:,:,1) = obj.SavedFrame{rr};
    end
    % Single-region
    %buffer.data(:,:,:,:,1) = obj.SavedFrame;
    buffer.frameNumbers(1) = obj.SavedFrameNumber;
end

for cc = 1:numChunks
    n0 = chunkStarts(cc);
    n1 = chunkEnds(cc);
    chunkLength = n1-n0+1;
    
    % modify here to read multiple frames.
    % it looks like i need buffer.data to be a cell array.  but then i
    % have to copy the data region by region.  that's sloooowwww.
    % (Is it slow?  Measure first.)
    
    % Multi-region
    dataFromFile = subReadData(obj, n1-n0+1);
    for rr = 1:obj.numRegions()
        buffer.data{rr}(:,:,:,:,2:1+chunkLength) = dataFromFile{rr};
    end
    % Single-region
    %cellElem = @(A) A{1};
    %buffer.data(:,:,:,:,2:1+chunkLength) = cellElem(readFromFile(obj, n1-n0+1));
    
    buffer.frameNumbers(2:1+chunkLength) = n0:n1;
    
    for ff = 1:numFields
        bufferTimes = n2t(buffer.frameNumbers, ff);
        
        %fprintf('Field %i, chunk %i (n = %i:%i)', ff, cc, n0, n1);
        %fprintf('\tBuffer times %2.2f to %2.2f\n', bufferTimes([1 end]));
        
        tOutFromChunk = X.Times(fields(ff).firstSampleInChunk(cc):fields(ff).lastSampleInChunk(cc));
        
        if ~isempty(tOutFromChunk)
            %fprintf('\tThere are %i outputs in this chunk:', length(tOutFromChunk));
            %fprintf(' %2.2f', tOutFromChunk); fprintf('\n');
            
            [N, bins, indices] = binSamples(tOutFromChunk, ...
                bufferTimes(1:chunkLength+1));
            
            if 0
                occupiedBins = find(N);
                for ww = occupiedBins
                    sampleTimes = tOutFromChunk(indices(ww):indices(ww+1)-1);
                    t0 = bufferTimes(ww);
                    t1 = bufferTimes(ww+1);
                    fprintf('\tTimes ');
                    fprintf('%2.2f ', sampleTimes);
                    fprintf(' between %2.2f and %2.2f\n', t0, t1);
                end
            end
            
            distLeft = tOutFromChunk - bufferTimes(bins);
            distRight = bufferTimes(bins+1) - tOutFromChunk;

            weightLeft = reshape(distRight ./ (distRight + distLeft), ...
                [ones(size(frameSizes(1,:))), length(distRight)]);
            weightRight = 1 - weightLeft;
            
            %fprintf('Weights %2.2f, %2.2f\n', weightLeft, weightRight);
            
            % modify here to handle multiple regions.
            % this will have to be a loop over the buffer cells and dataOut
            % cells.
            
            nOut = fields(ff).firstSampleInChunk(cc):fields(ff).lastSampleInChunk(cc);
            
            % Multi-region
            for rr = 1:obj.numRegions()
                data{rr}(:,:,:,ff,nOut) = ...
                    bsxfun(@times, weightLeft, buffer.data{rr}(:,:,:,ff,bins)) + ...
                    bsxfun(@times, weightRight, buffer.data{rr}(:,:,:,ff,bins+1));
            end
            % Single-region
            %dataOut = bsxfun(@times, weightLeft, buffer.data(:,:,:,ff,bins)) + ...
            %    bsxfun(@times, weightRight, buffer.data(:,:,:,ff,bins+1));
            %data(:,:,:,ff,nOut) = dataOut;
            
        end
    end
    
    % Multi-region
    for rr = 1:obj.numRegions()
        buffer.data{rr}(:,:,:,:,1) = buffer.data{rr}(:,:,:,:,1+chunkLength);
    end
    % Single-region
    %buffer.data(:,:,:,:,1) = buffer.data(:,:,:,:,1+chunkLength);
    buffer.frameNumbers(1) = buffer.frameNumbers(1+chunkLength);
end

%% Keep the extra frame from the buffer
% Multi-region
for rr = 1:obj.numRegions()
    obj.SavedFrame{rr} = buffer.data{rr}(:,:,:,:,1+chunkLength);
end
% Single-region
%obj.SavedFrame = buffer.data(:,:,:,:,1+chunkLength);

obj.SavedFrameNumber = buffer.frameNumbers(1+chunkLength);


%fprintf('Saving %i\n', obj.SavedFrameNumber);

%% Single-region outputs don't return cell arrays.  (Right?)
if obj.numRegions() == 1
    data = data{1};
end




