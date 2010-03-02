function data = readFrames(obj, varargin)

if obj.FileHandle == -1
    error('No data file is open.  Try open()?');
end

if nargin > 1
    numFrames = varargin{1};
else
    numFrames = 1;
end

% Some outputs may be rotated x->y->z by a click or two.  Data for a region of
% nominal size [nx ny nz] will first be read into an array of size
% circshift([nx ny nz], [0 -permute]) where permute is defined below.  Then the
% array will be transformed to the right shape as A = shiftdim(A, 3-permute).
if obj.UnitVectors(1,:) == [1 0 0]
    permNum = 0;
elseif obj.UnitVectors(1,:) == [0 1 0]
    permNum = 1;
elseif obj.UnitVectors(1,:) == [0 0 1];
    permNum = 2;
else
    error('Invalid unit vectors.');
end

% There are always two cases:
%   1.  The file contains multiple Regions
%   2.  The file contains one Region
% If there is one Region, then the return value is an array with dimensions
% [cellsX cellsY cellsZ numFields].  If there are multiple regions, then the
% return value is a cell array of arrays, each with size [cellsX ...] etc. as
% above.
% GridReports have a HalfCells attribute instead, so we'll build a table of
% sizes here first.

sizes = {};
if ~isempty(obj.Regions)
    for nn = 1:length(obj.Regions)
        sizes{nn} = obj.Regions{nn}.Size;
    end
    numFields = length(obj.Fields);
elseif ~isempty(obj.HalfCells)
    for nn = 1:length(obj.HalfCells)
        sizes{nn} = obj.HalfCells{nn}.HalfCells(4:6) - ...
            obj.HalfCells{nn}.HalfCells(1:3) + 1;
    end
    numFields = 1;
else
    error('File has no Regions and no HalfCells; probably corrupt.');
end

% Provided that obj.Precision is 'float32' or 'float64' or something like that,
% this string will have the proper form 'float32=>float32' for fread().
precisionString = [obj.Precision, '=>', obj.Precision];

if length(sizes) == 1
    readLength = obj.FrameSize * numFrames;
    data = fread(obj.FileHandle, readLength, precisionString);
    data = reshape(data, [sizes{1}, numFields, numFrames]);
    
elseif length(sizes) > 1
    readLength = obj.FrameSize * numFrames;
    [rawdata, count] = fread(obj.FileHandle, readLength, precisionString);
    
    if count ~= numFrames*obj.FrameSize
        error('Could not read %i frames.', numFrames);
    end
    
    % Resize the data array.  Out of memory errors should happen here and not
    % elsewhere.
    
    data = cell([length(sizes), 1]);
    for rr = 1:length(sizes)
        regionSize = sizes{rr}.NumYeeCells*numFields*numFrames;
        data{rr} = zeros([regionSize, 1]);
%        disp(sprintf('Chunk %i is size %i', rr, length(data{rr})));
    end
    
    % Copy the raw data into the output array.  The hard part is that the file
    % is ordered
    %
    % X Y Z Field Region Timestep
    % 
    % but the output array is ordered
    %
    % Region X Y Z Field Timestep
    % 
    % so I need to rescatter the data into proper locations. Furthermore, some
    % files are not ordered XYZ but maybe YZX or ZXY, which is crudely encoded
    % as 'unitVector0', 'unitVector1', and 'unitVector2' in the output spec
    % file.
    
    readFieldSize = obj.FrameSize / numFields;
    assert(readFieldSize == int64(readFieldSize));
    
    for frameNumber = 1:numFrames
%        disp(sprintf('Frame %i', frameNumber))
        for fieldNumber = 1:numFields
            readFieldBegin = (fieldNumber-1)*readFieldSize + ...
                (frameNumber-1)*obj.FrameSize + 1;
            for rr = 1:length(sizes)
                regionFieldSize = prod(sizes{rr});
                
                readFirst = readFieldBegin + obj.RegionOffsetsInFields(rr);
                readLast = readFirst + regionFieldSize - 1;
                
%                disp(sprintf('Read from %i to %i for region %i field %i', ...
%                    readFirst, readLast, rr, fieldNumber));
                
                writeFirst = (fieldNumber-1)*regionFieldSize + ...
                    (frameNumber-1)*regionFieldSize*numFields + 1;
                writeLast = writeFirst + regionFieldSize - 1;
                
                data{rr}(writeFirst:writeLast) = rawdata(readFirst:readLast);
            end
        end
    end
    
    % Reshape the data cells to [X Y Z Field Timestep].
    % This has two steps:
    %   1.  Reshape the data to its permuted size, based on UnitVectors.
    %   2.  Shiftdim the data into XYZ order again.
    % In the file, each region has XYZ size regionSizeWritten, and we'll
    % permute it back to regionSizeHere.
    for rr = 1:length(sizes)
        regionSizeWritten = circshift(sizes{rr}, [0, -permNum]);
        regionSizeHere = sizes{rr};
        data{rr} = reshape(data{rr}, ...
            [regionSizeWritten, numFields, numFrames]);
        permVec = [circshift(1:3, [0, permNum]), 4:ndims(data{rr})];
        data{rr} = permute(data{rr}, permVec);
    end
    
else
    error('OutputFile has %d regions, how weird.', length(sizes));
end

