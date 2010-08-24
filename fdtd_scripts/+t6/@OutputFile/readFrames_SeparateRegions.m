function data = readFrames_SeparateRegions(obj, numFrames)

% There are always two cases:
%   1.  The file contains multiple Regions
%   2.  The file contains one Region
% If there is one Region, then the return value is an array with dimensions
% [cellsX cellsY cellsZ numFields].  If there are multiple regions, then the
% return value is a cell array of arrays, each with size [cellsX ...] etc. as
% above.

sizes = {};
if ~isempty(obj.Regions)
    for nn = 1:length(obj.Regions)
        sizes{nn} = obj.Regions{nn}.Size;
    end
    numFields = length(obj.Fields);
else
    error('File has no Regions.');
end

% Provided that obj.Precision is 'float32' or 'float64' or something like that,
% this string will have the proper form 'float32=>float32' for fread().
precisionString = [obj.Precision, '=>', obj.Precision];

%if length(sizes) == 1
if length(obj.Regions) == 1
    readLength = obj.FrameSize * numFrames;
    data = fread(obj.FileHandle, readLength, precisionString);
    data = reshape(data, [sizes{1}, numFields, numFrames]);
    
elseif length(obj.Regions) > 1
    readLength = obj.FrameSize * numFrames;
    [rawdata, count] = fread(obj.FileHandle, readLength, precisionString);
    
    if count ~= numFrames*obj.FrameSize
        error('Could not read %i frames.', numFrames);
    end
    
    % Resize the data array.  Out of memory errors should happen here and not
    % elsewhere.
    data = cell([length(obj.Regions), 1]);
    for rr = 1:length(obj.Regions)
        regionSize = prod(sizes{rr})*numFields*numFrames;
        data{rr} = zeros([regionSize, 1]);
%        disp(sprintf('Chunk %i is size %i', rr, length(data{rr})));
    end
    
    % Copy the raw data into the output array.  The hard part is that the file
    % is ordered
    %
    % (TROGDOR 5) X Y Z Field Region Timestep
    % (TROGDOR 6) X Y Z Region Field Timestep
    % 
    % but the output array is ordered
    %
    % Region X Y Z Field Timestep
    % 
    % so I need to rescatter the data into proper locations.
    
    readFieldSize = obj.FrameSize / numFields; % like prod(x y z regions)
    assert(readFieldSize == int64(readFieldSize));
    
    for frameNumber = 1:numFrames
%        disp(sprintf('Frame %i', frameNumber))
        for fieldNumber = 1:numFields
            readFieldBegin = (fieldNumber-1)*readFieldSize + ...
                (frameNumber-1)*obj.FrameSize + 1;
            
            % At this point, if I read from readFieldBegin for readFieldSize
            % floats, I will read all regions, all x y z, for one field, for
            % one timestep.
            %
            % Instead, go through each region in turn.
            
            for rr = 1:length(obj.Regions)
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
    for rr = 1:length(obj.Regions)
        %regionSizeWritten = circshift(sizes{rr}, [0, -permNum]);
        regionSizeWritten = sizes{rr};
        regionSizeHere = sizes{rr};
        data{rr} = reshape(data{rr}, ...
            [regionSizeWritten, numFields, numFrames]);
        %permVec = [circshift(1:3, [0, permNum]), 4:ndims(data{rr})];
        %data{rr} = permute(data{rr}, permVec);
    end
    
else
    error('OutputFile has %d regions, how weird.', length(sizes));
end

