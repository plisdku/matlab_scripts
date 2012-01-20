function data = readFromFile(obj, numFrames, dataXYZ, outXYZ)

numFields = length(obj.Fields);

% Provided that obj.Precision is 'float32' or 'float64' or something like that,
% this string will have the proper form 'float32=>float32' for fread().
precisionString = [obj.Precision, '=>', obj.Precision];

readLength = obj.FrameSize * numFrames;
[rawdata, count] = fread(obj.FileHandle, readLength, precisionString);

if count ~= numFrames*obj.FrameSize
    error('Could not read %i frames.', numFrames);
end
obj.NextFrameNumber = obj.NextFrameNumber + numFrames;

data = cell([obj.numRegions(), 1]);
for rr = 1:obj.numRegions()
    data{rr} = zeros([obj.Regions.Size(rr,:), numFields, numFrames]);
end

%% Copy rawdata into data{:}

frameFieldVals = obj.FrameSize / numFields;
assert(frameFieldVals == int64(frameFieldVals));

for frameNumber = 1:numFrames
    for fieldNumber = 1:numFields
        readFieldBegin = (fieldNumber-1)*frameFieldVals +...
            (frameNumber-1)*obj.FrameSize + 1;
        
        % At this point, if I read from readFieldBegin for readFieldSize
        % floats, I will read all regions, all x y z, for one field, for
        % one timestep.
        %
        % Instead, go through each region in turn.
        
        for rr = 1:obj.numRegions()
            regionFieldSize = prod(obj.Regions.Size(rr,:));
            read0 = readFieldBegin + obj.RegionOffsetsInFields(rr);
            read1 = read0 + regionFieldSize - 1;
            
            write0 = (fieldNumber-1)*regionFieldSize + ...
                (frameNumber-1)*regionFieldSize*numFields + 1;
            write1 = write0 + regionFieldSize - 1;
            
            data{rr}(write0:write1) = rawdata(read0:read1);
        end
    end
end

%% Interpolate if desired

if nargin == 4
    %fprintf('Interpolating too');
    
    for rr = 1:obj.numRegions()
        outData = zeros(numel(outXYZ{rr,1}), numel(outXYZ{rr,2}), ...
            numel(outXYZ{rr,3}), numFields, size(data{rr},5));
        
        if numel(outData) > 0
        for ff = 1:numFields
            outData(:,:,:,ff,:) = gridInterp(...
                dataXYZ{rr,ff}{1}, dataXYZ{rr,ff}{2}, dataXYZ{rr,ff}{3}, ...
                data{rr}(:,:,:,ff,:), ...
                outXYZ{rr,1}, outXYZ{rr,2}, outXYZ{rr,3});
            %outData(:,:,:,ff,:) = blah;
        end
        end
        data{rr} = outData;
    end
end



