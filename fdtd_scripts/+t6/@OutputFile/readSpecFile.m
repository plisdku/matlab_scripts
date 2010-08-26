% readSpecFile
function readSpecFile(obj)

isTrogdor5 = 0;

fid = fopen(obj.SpecFileName, 'r');

try % everything else, but close file before rethrow
    done = 0;
    while (done == 0)
        lineFromFile = fgets(fid);
        if ~ischar(lineFromFile)
            %disp('Done.');
            done = 1;
        else
           [token, remainder] = strtok(lineFromFile);
           
            switch token
            case 'trogdor5data'
                isTrogdor5 = 1;
            case 'trogdorVersionNumber'
                obj.TrogdorVersionString = remainder;
            case 'trogdorBuildDate'
            case 'specfile'
            case 'datafile'
            case 'materialfile'
            case 'precision'
                precision = sscanf(remainder, '%s');
                if ~strcmp(precision, 'float32') && ...
                    ~strcmp(precision, 'float64')
                    error(sprintf('Invalid precision %s', precision));
                end
                obj.Precision = precision;
                
                if strcmp(precision, 'float32')
                    obj.BytesPerValue = 4;
                elseif strcmp(precision, 'float64')
                    obj.BytesPerValue = 8;
                end
                                
            case 'date'
                obj.DateString = remainder;
            case 'dxyz'
                [dat, count] = sscanf(remainder, ' (%f, %f, %f)');
                if count == 3
                    obj.Dxyz = dat';
                else
                    error('Cannot parse line %s', lineFromFile);
                end
            case 'dt'
                [dat, count] = sscanf(remainder, ' %f');
                if count == 1
                    obj.Dt = dat;
                else
                    error('Cannot parse line %s', lineFromFile);
                end
            case 'field'
                [fieldName, numbers] = strtok(remainder);
                [dat, count] = sscanf(numbers, ' (%f, %f, %f) %f');
                if count == 4
                    field = [];
                    field.Name = fieldName;
                    field.Offset = dat';                    
                    obj.Fields = {obj.Fields{:}, field};
                else
                    error('Cannot parse line %s', lineFromFile);
                end
            case 'unitVector0'
                [dat, count] = sscanf(remainder, ' (%f, %f, %f)');
                if count == 3
                    obj.UnitVectors(1,:) = dat';
                else
                    error('Cannot parse line %s', lineFromFile);
                end
            case 'unitVector1'
                [dat, count] = sscanf(remainder, ' (%f, %f, %f)');
                if count == 3
                    obj.UnitVectors(2,:) = dat';
                else
                    error('Cannot parse line %s', lineFromFile);
                end
            case 'unitVector2'
                [dat, count] = sscanf(remainder, ' (%f, %f, %f)');
                if count == 3
                    obj.UnitVectors(3,:) = dat';
                else
                    error('Cannot parse line %s', lineFromFile);
                end
            case 'runlineDirection'
                [dat, count] = sscanf(remainder, ' %f');
                if count == 1
                    obj.RunlineDirection = dat;
                else
                    error('Cannot parse line %s', lineFromFile);
                end
            % For field outputs
            case 'region'
                [dat, count] = sscanf(remainder, ...
                    ' [(%f, %f, %f), (%f, %f, %f)] stride (%f, %f, %f)');
                if count == 9
                    region = struct;
                    region.YeeCells = dat(1:6)';
                    region.Size = ceil( (dat(4:6)-dat(1:3)+1) ./ dat(7:9) )';
                    region.Stride = dat(7:9)';
                    region.NumYeeCells = prod(region.Size);
                    obj.Regions = {obj.Regions{:}, region};
                else
                    error('Cannot parse line %s', lineFromFile);
                end
            case 'duration'
                [dat, count] = sscanf(remainder, ' from %f to %f period %f');
                if count == 3
                    duration = struct;
                    duration.First = dat(1);
                    duration.Last = dat(2);
                    duration.Period = dat(3);
                    duration.NumTimesteps = floor( (dat(2)-dat(1)+1)/dat(3) );
                    obj.Durations = {obj.Durations{:}, duration};
                else
                    error('Cannot parse line %s', lineFromFile);
                end
            % For grid reports
            case 'material'
                [materialNumber, materialName] = strtok(remainder);
                obj.Materials{str2num(materialNumber)+1} = materialName;
            case 'halfCells'
                [dat, count] = sscanf(remainder, ...
                    ' [(%f, %f, %f), (%f, %f, %f)]');
                if count == 6
                    halfCells = struct;
                    halfCells.HalfCells = dat';
                    cellSpan = dat(4:6)-dat(1:3)+1;
                    halfCells.NumHalfCells = prod(cellSpan);
                    obj.HalfCells = {obj.HalfCells{:}, halfCells};
                else
                    error('Cannot parse line %s', lineFromFile);
                end
            otherwise
                % nothing.
            end
        end
    end
catch exception
    fclose(fid)
    rethrow(exception)
end

fclose(fid);

% Now work on the hidden properties: obj.RegionOffsetsInFields,
% and obj.FieldOffsetsInFrames.

totalYeeCells = 0;
regionOffset = 0;
for nn = 1:length(obj.Regions)
    totalYeeCells = totalYeeCells + obj.Regions{nn}.NumYeeCells;
    obj.RegionOffsetsInFields(nn) = regionOffset;
    regionOffset = regionOffset + obj.Regions{nn}.NumYeeCells;
end
obj.FrameSize = totalYeeCells * length(obj.Fields);

fieldOffset = 0;
for nn = 1:length(obj.Fields)
    obj.FieldOffsetsInFrames(nn) = fieldOffset;
    fieldOffset = fieldOffset + totalYeeCells;
end

    
