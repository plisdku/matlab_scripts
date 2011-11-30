% readSpecFile
function readSpecFile(obj)

isTrogdor5 = 0;


% Figure out how many regions there will be, in advance
numRegions = 0;
numDurations = 0;

fid = fopen(obj.SpecFileName, 'r');
try
    done = 0;
    while done == 0
        lineFromFile = fgets(fid);
        if ~ischar(lineFromFile)
            done = 1;
        else
            [token, remainder] = strtok(lineFromFile);
            switch token
                case 'region'
                    numRegions = numRegions + 1;
                case 'duration'
                    numDurations = numDurations + 1;
                otherwise
                    % nothing
            end
        end
    end
catch err
    error('How does this fail?');
end
fclose(fid);

%fprintf('%i regions, %i durations\n', numRegions, numDurations);

fid = fopen(obj.SpecFileName, 'r');

linesRead = 0;

obj.Regions.YeeCells = zeros(numRegions, 6);
obj.Regions.Size = zeros(numRegions, 3);
obj.Regions.Stride = zeros(numRegions, 3);
obj.Regions.NumYeeCells = zeros(numRegions, 1);
obj.Regions.Bounds = nan(numRegions, 6);

obj.Durations = cell(numDurations, 1);
nRegion = 1;
nDuration = 1;

%t0 = cputime;

try % everything else, but close file before rethrow
    done = 0;
    while (done == 0)
        linesRead = linesRead + 1;
        %fprintf('On line %i\n', linesRead);
        %if mod(linesRead, 100) == 0
        %    fprintf('%2.2f lines/time\n', linesRead / (cputime - t0));
        %end
        
        %if linesRead == 8000
        %    done = 1;
        %end
        
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
                [dat, count] = sscanf(remainder, ' [%f, %f, %f]');
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
            case 'origin'
                [dat, count] = sscanf(remainder, ' [%f, %f, %f]');
                if count == 3
                    obj.Origin = dat';
                else
                    error('Cannot parse line %s', lineFromFile);
                end
            case 'field'
                [fieldName, numbers] = strtok(remainder);
                [dat, count] = sscanf(numbers, ' [%f, %f, %f] %f');
                if count == 4
                    field = [];
                    field.Name = fieldName;
                    field.Offset = dat';                    
                    obj.Fields = [obj.Fields {field}];
                else
                    error('Cannot parse line %s', lineFromFile);
                end
            case 'unitVector0'
                [dat, count] = sscanf(remainder, ' [%f, %f, %f]');
                if count == 3
                    obj.UnitVectors(1,:) = dat';
                else
                    error('Cannot parse line %s', lineFromFile);
                end
            case 'unitVector1'
                [dat, count] = sscanf(remainder, ' [%f, %f, %f]');
                if count == 3
                    obj.UnitVectors(2,:) = dat';
                else
                    error('Cannot parse line %s', lineFromFile);
                end
            case 'unitVector2'
                [dat, count] = sscanf(remainder, ' [%f, %f, %f]');
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
                    ' [[%f, %f, %f], [%f, %f, %f]] stride [%f, %f, %f] bounds [[%f, %f, %f], [%f, %f, %f]]');
                if count == 9
                    obj.Regions.YeeCells(nRegion,:) = dat(1:6)';
                    obj.Regions.Size(nRegion,:) = ceil( (dat(4:6)-dat(1:3)+1) ./ dat(7:9))';
                    obj.Regions.Stride(nRegion,:) = dat(7:9)';
                    obj.Regions.NumYeeCells(nRegion) = prod(obj.Regions.Size(nRegion,:));
                    nRegion = nRegion + 1;
                elseif count == 15
                    obj.Regions.YeeCells(nRegion,:) = dat(1:6)';
                    obj.Regions.Size(nRegion,:) = ceil( (dat(4:6)-dat(1:3)+1) ./ dat(7:9))';
                    obj.Regions.Stride(nRegion,:) = dat(7:9)';
                    obj.Regions.NumYeeCells(nRegion) = prod(obj.Regions.Size(nRegion,:));
                    obj.Regions.Bounds(nRegion,:) = dat(10:15)';
                    nRegion = nRegion + 1;
                else
                    count
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
                    obj.Durations{nDuration} = duration;
                    nDuration = nDuration + 1;
                    %obj.Durations = {obj.Durations{:}, duration};
                else
                    error('Cannot parse line %s', lineFromFile);
                end
            % For grid reports
            case 'material'
                [materialNumber, materialName] = strtok(remainder);
                obj.Materials{str2num(materialNumber)+1} = materialName;
            case 'halfCells'
                [dat, count] = sscanf(remainder, ...
                    ' [[%f, %f, %f], [%f, %f, %f]]');
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

obj.RegionOffsetsInFields = zeros(numRegions,1);

totalYeeCells = 0;
regionOffset = 0;
for nn = 1:numRegions
    totalYeeCells = totalYeeCells + obj.Regions.NumYeeCells(nn,:);
    obj.RegionOffsetsInFields(nn) = regionOffset;
    regionOffset = regionOffset + obj.Regions.NumYeeCells(nn,:);
end
obj.FrameSize = totalYeeCells * length(obj.Fields);

fieldOffset = 0;
for nn = 1:length(obj.Fields)
    obj.FieldOffsetsInFrames(nn) = fieldOffset;
    fieldOffset = fieldOffset + totalYeeCells;
end

    
