function [numer, denom] = readPermittivity(filePrefix)
% readPermittivity    Load numerator and denominator data for a whole grid
%
% Usage:
%
% [numer, denom] = readPermittivity(filePrefix)
%
% e.g. readPermittivity('permittivity');
%
% size(numer) = [numX numY numZ 3 numLags]
% size(denom) = [numX numY numZ 3 numLags]

header = [filePrefix, '.txt'];

fid = fopen(header, 'r');

precision = [];
yeeCells = [];
yeeSize = [];
stride = [];
numYeeCells = [];
numerOrder = [];
denomOrder = [];

try
done = 0;
while (done == 0)
    lineFromFile = fgets(fid);
    if ~ischar(lineFromFile)
        done = 1;
    else
        [token, remainder] = strtok(lineFromFile);
        switch token
            case 'precision'
                precision = sscanf(remainder, '%s');
            case 'region'
                [dat, count] = sscanf(remainder, ...
                    ' [[%f, %f, %f], [%f, %f, %f]] stride [%f, %f, %f]');
                if count == 9
                    yeeCells = dat(1:6)';
                    yeeSize = ceil( (dat(4:6)-dat(1:3)+1) ./ dat(7:9) )';
                    stride = dat(7:9)';
                    numYeeCells = prod(yeeSize);
                else
                    error('Cannot parse line %s', lineFromFile);
                end
            case 'numerator'
                numerOrder = sscanf(remainder, ' %i lags');
            case 'denominator'
                denomOrder = sscanf(remainder, ' %i lags');
            otherwise
                %nothing
        end
    end
end
catch err
    error('Flail.');
end

fclose(fid);

if isempty(precision) || isempty(yeeCells) || isempty(yeeSize) ...
    || isempty(stride) || isempty(numYeeCells) || isempty(numerOrder)...
    || isempty(denomOrder)
    error('Something is missing.');
end

dims = [yeeSize 3 3 numerOrder+denomOrder];

fid = fopen(filePrefix);

data = fread(fid, inf, precision);
fclose(fid);

p = reshape(data, dims);

numer = p(:,:,:,:,:,1:numerOrder);
denom = p(:,:,:,:,:,numerOrder+1:end);

