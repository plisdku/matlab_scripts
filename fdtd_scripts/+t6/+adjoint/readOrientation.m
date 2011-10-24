function o = readOrientation(filePrefix)
% readOrientation    Read orientation matrices from file
%
% Usage:
% o = readOrientation(filePrefix);
%
% size(o) = [numX numY numZ orientation_i orientation_j eps_i eps_j]
%
% Typically the orientation report will contain the orientation matrices
% for Ex, Ey and Ez in order (the fieldXYZ index).

header = [filePrefix, '.txt'];

fid = fopen(header, 'r');

precision = [];
dims = [];

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
            case 'dims'
                [dat, count] = sscanf(remainder, ...
                    ' %i %i %i %i %i %i');
                if count == 7
                    dims = dat;
                elseif count == 6
                    dims = dat(1:6);
                else
                    error('Cannot parse line %s', lineFromFile);
                end
            otherwise
                %nothing
        end
    end
end
catch err
    error(err.message);
end

fclose(fid);

if isempty(precision) || isempty(dims)
    error('Something is missing.');
end

fid = fopen(filePrefix);

data = fread(fid, inf, precision);
fclose(fid);

o = reshape(data, dims');


