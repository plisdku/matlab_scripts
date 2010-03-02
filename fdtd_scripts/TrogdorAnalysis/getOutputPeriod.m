function period = getOutputPeriod(fileprefix)
% period = getOutputPeriod(filePrefix) returns the sampling period of a
%   Trogdor output, if applicable.  The default period is 1, i.e. one
%   sample per timestep.
%
%   See also: getOutputDimensions, getOutputBounds
%
%   version 4.5
%   July 29, 2008


specfile = [fileprefix, '.txt'];

specFID = fopen(specfile, 'rt');
if (specFID == -1)
    error('Could not open spec file %s.', specfile);
end

done = 0;

while done == 0
    line = fgetl(specFID);
    if ~ischar(line)
        break
    else
        [token, remainder] = strtok(line);
        if strcmp(token, 'period')
            [stuff, count] = sscanf(remainder, '%f');
            if count == 1
                period = stuff;
                done = 1;
            end
        end
    end
end

if ~done
    error('Could not find bounds in spec file.');
end
