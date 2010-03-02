function bounds = getOutputBounds(fileprefix)
% bounds = getOutputBounds(filePrefix) returns the 3D bounds of a Trogdor
%   output region, if applicable.  The bounds will be in the standard
%   Trogdor bounding box format,
%       [x1 y1 z1 x2 y2 z2]
%   and can be plotted directly with trogRect, e.g.
%       trogRect(getOutputBounds('outEField'));
%
%   See also: getOutputDimensions, getOutputPeriod, trogRect
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
        if strcmp(token, 'bounds')
            [stuff, count] = sscanf(remainder, '%f %f %f %f %f %f');
            if count == 6
                bounds = stuff;
                done = 1;
            end
        end
    end
end

if ~done
    error('Could not find bounds in spec file.');
end
