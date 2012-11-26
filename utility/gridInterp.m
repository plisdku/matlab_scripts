function B = gridInterp(varargin)
%gridInterp Faster, memory-efficient interpn for regular cartesian grids
%   v = gridInterp(X1,X2,X3,...,A,Y1,Y2,Y3,...) interpolates to find v, the
%   values of the underlying N-D function A at all points (y1,y2,y3,...)
%   where y1 is in Y1, y2 is in Y2, etc.  Arrays X1, X2, etc. are 1D arrays
%   of coordinates at which the data in A are given.  Arrays Y1, Y2, etc.
%   are 1D arrays of coordinates at which the data in A is desired.
%
%   If Xn for some n are empty arrays, then no interpolation
%   will be carried out in those dimensions.  This provides a way to select
%   which dimensions to interpolate.

numdim = floor(nargin/2);
if numdim > 10
    % The reason for this is that I wrote a cell array indicesAll
    % consisting of ten ':' strings.  It's much faster than using repmat on
    % cells, or any of my other strategies.
    error('gridInterp can only handle up to 10 dimensions');
end

if nargin ~= 2*numdim + 1
    error('Wrong number of arguments');
end

x = varargin(1:numdim);
A = varargin{numdim+1};
y = varargin(numdim+2:end);

%C = A;
%D = A;

for nn = 1:numdim
    if ~isempty(x{nn})
    if size(A,nn) ~= length(x{nn})
        error('Length of X%i must equal size(A,%i)', nn, nn');
    end
    end
end

B = A;

% This shaves a lot of time off of calls with SMALL arrays going in.
%indicesAll = repmat({':'}, 1, numdim); % useful variable factored from loop
indicesAll = {':',':',':',':',':',':',':',':',':',':'};

for dim = 1:numdim
if size(A,dim) > 1 && ~isempty(x{dim})
    xIn = x{dim};
    xOut = y{dim};
    
    % Clamp the output positions to not exceed the provided input
    % positions.  Extrapolation to outside values will consequently be by
    % clamping of course.  I don't need anything else right now.
    xOut(xOut < xIn(1)) = xIn(1);
    xOut(xOut > xIn(end)) = xIn(end);
    
    dx = (xIn(end)-xIn(1))/(length(xIn)-1);
    
    leftCell = 1+floor( (xOut-xIn(1))/dx ); % these are all valid indices
    rightCell = min(leftCell+1, size(A,dim)); % clamp to valid indices
    
    distLeft = xOut - xIn(leftCell);
    distRight = dx - distLeft;
    
    szWeight = ones(1,ndims(A)); szWeight(dim) = length(distRight);
    weightLeft = zeros(szWeight);
    weightLeft(:) = distRight(:)/dx;
    weightRight = 1 - weightLeft;
    
    indicesLeft = indicesAll;
    indicesRight = indicesAll;
    
    indicesLeft{dim} = leftCell;
    indicesRight{dim} = rightCell;
    
    B = bsxfun(@times, weightLeft, B(indicesLeft{:})) + ...
        bsxfun(@times, weightRight, B(indicesRight{:}));
end
end


