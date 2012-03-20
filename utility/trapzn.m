function A = trapzn(xyzCoordinates, A)
% z = trapzn(x, y) computes the integral of y with respect to x along all its
% non-singleton dimensions using Matlab's built-in trapz() function.
%
%   y   N-D array
%   x   cell array of N elements, where length(x{n}) == size(y, n).
%
% Example: compute the integral of 1 over the x-z plane.  The domain of
% integration is [1 10] x 1 x [1 10].  Since the y-direction is a singleton
% dimension of the integrand, there is no integration in that direction. 
%
%   A = ones(10,1,10);
%   x = {1:10, 1, 1:10};
%   trapzn(x,A)
%
% ans =
%
%    81
%

for xyz = 1:length(xyzCoordinates)
    if numel(xyzCoordinates{xyz}) > 1
        A = trapz(xyzCoordinates{xyz}, A, xyz);
    end
end
