function dt = courant(varargin)
%courant Return the maximum timestep under the Courant-Friedrichs-Lewy condition
%   dt = courant(dx) returns dx/c
%
%   dt = courant(dx, dy)
%   dt = courant([dx dy]) returns dx/c/sqrt(1/dx^2 + 1/dy^2)
%
%   dt = courant(dx, dy, dz)
%   dt = courant([dx dy dz]) returns dx/c/sqrt(1/dx^2 + 1/dy^2 + 1/dz^2)
%
% You may want to multiply dt by 0.99 to guarantee stability.

if nargin == 1
    if length(varargin{1}) == 1
        dxyz = [varargin{1}];
    elseif length(varargin{1}) == 2
        dxyz = [varargin{1}(1) varargin{1}(2)];
    elseif length(varargin{1}) == 3
        dxyz = [varargin{1}(1) varargin{1}(2) varargin{1}(3)];
    else
        error('Invalid cell size.');
    end
elseif nargin == 2
    dxyz = [varargin{1} varargin{2}];
elseif nargin == 3
    dxyz = [varargin{1} varargin{2} varargin{3}];
else
    error('Invalid cell size.');
end

dt = 1 / 3e8 / sqrt(sum(1./dxyz.^2));
