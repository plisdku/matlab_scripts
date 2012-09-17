function [zNumerator, zDenominator] = newMaterial(name, varargin)
%newMaterial Declare a new permittivity or permeability type
%
% Materials are IIR filters describing permittivity or permeability.  They
% may be expressed as rational functions in z or s, where s = jw and z =
% exp(s*dt).  Allowed arguments:
%
%   Numerator, Denominator    Descending-order coefficients of polynomials in s
%   Let Numerator = [n0 n1 n2].  Its interpretation is a polynomial
%       N(s) = n0*s^2 + n1*s + n2.
%   Denominator, likewise, is a vector of polynomial coefficients.
%   
%   ZNumerator, ZDenominator  Descending-order coefficients of polynomials in z
%   Let ZNumerator = [n0 n1 n2].  Its interpretation is a polynomial
%       N(z) = n0*z^2 + n1*z + n2.
%   Denominator, likewise, is a vector of polynomial coefficients.
%   
% Examples:
% 
% newMaterial('Air', 'Numerator', 1, 'Denominator', 1) creates free space.
%
% newMaterial('Dielectric', 'Numerator', 4, 'Denominator', 1) creates a
% nondispersive dielectric with permittivity (or permeability) 4.
% 
% Drude material: relative permittivity is
%   epsr = epsinf - wp^2/(w^2 + j*w*gamma).
% Convert to s = j*w:
%   epsr = epsinf - wp^2/(-s^2 + s*gamma)
%

sim = t6.TrogdorSimulation.instance();

X.Numerator = [];
X.Denominator = [];
X.ZNumerator = [];
X.ZDenominator = [];

if numel(varargin) == 2 && all(cellfun(@isnumeric, varargin))
    X.ZNumerator = varargin{1};
    X.ZDenominator = varargin{2};
else
    X = parseargs(X, varargin{:});
end

validateInput(X);

% Numerator and Denominator are polynomials in the frequency.  They need to
% get turned into polynomials in z = exp(jwt).  Test whether we have
% frequency or z descriptions, and convert appropriately.
if ~isempty(X.Numerator)
    
    if numel(X.Numerator) > 1
        dt = sim.Dt();
        [zNumerator, zDenominator] = bilinear(X.Numerator, X.Denominator, 1/dt);
    else
        zNumerator = X.Numerator;
        zDenominator = X.Denominator;
    end
    
    % bilinear takes X.Numerator and X.Denominator to be in decreasing
    % powers of s.  Blarg.
else
    zNumerator = X.ZNumerator;
    zDenominator = X.ZDenominator;
end

if ~t6.testStability(zNumerator, zDenominator, sim.Dt, sim.Dxyz, sim.NumCells)
    warning('Material %s appears to be unstable as a permittivity', name);
end

material = struct('name', name, 'numerator', zNumerator, ...
    'denominator', zDenominator);

sim.Materials = {sim.Materials{:}, material};


function validateInput(X)

if ~isempty(X.Numerator) && ~isempty(X.Denominator)
    if ~isempty(X.ZNumerator) || ~isempty(X.ZDenominator)
        error(['Specify material only with Numerator and Denominator or '...
            'ZNumerator and ZDenominator']);
    end

    validateattributes(X.Numerator, {'numeric'}, {'vector'}, mfilename,...
        'Numerator');
    validateattributes(X.Denominator, {'numeric'}, {'vector'}, mfilename,...
        'Numerator');
    
    % Feb 16 2012: the requirement of same order for numerator and
    % denominator is internal to Trogdor, partly, since it simplified
    % programming.  Additionally it makes it ok to just send the output of
    % Matlab's bilinear() function in to Trogdor as material coefficients.
    % bilinear() produces a polynomial in z, but Trogdor takes a polynomial
    % in 1/z.  If the numerator and denominator are the same order, then
    % the coefficients for the z and 1/z rational functions are in fact the
    % same.
    if numel(X.Numerator) ~= numel(X.Denominator)
        error(['Presently we require the numerator and denominator to have '...
            'the same order']);
    end
    
elseif ~isempty(X.ZNumerator) && ~isempty(X.ZDenominator)
    if ~isempty(X.Numerator) || ~isempty(X.Denominator)
        error(['Specify material only with Numerator and Denominator or '...
            'ZNumerator and ZDenominator']);
    end

    validateattributes(X.ZNumerator, {'numeric'}, {'vector'}, mfilename,...
        'Numerator');
    validateattributes(X.ZDenominator, {'numeric'}, {'vector'}, mfilename,...
        'Numerator');
    
    if numel(X.Numerator) ~= numel(X.Denominator)
        error(['Presently we require the numerator and denominator to have '...
            'the same order']);
    end
else
    error(['Expected Numerator and Denominator or ZNumerator and '...
        'ZDenominator arguments'])
end



