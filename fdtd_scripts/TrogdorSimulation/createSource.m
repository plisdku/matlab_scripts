% createSource  Initialize a Source object for a Trogdor simulation
%   createSource(field, polarization, extent, data) sets up source data
%   to be written into a Trogdor grid.  The given data is written to a file and Trogdor loads it
%   through an Input object at runtime.  The data should be big enough to fill
%   extent for all required timesteps, so size(data) = [Nx Ny Nz Nt];  However
%   for low-dimensional sources it is okay to leave out singular dimensions, as
%   in the example.
%
%   Example: make a sinusoidal source at the origin
%   
%   numT = 1000;
%   
%   s = createSource('electric', [0 0 1], [0 0 0 0 0 0], sin((1:numT)/20));
%   g = createGrid(..., s);
%   
%   createSource(field, polarization, extent, formula) sets up a source
%   to be written into a Trogdor grid.  The formula, as a string, is passed to
%   Trogdor for parsing.
%
%   The formula may use ordinary functions like sin(), cos(), abs(), exp(),
%   etc. and may make use of the variables n (timestep number) and t (t = n*dt).
%   Currently only real numbers are supported.
%
%   Example: make a sinusoidal Ez source at the origin
%   
%   s = createSource('electric', [0 0 1], [0 0 0 0 0 0], 'sin(n/20)');
%   g = createGrid(..., s);
%
%   createSource(field, polarization, extent, data, 'soft') makes a "soft
%   source" that adds E field to the given region instead of fixing the field.
%   This is equivalent to a current source provided that you know the dispersive
%   properties of the material at the source.  If you place the soft source in
%   a region of permittivity epsr*eps0, then
%
%       E_soft = dE/dt = - J_soft / epsr / eps0
%
%   is the relation between the amplitude E_soft of the soft source and the
%   current of the equivalent current source.  That is, your source field is now
%	an additional component to the time derivative of the field.
%
%   See also: createPlaneWaveSource, createOneFieldInput, createLink
%
%   version 4.5
%   July 29, 2008

function source = createSource(varargin);
% old style: (fieldName, dimensions, data)
% new styles: (fieldName, polarization, dimensions, data)
%             (fieldName, polarization, dimensions, formula)

source.type = 'Source';
source.isSoft = 0;

if (nargin == 3)   % Handle LEGACY SOURCES
    warning('Script uses deprecated arguments for createSource.  Proceeding.');
    fieldName = varargin{1};
    dimensions = varargin{2};
    data = varargin{3};
    
    if (fieldName(1) == 'e')
        source.fieldName = 'electric';
    else
        source.fieldName = 'magnetic';
    end
    
    if (fieldName(2) == 'x')
        source.polarization = [1 0 0];
    elseif (fieldName(2) == 'y')
        source.polarization = [0 1 0];
    else
        source.polarization = [0 0 1];
    end
    
    source.dimensions = dimensions;
    source.data = data;
else
    fieldName = varargin{1};
    polarization = varargin{2};
    dimensions = varargin{3};
    dataOrFormula = varargin{4};
    
    source.fieldName = fieldName;
    source.polarization = polarization;
    source.dimensions = dimensions;
    
    if ischar(dataOrFormula)
        source.formula = dataOrFormula;
    else
        source.data = dataOrFormula;
    end
    
    if (nargin >= 5)
        if (strcmp(varargin{5}, 'soft'))
            source.isSoft = 1;
        end
    end
end



