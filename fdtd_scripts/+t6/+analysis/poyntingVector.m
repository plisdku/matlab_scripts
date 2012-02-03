function [pv, timeOrFreq, positions] = poyntingVector(fileName, varargin)
% poyntingVector Calculate the Poynting vector (time or frequency domain)
%   from a Trogdor output file.
%
% The Trogdor output file must contain the fields ex, ey, ez, hx, hy and
% hz.  Such output files may be conveniently added to the simulation with
% addPoyntingPlane(), addPoyntingSurface(), and addPoyntingVolume().
%
% [pv, freqs, positions] = poyntingVector(filename) returns the complex
% Poynting vector 0.5*cross(E, conj(H)) from an output file at all
% frequencies used by fft().  The spectra of E and H are obtained using
% spectrum() and include phase correction.
%
% [pv, freqs, positions] = poyntingVector(filename, 'Frequency',
% frequencies) returns the complex Poynting vector at the given
% frequencies only.
%
% [pv, freqs, positions] = poyntingVector(filename, 'SteadyStateFrequency',
% frequency) returns the complex Poynting vector at a single frequency from
% a simulation with one dominant frequency (e.g. a simulation illuminated
% by a single wavelength and run to steady state).  See spectrum() for
% details on the steady state calculation.
%
% [pv, time, positions] = poyntingVector(filename, 'Time') returns the
% Poynting vector cross(E, H) from an output file.  The time corresponding
% to each timestep and the position of each sample in space are returned as
% well.  Due to FDTD time interleaving, the E and H fields are measured a
% half timestep apart; this is not corrected for in any way. Note as well
% that the time domain Poynting vector cross(E,H) does not have the same
% prefactor as used for the frequency domain Poynting vector, 0.5*cross(E,
% conj(H)).
%
% To obtain the flux of the Poynting vector through a plane or closed box,
% consider using poyntingVectorFlux().

if nargin == 2 && strcmpi(varargin{1}, 'Time')
    varargin{2} = 1;
end

X.SteadyStateFrequency = [];
X.Frequency = [];
X.Time = [];
X = parseargs(X, varargin{:});


file = t6.OutputFile(fileName);
if numel(file.Fields) ~= 6
    error('File must store ex, ey, ez, hx, hy and hz');
end

neededFields = {'ex', 'ey', 'ez', 'hx', 'hy', 'hz'};

for nn = 1:6
    if ~strcmp(file.Fields{nn}.Name, neededFields{nn})
        error('File must store ex, ey, ez, hx, hy and hz');
    end
end

if X.Time
    % time-domain PV.
    % timeOrFreq = timesteps of E, and will be understood to be a little
    % wrong
    
    [data, positions] = t6.readOutputFile(fileName);
    
    if iscell(data)
        pv = cell(size(file.Regions));
        for rr = 1:length(file.Regions)
            pv{rr} = cross(data{rr}(:,:,:,1:3,:), data{rr}(:,:,:,4:6,:), 4);
        end
    else
        pv = cross(data(:,:,:,1:3,:), data(:,:,:,4:6,:), 4);
    end
    
    timeOrFreq = file.times();
    if iscell(timeOrFreq)
        timeOrFreq = timeOrFreq{1};
    end
    
else
    % Frequency domain!
    
    [data, timeOrFreq] = t6.analysis.spectrum(fileName, varargin{:});
    
    if iscell(data)
        pv = cell(size(file.Regions));
        for rr = 1:length(file.Regions)
            pv{rr} = 0.5*cross(data{rr}(:,:,:,1:3,:), ...
                conj(data{rr}(:,:,:,4:6,:)), 4);
        end
    else
        pv = 0.5*cross(data(:,:,:,1:3,:), conj(data(:,:,:,4:6,:)), 4);
    end
end

if nargout > 2
    [posX, posY, posZ] = file.positions();
    
    aField = 1;
    if iscell(data)
        positions = cell(numel(file.Regions),3);
        
        for rr = 1:numel(file.Regions)
            positions{rr,1} = posX{rr,aField};
            positions{rr,2} = posY{rr,aField};
            positions{rr,3} = posZ{rr,aField};
        end
        
    else
        positions{1} = posX{aField};
        positions{2} = posY{aField};
        positions{3} = posZ{aField};
    end
end

