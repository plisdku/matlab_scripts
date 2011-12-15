function [E, H, timesteps] = solve(o, posE, posH, varargin)
% solve  Return E(x,freqs) and H(x,freqs) in the multilayer
%   [E H] = solve(posE, posH) calculates the fields
%       in the multilayer at the positions given in posE and posH.  The source
%       must have been specified with addSource.
%
%   [E H] = solve(posE, posH, cutoffFrequency) will
%       omit from consideration all frequency components in the source that
%       exceed the given cutoff.  Typically the majority of the spectrum is
%       not needed and using a cutoff can accelerate the calculation.  The
%       default cutoff frequency corresponds to a wavelength of roughly
%       three Yee cells.
%
%   Arguments:
%       posE, posH      Positions at which to return E and H fields.  If they
%                       are arrays of dimension [N 3], then the E fields
%                       (H fields) Ex, Ey and Ez will all be returned for
%                       each position.  They may also be cell arrays containing
%                       different arrays of positions for each field
%                       component, with dimensions {[N 3], [M 3], [P 3]}.
%   Return values:
%       E, H            Cell arrays of electromagnetic fields.
%                       E = {Ex, Ey, Ez} where Ex (Ey, Ez) has size
%                       [N 3 numTimesteps], for N samples in posE or posE{1}

cutoffFrequency = 2*pi*3e8/3/norm(o.dxyz); % cutoff wavelength is 3*dx
if nargin > 3
    if ~isempty(varargin{1})
        cutoffFrequency = varargin{1};
        if cutoffFrequency > pi/o.dt
            warning('Cutoff frequency is above Nyquist frequency; this may cause problems');
        end
    end
end

if ~iscell(posE)
    posE = {posE posE posE};
end
if ~iscell(posH)
    posH = {posH posH posH};
end

for xyz = 1:3
    assert(size(posE{xyz},2) == 3);
    assert(size(posH{xyz},2) == 3);
end

dx = o.dxyz(1);
dy = o.dxyz(2);
dz = o.dxyz(3);

timesteps = 1:size(o.source,2);
numOutputSteps = length(timesteps);

omegas = 2*pi*(0:numOutputSteps-1)/(o.dt*numOutputSteps);
freqIndices = find(omegas < cutoffFrequency);

[fieldPhasors, z, zParallel, zNormal] = o.phasors(omegas(freqIndices));

%%

boundsE = [-inf, o.boundariesE, inf];
boundsH = [-inf, o.boundariesH, inf];

indicesE = cell(3, o.numLayers(), 1);
indicesH = indicesE;

for nn = 1:o.numLayers()
    for xyz = 1:3
        indicesE{xyz,nn} = find(posE{xyz}(:,3) > boundsE(nn) & ...
            posE{xyz}(:,3) <= boundsE(nn+1));
        indicesH{xyz,nn} = find(posH{xyz}(:,3) > boundsH(nn) & ...
            posH{xyz}(:,3) <= boundsH(nn+1));
    end
end


%% Fourier representation of fields

E = cell(3,1);
H = cell(3,1);
for xyz = 1:3
    E{xyz} = zeros(length(posE{xyz}), numel(freqIndices));
    H{xyz} = zeros(length(posH{xyz}), numel(freqIndices));
end

localPosE = cell(3,1);
localPosH = cell(3,1);

tic;
for nn = 1:o.numLayers()
    
    for xyz = 1:3
        localPosE{xyz} = posE{xyz}(indicesE{xyz,nn},:);
        localPosE{xyz}(:,3) = localPosE{xyz}(:,3) - o.layerOrigins(nn);
        
        localPosH{xyz} = posH{xyz}(indicesH{xyz,nn},:);
        localPosH{xyz}(:,3) = localPosH{xyz}(:,3) - o.layerOrigins(nn);
    end
    
    for ff = freqIndices % select lowest frequencies
        %fprintf('Working on frequency %i\n', ff);
        for xyz = 1:3
            transverseFactorE = ...
                zParallel(1,ff).^(localPosE{xyz}(:,1)/dx) .*...
                zParallel(2,ff).^(localPosE{xyz}(:,2)/dy);
            forwardFactorE = zNormal(nn,ff).^(localPosE{xyz}(:,3)/dz);
            
            E{xyz}(indicesE{xyz,nn},ff) = transverseFactorE.* (...
                fieldPhasors(o.eForward(xyz,nn),ff).*forwardFactorE +...
                fieldPhasors(o.eBackward(xyz,nn),ff)./forwardFactorE );
            
            transverseFactorH = ...
                zParallel(1,ff).^(localPosH{xyz}(:,1)/dx) .*...
                zParallel(2,ff).^(localPosH{xyz}(:,2)/dy);
            forwardFactorH = zNormal(nn,ff).^(localPosH{xyz}(:,3)/dz);
            
            H{xyz}(indicesH{xyz,nn},ff) = transverseFactorH.* (...
                fieldPhasors(o.hForward(xyz,nn),ff).*forwardFactorH +...
                fieldPhasors(o.hBackward(xyz,nn),ff)./forwardFactorH );
        end
    end
end
evaluateFieldsTime = toc;

fprintf('Time to calculate field phasors at all output points: %2.2f s\n',...
    evaluateFieldsTime);

% NaNs can appear in the analysis sometimes but it doesn't seem to hurt to
% throw them out...
for xyz = 1:3
    E{xyz}(isnan(E{xyz})) = 0;
    H{xyz}(isnan(H{xyz})) = 0;
end







