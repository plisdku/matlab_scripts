function [E, H] = phasors(direction, omega, dielectricEpsilon, modeIndex, ...
    posE, posH)
% phasors   Obtain frequency-dependent electromagnetic fields for waveguide mode
% Usage: phasors(direction, frequency, dielectricEpsilon)
%        phasors(direction, frequency, dielectricEpsilon, modeIndex)
%

import t5.dev.multilayerMode.*

tfsfreq_0;

%posE = {afp.posEx, afp.posEy, afp.posEz};
%posH = {afp.posHx, afp.posHy, afp.posHz};

lambda = 2*pi*3e8./omega;
k0 = 2*pi/lambda;
kMin = 1.01*k0;
kMax = 3*k0;

epsr = dielectricEpsilon;
mur = ones(size(epsr));

%% Determine symmetry direction

gridSampleSize = afp.sampleHalfCells(4:6) - afp.sampleHalfCells(1:3);
if sum(gridSampleSize > 0) < 1
    error('No dielectric structure is apparent in this grid');
elseif sum(gridSampleSize > 0) > 1
    error('This grid has less than two translationally-symmetric directions')
end

multilayerXYZ = find(gridSampleSize);
layerNormal = [1 1 1] .* (gridSampleSize > 0);
%directionXYZ = find(direction);

if dot(layerNormal, direction) ~= 0
    error('Propagation direction must be tangent to multilayer structure');
end

% Determination of modes just requires the modal character (TE or TM) and
% the frequency and the structure.  We can save all the coordinate
% transformation issues for later.

xyz = multilayerXYZ;

%% Prepare for coordinate transformations
% I have the transformation
%
% x_AFP = A * x_TMM
%
% and for vector fields like the electric field E,
%
% E_AFP = A * E_TMM.
%
% The TMM code calculates Ex and Ey only.  So, my transformation matrices
% for E and H are actually a little smaller.

tmmX = reshape(direction/norm(direction), [], 1);
tmmY = reshape(layerNormal/norm(layerNormal), [], 1);
tmmZ = cross(tmmX, tmmY);

A = [tmmX tmmY tmmZ];

AE = A(:,1:2);     % E_afp = AE * E_tmm
AH = A(:,3);       % H_afp = AH * H_tmm

%% Find the positions of material boundaries

yy = linspace(afp.sampleBounds(xyz), afp.sampleBounds(xyz+3), numel(afp.grid));
jumpIndex = find(diff(afp.grid)); % the bottom of each jump.
jumpPos = spline(1:numel(yy), yy, jumpIndex+0.5); % Position of boundaries in m

waveguideBoundaries = jumpPos(2:end-1);

%% Now use modesolver to obtain waveguide mode

[kParallel, t22, kvals] = tmm.modesTM(waveguideBoundaries, epsr, mur, omega, ...
    kMin, kMax);
fprintf('%i possible modes detected\n', numel(kParallel));

%% Get the fields out

kx = kParallel(modeIndex);

yVals = unique([posE{1}(:,xyz); posE{2}(:,xyz); posE{3}(:,xyz); ...
    posH{1}(:,xyz); posH{2}(:,xyz); posH{3}(:,xyz)]);

%[hz ex ey] = tmm.solveTM(waveguideBoundaries, epsr, mur, omega, kx, ...
%    {afp.posHz(:,xyz), afp.posEx(:,xyz), afp.posEy(:,xyz)}, true);
[hz ex ey] = tmm.solveTM(waveguideBoundaries, epsr, mur, omega, kx, yVals, ...
    true);

%% Now we need to transform the fields the right way.
% There's a coordinate transformation and a vector transformation.

E = {[], [], []};
H = {[], [], []};

% Do each field separately.
for exyz = 1:3
    
    % 1. Transform the coordinates from AFP to TMM
    posTMM = (A \ posE{exyz}')';
    
    % 2. Get the fields in TMM-land
    
    tmmE = [subSafeSpline(yVals, ex, posTMM(:,2)') .* exp(1i*kx*posTMM(:,1)'); ...
        subSafeSpline(yVals, ey, posTMM(:,2)') .* exp(1i*kx*posTMM(:,1)')];
    
    % 3. Get the fields for AFP
    E{exyz} = reshape(AE(exyz,:) * tmmE, [], 1);
end

for hxyz = 1:3
    
    % 1. Transform the coordinates from AFP to TMM
    posTMM = (A \ posH{exyz}')';
    
    % 2. Get the fields in TMM-land
    tmmH = subSafeSpline(yVals, hz, posTMM(:,2)') .* exp(1i*kx*posTMM(:,1)');
    
    % 3. Get the fields for AFP
    H{hxyz} = reshape(AH(hxyz,:) * tmmH, [], 1);
end

function yy = subSafeSpline(x, y, xx)
% The reason for this function to exist is that Matlab's spline() function
% will gag if there is only one data point.  Likewise interp1() will bomb
% if there is only one data point, EVEN IF I have extrapolation turned on,
% EVEN IF the method I selected is "nearest neighbor".  It's just hardcoded
% that way, which means it's not robust to cases like spline(0, 1, 0), a
% sort of thing which could very easily happen in my case.  I just want a
% table lookup really; why is this so hard?

if numel(y) > 1
    yy = spline(x,y,xx);
else
    yy = repmat(y, size(xx));
end


