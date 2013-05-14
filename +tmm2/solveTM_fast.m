function outStruct = solveTM_fast(boundaryZ, epsr, mur, omega, ky, ...
    sourceHxLeft, sourceHxRight, ...
    Mx, Jy, ...
    posHx, posHxy, posHxz, ...
    posEy, posEyy, posEyz, ...
    posEz, posEzy, posEzz, ...
    forceBoundModes)

import tmm.*;

numLayers = length(epsr);
numBoundaries = numLayers-1;

rowVec = @(A) reshape(A, 1, []);
JM = [rowVec(Jy); rowVec(Mx)];

n = sqrt(epsr.*mur);

ks = sqrt(omega^2*n.^2 - ky^2);
ks(imag(ks) < 0) = -ks(imag(ks) < 0); % decay goes the right way now

% Make matrices to convert [H+, H-] in layer n to [H+, H-] in layer n+1,
% and vice-versa.
% 
% We will use a forward step and a backward step to find all the forward
% and backward wave amplitudes.

% H(n+1) = forwardMatrix{n}*H(n) - sourceMatrix{n}*JM(n)
% H(n) = backwardMatrix{n}*H(n+1) + sourceMatrix{n}*JM(n)

forwardMatrix = cell(length(boundaryZ), 1);
sourceMatrix = forwardMatrix;

for nLayer = 1:numBoundaries
    forwardMatrix{nLayer} = ...
        matrixHE2HH(omega, ks(nLayer+1), boundaryZ(nLayer), epsr(nLayer+1)) * ...
        matrixHH2HE(omega, ks(nLayer), boundaryZ(nLayer), epsr(nLayer));
    
    sourceMatrix{nLayer} = ...
        matrixHE2HH(omega, ks(nLayer+1), boundaryZ(nLayer), epsr(nLayer+1));
end

% The total source matrices are what we multiply source J and M by to get
% the source term for the final, like, whatever it is.

totalJM = [0; 0];

if numBoundaries > 0
    totalSourceMatrix = cell(numBoundaries, 1);
    totalSourceMatrix{numBoundaries} = ...
        matrixHE2HH(omega, ks(numLayers), boundaryZ(numLayers-1), epsr(numLayers));
    for nLayer = numBoundaries-1:-1:1
        totalSourceMatrix{nLayer} = totalSourceMatrix{nLayer+1} * ...
            matrixHH2HE(omega, ks(nLayer+1), boundaryZ(nLayer+1), epsr(nLayer+1)) * ...
            matrixHE2HH(omega, ks(nLayer+1), boundaryZ(nLayer), epsr(nLayer+1));
    end

    for nLayer = 1:numBoundaries
        totalJM = totalJM - totalSourceMatrix{nLayer}*JM(:,nLayer);
    end
end

%% Create the transfer matrix.

% This will the transfer matrix that provides the last layer's E in terms
% of the first layer's E.
transferHH = eye(2);

% E_N = transferLayer{N}*E_first;
transferLayer = cell(length(boundaryZ)+1, 1); % for intermediate layers

transferLayer{1} = transferHH;
for nLayer = 1:length(forwardMatrix)
    transferHH = forwardMatrix{nLayer}*transferHH;
    transferLayer{nLayer+1} = transferHH;
end

%% Create scattering matrix-like system:

A = [1, -transferHH(1,2); 0, -transferHH(2,2)];
B = [0, -transferHH(1,1); 1, -transferHH(2,1)];

% Solve A*H_outward = totalJM - B*H_inward.
% H_outward = [HN+; H1-]
% H_inward = [HN-; H1+]

H_outward = A \ (totalJM - B*[sourceHxRight; sourceHxLeft]);

H0 = [ sourceHxLeft; H_outward(2) ];

%% For mode solutions:

if forceBoundModes
    H0(1) = 0;
end

%% Get the forward and backward E in each layer (use transferLayer)

outStruct = struct;

intervals = [-inf, boundaryZ(:)', inf];

%% Forward/backward coefficients.  Get them all...
H = zeros(2, numLayers);
H(:,1) = transferLayer{1}*H0;
for nLayer = 2:numLayers
    H(:,nLayer) = forwardMatrix{nLayer-1}*H(:,nLayer-1) - ...
        sourceMatrix{nLayer-1} * JM(:, nLayer-1);
end

%% Outputs!
% First allocate them (YAWN)

if ~isempty(posHx) outStruct.Hx = 0*posHx; end
if ~isempty(posHxz) outStruct.Hxz = 0*posHxz; end
if ~isempty(posHxy) outStruct.Hxy = 0*posHxy; end
if ~isempty(posEy) outStruct.Ey = 0*posEy; end
if ~isempty(posEyz) outStruct.Eyz = 0*posEyz; end
if ~isempty(posEyy) outStruct.Eyy = 0*posEyy; end
if ~isempty(posEz) outStruct.Ez = 0*posEz; end
if ~isempty(posEzz) outStruct.Ezz = 0*posEzz; end
if ~isempty(posEzy) outStruct.Ezy = 0*posEzy; end

%ddy = diag([1i*ky, 1i*ky]);
ddy = 1i*ky;
            
for nLayer = 1:numLayers
    
    ddz = diag([1i*ks(nLayer), -1i*ks(nLayer)]);
    
    if ~isempty(posHx)
        indices = find( posHx > intervals(nLayer) & ...
            posHx <= intervals(nLayer+1));
        
        for ii = indices
            z = posHx(ii);
            hh2he = matrixHH2HE(omega, ks(nLayer), z, epsr(nLayer));
            HE = hh2he*H(:,nLayer);
            outStruct.Hx(ii) = HE(1);
        end
    end
    if ~isempty(posHxz)
        indices = find( posHxz > intervals(nLayer) & ...
            posHxz <= intervals(nLayer+1));
        
        for ii = indices
            z = posHxz(ii);
            hh2he = matrixHH2HE(omega, ks(nLayer), z, epsr(nLayer));
            HE = hh2he*ddz*H(:,nLayer);
            outStruct.Hxz(ii) = HE(1);
        end
    end
    if ~isempty(posHxy)
        indices = find( posHxy > intervals(nLayer) & ...
            posHxy <= intervals(nLayer+1));
        
        for ii = indices
            z = posHxy(ii);
            hh2he = matrixHH2HE(omega, ks(nLayer), z, epsr(nLayer));
            HE = hh2he*ddy*H(:,nLayer);
            outStruct.Hxy(ii) = HE(1);
        end
    end
    
    if ~isempty(posEy)
        indices = find(posEy > intervals(nLayer) & ...
         posEy <= intervals(nLayer+1));
    
        for ii = indices
            z = posEy(ii);
            hh2he = matrixHH2HE(omega, ks(nLayer), z, epsr(nLayer));
            HE = hh2he*H(:,nLayer);
            outStruct.Ey(ii) = HE(2);
        end
    end
    if ~isempty(posEyz)
        indices = find(posEyz > intervals(nLayer) & ...
         posEyz <= intervals(nLayer+1));
    
        for ii = indices
            z = posEyz(ii);
            hh2he = matrixHH2HE(omega, ks(nLayer), z, epsr(nLayer));
            HE = hh2he*ddz*H(:,nLayer);
            outStruct.Eyz(ii) = HE(2);
        end
    end
    if ~isempty(posEyy)
        indices = find(posEyy > intervals(nLayer) & ...
         posEyy <= intervals(nLayer+1));
    
        for ii = indices
            z = posEyy(ii);
            hh2he = matrixHH2HE(omega, ks(nLayer), z, epsr(nLayer));
            HE = hh2he*ddy*H(:,nLayer);
            outStruct.Eyy(ii) = HE(2);
        end
    end
    
    if ~isempty(posEz)
        indices = find(posEz > intervals(nLayer) & ...
         posEz <= intervals(nLayer+1));
     
        for ii = indices
            z = posEz(ii);
            hh2he = matrixHH2HE(omega, ks(nLayer), z, epsr(nLayer));
            HE = hh2he*H(:,nLayer);
            outStruct.Ez(ii) = HE(1)*ky/omega/epsr(nLayer);
        end
    end
    if ~isempty(posEzz)
        indices = find(posEzz > intervals(nLayer) & ...
         posEzz <= intervals(nLayer+1));
     
        for ii = indices
            z = posEzz(ii);
            hh2he = matrixHH2HE(omega, ks(nLayer), z, epsr(nLayer));
            HE = hh2he*ddz*H(:,nLayer);
            outStruct.Ezz(ii) = HE(1)*ky/omega/epsr(nLayer);
        end
    end
    if ~isempty(posEzy)
        indices = find(posEzy > intervals(nLayer) & ...
         posEzy <= intervals(nLayer+1));
     
        for ii = indices
            z = posEzy(ii);
            hh2he = matrixHH2HE(omega, ks(nLayer), z, epsr(nLayer));
            HE = hh2he*ddy*H(:,nLayer);
            outStruct.Ezy(ii) = HE(1)*ky/omega/epsr(nLayer);
        end
    end
    
end

outStruct.H = H;