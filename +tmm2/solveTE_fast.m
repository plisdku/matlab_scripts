function outStruct = solveTE_fast(boundaryX, epsr, mur, omega, ky, ...
    sourceEzLeft, sourceEzRight, ...
    Jz, Mx, ...
    posEz, posEzx, posEzy, ...
    posHx, posHxx, posHxy, ...
    posHy, posHyx, posHyy, ...
    forceBoundModes)

import tmm2.*;
outStruct = struct;

numLayers = length(epsr);
numBoundaries = numLayers-1;

rowVec = @(A) reshape(A, 1, []);
MJ = [rowVec(Mx); rowVec(Jz)];

n = sqrt(epsr.*mur);

ks = sqrt(omega^2*n.^2 - ky^2);
ks(imag(ks) < 0) = -ks(imag(ks) < 0); % decay goes the right way now
outStruct.kx = ks;

% Make matrices to convert [E+, E-] in layer n to [E+, E-] in layer n+1,
% and vice-versa.
% 
% We will use a forward step and a backward step to find all the forward
% and backward wave amplitudes.

% E(n+1) = forwardMatrix{n}*E(n) + sourceMatrix{n}*MJ(n)
% E(n) = backwardMatrix{n}*E(n+1) - sourceMatrix{n}*MJ(n)

forwardMatrix = cell(length(boundaryX), 1);
sourceMatrix = forwardMatrix;

for nLayer = 1:numBoundaries
    forwardMatrix{nLayer} = ...
        matrixEH2EE(omega, ks(nLayer+1), boundaryX(nLayer), mur(nLayer+1)) * ...
        matrixEE2EH(omega, ks(nLayer), boundaryX(nLayer), mur(nLayer));
    
    sourceMatrix{nLayer} = ...
        matrixEH2EE(omega, ks(nLayer+1), boundaryX(nLayer), mur(nLayer+1));
end

% The total source matrices are what we multiply source J and M by to get
% the source term for the final, like, whatever it is.

totalMJ = [0; 0];

if numBoundaries > 0
    totalSourceMatrix = cell(numBoundaries, 1);
    totalSourceMatrix{numBoundaries} = ...
        matrixEH2EE(omega, ks(numLayers), boundaryX(numLayers-1), mur(numLayers));
    for nLayer = numBoundaries-1:-1:1
        totalSourceMatrix{nLayer} = totalSourceMatrix{nLayer+1} * ...
            matrixEE2EH(omega, ks(nLayer+1), boundaryX(nLayer+1), mur(nLayer+1)) * ...
            matrixEH2EE(omega, ks(nLayer+1), boundaryX(nLayer), mur(nLayer+1));
    end

    for nLayer = 1:numBoundaries
        totalMJ = totalMJ + totalSourceMatrix{nLayer}*MJ(:,nLayer);
    end
end

%% Create the transfer matrix.

% This will the transfer matrix that provides the last layer's E in terms
% of the first layer's E.
transferEE = eye(2);

% E_N = transferLayer{N}*E_first;
transferLayer = cell(length(boundaryX)+1, 1); % for intermediate layers

transferLayer{1} = transferEE;
for nLayer = 1:length(forwardMatrix)
    transferEE = forwardMatrix{nLayer}*transferEE;
    transferLayer{nLayer+1} = transferEE;
end

%% Create scattering matrix-like system:

A = [1, -transferEE(1,2); 0, -transferEE(2,2)];
B = [0, -transferEE(1,1); 1, -transferEE(2,1)];

% Solve A*E_outward = totalMJ - B*E_inward.
% E_outward = [EN+; E1-]
% E_inward = [EN-; E1+]

E_outward = A \ (totalMJ - B*[sourceEzRight; sourceEzLeft]);

E0 = [ sourceEzLeft; E_outward(2) ];

%% For mode solutions:

if forceBoundModes
    E0(1) = 0;
end

%% Get the forward and backward E in each layer (use transferLayer)

intervals = [-inf, boundaryX(:)', inf];

%% Forward/backward coefficients.  Get them all...
E = zeros(2, numLayers);
E(:,1) = transferLayer{1}*E0;
for nLayer = 2:numLayers
    E(:,nLayer) = forwardMatrix{nLayer-1}*E(:,nLayer-1) + ...
        sourceMatrix{nLayer-1} * MJ(:, nLayer-1);
end

%% Outputs!
% First allocate them (YAWN)

if ~isempty(posEz) outStruct.Ez = 0*posEz; end
if ~isempty(posEzy) outStruct.Ezy = 0*posEzy; end
if ~isempty(posEzx) outStruct.Ezx = 0*posEzx; end
if ~isempty(posHx) outStruct.Hx = 0*posHx; end
if ~isempty(posHxy) outStruct.Hxy = 0*posHxy; end
if ~isempty(posHxx) outStruct.Hxx = 0*posHxx; end
if ~isempty(posHy) outStruct.Hy = 0*posHy; end
if ~isempty(posHyy) outStruct.Hyy = 0*posHyy; end
if ~isempty(posHyx) outStruct.Hyx = 0*posHyx; end

%ddy = diag([1i*ky, 1i*ky]);
ddy = 1i*ky;

for nLayer = 1:numLayers
    
    ddx = diag([1i*ks(nLayer), -1i*ks(nLayer)]);
    
    if ~isempty(posEz)
        indices = find( posEz > intervals(nLayer) & ...
            posEz <= intervals(nLayer+1));
        
        for ii = indices
            x = posEz(ii);
            ee2eh = matrixEE2EH(omega, ks(nLayer), x, mur(nLayer));
            EH = ee2eh*E(:,nLayer);
            outStruct.Ez(ii) = EH(1);
        end
    end
    if ~isempty(posEzy)
        indices = find( posEzy > intervals(nLayer) & ...
            posEzy <= intervals(nLayer+1));
        
        for ii = indices
            x = posEzy(ii);
            ee2eh = matrixEE2EH(omega, ks(nLayer), x, mur(nLayer));
            EH = ee2eh*ddy*E(:,nLayer);
            outStruct.Ezy(ii) = EH(1);
        end
    end
    if ~isempty(posEzx)
        indices = find( posEzx > intervals(nLayer) & ...
            posEzx <= intervals(nLayer+1));
        
        for ii = indices
            x = posEzx(ii);
            ee2eh = matrixEE2EH(omega, ks(nLayer), x, mur(nLayer));
            EH = ee2eh*ddx*E(:,nLayer);
            outStruct.Ezx(ii) = EH(1);
        end
    end
    
    if ~isempty(posHx)
        indices = find(posHx > intervals(nLayer) & ...
         posHx <= intervals(nLayer+1));
    
        for ii = indices
            x = posHx(ii);
            ee2eh = matrixEE2EH(omega, ks(nLayer), x, mur(nLayer));
            EH = ee2eh*E(:,nLayer);
            outStruct.Hx(ii) = EH(1)*ky/omega/mur(nLayer);
        end
    end
    if ~isempty(posHxy)
        indices = find(posHxy > intervals(nLayer) & ...
         posHxy <= intervals(nLayer+1));
    
        for ii = indices
            x = posHxy(ii);
            ee2eh = matrixEE2EH(omega, ks(nLayer), x, mur(nLayer));
            EH = ee2eh*ddy*E(:,nLayer);
            outStruct.Hxy(ii) = EH(1)*ky/omega/mur(nLayer);
        end
    end
    if ~isempty(posHxx)
        indices = find(posHxx > intervals(nLayer) & ...
         posHxx <= intervals(nLayer+1));
    
        for ii = indices
            x = posHxx(ii);
            ee2eh = matrixEE2EH(omega, ks(nLayer), x, mur(nLayer));
            EH = ee2eh*ddx*E(:,nLayer);
            outStruct.Hxx(ii) = EH(1)*ky/omega/mur(nLayer);
        end
    end
    
    if ~isempty(posHy)
        indices = find(posHy > intervals(nLayer) & ...
         posHy <= intervals(nLayer+1));
     
        for ii = indices
            x = posHy(ii);
            ee2eh = matrixEE2EH(omega, ks(nLayer), x, mur(nLayer));
            EH = ee2eh*E(:,nLayer);
            outStruct.Hy(ii) = EH(2);
        end
    end
    if ~isempty(posHyy)
        indices = find(posHyy > intervals(nLayer) & ...
         posHyy <= intervals(nLayer+1));
     
        for ii = indices
            x = posHyy(ii);
            ee2eh = matrixEE2EH(omega, ks(nLayer), x, mur(nLayer));
            EH = ee2eh*ddy*E(:,nLayer);
            outStruct.Hyy(ii) = EH(2);
        end
    end
    if ~isempty(posHyx)
        indices = find(posHyx > intervals(nLayer) & ...
         posHyx <= intervals(nLayer+1));
     
        for ii = indices
            x = posHyx(ii);
            ee2eh = matrixEE2EH(omega, ks(nLayer), x, mur(nLayer));
            EH = ee2eh*ddx*E(:,nLayer);
            outStruct.Hyx(ii) = EH(2);
        end
    end
    
end

outStruct.E = E;
