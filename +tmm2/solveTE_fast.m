function outStruct = solveTE_fast(boundaryZ, epsr, mur, omega, ky, ...
    sourceExLeft, sourceExRight, ...
    Jx, My, ...
    posEx, posExy, posExz, ...
    posHy, posHyy, posHyz, ...
    posHz, posHzy, posHzz, ...
    forceBoundModes)

import tmm.*;

numLayers = length(epsr);
numBoundaries = numLayers-1;

rowVec = @(A) reshape(A, 1, []);
MJ = [rowVec(My); rowVec(Jx)];

n = sqrt(epsr.*mur);

ks = sqrt(omega^2*n.^2 - ky^2);
ks(imag(ks) < 0) = -ks(imag(ks) < 0); % decay goes the right way now

% Make matrices to convert [E+, E-] in layer n to [E+, E-] in layer n+1,
% and vice-versa.
% 
% We will use a forward step and a backward step to find all the forward
% and backward wave amplitudes.

% E(n+1) = forwardMatrix{n}*E(n) + sourceMatrix{n}*MJ(n)
% E(n) = backwardMatrix{n}*E(n+1) - sourceMatrix{n}*MJ(n)

forwardMatrix = cell(length(boundaryZ), 1);
sourceMatrix = forwardMatrix;

for nLayer = 1:numBoundaries
    forwardMatrix{nLayer} = ...
        matrixEH2EE(omega, ks(nLayer+1), boundaryZ(nLayer), mur(nLayer+1)) * ...
        matrixEE2EH(omega, ks(nLayer), boundaryZ(nLayer), mur(nLayer));
    
    sourceMatrix{nLayer} = ...
        matrixEH2EE(omega, ks(nLayer+1), boundaryZ(nLayer), mur(nLayer+1));
end

% The total source matrices are what we multiply source J and M by to get
% the source term for the final, like, whatever it is.

totalMJ = [0; 0];

if numBoundaries > 0
    totalSourceMatrix = cell(numBoundaries, 1);
    totalSourceMatrix{numBoundaries} = ...
        matrixEH2EE(omega, ks(numLayers), boundaryZ(numLayers-1), mur(numLayers));
    for nLayer = numBoundaries-1:-1:1
        totalSourceMatrix{nLayer} = totalSourceMatrix{nLayer+1} * ...
            matrixEE2EH(omega, ks(nLayer+1), boundaryZ(nLayer+1), mur(nLayer+1)) * ...
            matrixEH2EE(omega, ks(nLayer+1), boundaryZ(nLayer), mur(nLayer+1));
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
transferLayer = cell(length(boundaryZ)+1, 1); % for intermediate layers

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

E_outward = A \ (totalMJ - B*[sourceExRight; sourceExLeft]);

E0 = [ sourceExLeft; E_outward(2) ];

%% For mode solutions:

if forceBoundModes
    E0(1) = 0;
end

%% Get the forward and backward E in each layer (use transferLayer)

outStruct = struct;

intervals = [-inf, boundaryZ(:)', inf];

%% Forward/backward coefficients.  Get them all...
E = zeros(2, numLayers);
E(:,1) = transferLayer{1}*E0;
for nLayer = 2:numLayers
    E(:,nLayer) = forwardMatrix{nLayer-1}*E(:,nLayer-1) + ...
        sourceMatrix{nLayer-1} * MJ(:, nLayer-1);
end

%% Outputs!
% First allocate them (YAWN)

if ~isempty(posEx) outStruct.Ex = 0*posEx; end
if ~isempty(posExz) outStruct.Exz = 0*posExz; end
if ~isempty(posExy) outStruct.Exy = 0*posExy; end
if ~isempty(posHy) outStruct.Hy = 0*posHy; end
if ~isempty(posHyz) outStruct.Hyz = 0*posHyz; end
if ~isempty(posHyy) outStruct.Hyy = 0*posHyy; end
if ~isempty(posHz) outStruct.Hz = 0*posHz; end
if ~isempty(posHzz) outStruct.Hzz = 0*posHzz; end
if ~isempty(posHzy) outStruct.Hzy = 0*posHzy; end

%ddy = diag([1i*ky, 1i*ky]);
ddy = 1i*ky;
            
for nLayer = 1:numLayers
    
    ddz = diag([1i*ks(nLayer), -1i*ks(nLayer)]);
    
    if ~isempty(posEx)
        indices = find( posEx > intervals(nLayer) & ...
            posEx <= intervals(nLayer+1));
        
        for ii = indices
            z = posEx(ii);
            ee2eh = matrixEE2EH(omega, ks(nLayer), z, mur(nLayer));
            EH = ee2eh*E(:,nLayer);
            outStruct.Ex(ii) = EH(1);
        end
    end
    if ~isempty(posExz)
        indices = find( posExz > intervals(nLayer) & ...
            posExz <= intervals(nLayer+1));
        
        for ii = indices
            z = posExz(ii);
            ee2eh = matrixEE2EH(omega, ks(nLayer), z, mur(nLayer));
            EH = ee2eh*ddz*E(:,nLayer);
            outStruct.Exz(ii) = EH(1);
        end
    end
    if ~isempty(posExy)
        indices = find( posExy > intervals(nLayer) & ...
            posExy <= intervals(nLayer+1));
        
        for ii = indices
            z = posExy(ii);
            ee2eh = matrixEE2EH(omega, ks(nLayer), z, mur(nLayer));
            EH = ee2eh*ddy*E(:,nLayer);
            outStruct.Exy(ii) = EH(1);
        end
    end
    
    if ~isempty(posHy)
        indices = find(posHy > intervals(nLayer) & ...
         posHy <= intervals(nLayer+1));
    
        for ii = indices
            z = posHy(ii);
            ee2eh = matrixEE2EH(omega, ks(nLayer), z, mur(nLayer));
            EH = ee2eh*E(:,nLayer);
            outStruct.Hy(ii) = EH(2);
        end
    end
    if ~isempty(posHyz)
        indices = find(posHyz > intervals(nLayer) & ...
         posHyz <= intervals(nLayer+1));
    
        for ii = indices
            z = posHyz(ii);
            ee2eh = matrixEE2EH(omega, ks(nLayer), z, mur(nLayer));
            EH = ee2eh*ddz*E(:,nLayer);
            outStruct.Hyz(ii) = EH(2);
        end
    end
    if ~isempty(posHyy)
        indices = find(posHyy > intervals(nLayer) & ...
         posHyy <= intervals(nLayer+1));
    
        for ii = indices
            z = posHyy(ii);
            ee2eh = matrixEE2EH(omega, ks(nLayer), z, mur(nLayer));
            EH = ee2eh*ddy*E(:,nLayer);
            outStruct.Hyy(ii) = EH(2);
        end
    end
    
    if ~isempty(posHz)
        indices = find(posHz > intervals(nLayer) & ...
         posHz <= intervals(nLayer+1));
     
        for ii = indices
            z = posHz(ii);
            ee2eh = matrixEE2EH(omega, ks(nLayer), z, mur(nLayer));
            EH = ee2eh*E(:,nLayer);
            outStruct.Hz(ii) = -EH(1)*ky/omega/mur(nLayer);
        end
    end
    if ~isempty(posHzz)
        indices = find(posHzz > intervals(nLayer) & ...
         posHzz <= intervals(nLayer+1));
     
        for ii = indices
            z = posHzz(ii);
            ee2eh = matrixEE2EH(omega, ks(nLayer), z, mur(nLayer));
            EH = ee2eh*ddz*E(:,nLayer);
            outStruct.Hzz(ii) = -EH(1)*ky/omega/mur(nLayer);
        end
    end
    if ~isempty(posHzy)
        indices = find(posHzy > intervals(nLayer) & ...
         posHzy <= intervals(nLayer+1));
     
        for ii = indices
            z = posHzy(ii);
            ee2eh = matrixEE2EH(omega, ks(nLayer), z, mur(nLayer));
            EH = ee2eh*ddy*E(:,nLayer);
            outStruct.Hzy(ii) = -EH(1)*ky/omega/mur(nLayer);
        end
    end
    
end

outStruct.E = E;
