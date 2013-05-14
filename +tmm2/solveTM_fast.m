function outStruct = solveTM_fast(boundaryX, epsr, mur, omega, ky, ...
    sourceHxLeft, sourceHxRight, ...
    Mz, Jy, ...
    posHz, posHzx, posHzy, ...
    posEx, posExx, posExy, ...
    posEy, posEyx, posEyy, ...
    forceBoundModes)

import tmm2.*;

numLayers = length(epsr);
numBoundaries = numLayers-1;

rowVec = @(A) reshape(A, 1, []);
JM = [rowVec(Jy); rowVec(Mz)];

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

forwardMatrix = cell(length(boundaryX), 1);
sourceMatrix = forwardMatrix;

for nLayer = 1:numBoundaries
    forwardMatrix{nLayer} = ...
        matrixHE2HH(omega, ks(nLayer+1), boundaryX(nLayer), epsr(nLayer+1)) * ...
        matrixHH2HE(omega, ks(nLayer), boundaryX(nLayer), epsr(nLayer));
    
    sourceMatrix{nLayer} = ...
        matrixHE2HH(omega, ks(nLayer+1), boundaryX(nLayer), epsr(nLayer+1));
end

% The total source matrices are what we multiply source J and M by to get
% the source term for the final, like, whatever it is.

totalJM = [0; 0];

if numBoundaries > 0
    totalSourceMatrix = cell(numBoundaries, 1);
    totalSourceMatrix{numBoundaries} = ...
        matrixHE2HH(omega, ks(numLayers), boundaryX(numLayers-1), epsr(numLayers));
    for nLayer = numBoundaries-1:-1:1
        totalSourceMatrix{nLayer} = totalSourceMatrix{nLayer+1} * ...
            matrixHH2HE(omega, ks(nLayer+1), boundaryX(nLayer+1), epsr(nLayer+1)) * ...
            matrixHE2HH(omega, ks(nLayer+1), boundaryX(nLayer), epsr(nLayer+1));
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
transferLayer = cell(length(boundaryX)+1, 1); % for intermediate layers

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

intervals = [-inf, boundaryX(:)', inf];

%% Forward/backward coefficients.  Get them all...
H = zeros(2, numLayers);
H(:,1) = transferLayer{1}*H0;
for nLayer = 2:numLayers
    H(:,nLayer) = forwardMatrix{nLayer-1}*H(:,nLayer-1) - ...
        sourceMatrix{nLayer-1} * JM(:, nLayer-1);
end

%% Outputs!
% First allocate them (YAWN)

if ~isempty(posHz) outStruct.Hz = 0*posHz; end
if ~isempty(posHzy) outStruct.Hzy = 0*posHzy; end
if ~isempty(posHzx) outStruct.Hzx = 0*posHzx; end
if ~isempty(posEx) outStruct.Ex = 0*posEx; end
if ~isempty(posExy) outStruct.Exy = 0*posExy; end
if ~isempty(posExx) outStruct.Exx = 0*posExx; end
if ~isempty(posEy) outStruct.Ey = 0*posEy; end
if ~isempty(posEyy) outStruct.Eyy = 0*posEyy; end
if ~isempty(posEyx) outStruct.Eyx = 0*posEyx; end

%ddy = diag([1i*ky, 1i*ky]);
ddy = 1i*ky;
            
for nLayer = 1:numLayers
    
    ddx = diag([1i*ks(nLayer), -1i*ks(nLayer)]);
    
    if ~isempty(posHz)
        indices = find( posHz > intervals(nLayer) & ...
            posHz <= intervals(nLayer+1));
        
        for ii = indices
            x = posHz(ii);
            hh2he = matrixHH2HE(omega, ks(nLayer), x, epsr(nLayer));
            HE = hh2he*H(:,nLayer);
            outStruct.Hz(ii) = HE(1);
        end
    end
    if ~isempty(posHzy)
        indices = find( posHzy > intervals(nLayer) & ...
            posHzy <= intervals(nLayer+1));
        
        for ii = indices
            x = posHzy(ii);
            hh2he = matrixHH2HE(omega, ks(nLayer), x, epsr(nLayer));
            HE = hh2he*ddy*H(:,nLayer);
            outStruct.Hzy(ii) = HE(1);
        end
    end
    if ~isempty(posHzx)
        indices = find( posHzx > intervals(nLayer) & ...
            posHzx <= intervals(nLayer+1));
        
        for ii = indices
            x = posHzx(ii);
            hh2he = matrixHH2HE(omega, ks(nLayer), x, epsr(nLayer));
            HE = hh2he*ddx*H(:,nLayer);
            outStruct.Hzx(ii) = HE(1);
        end
    end
    
    if ~isempty(posEx)
        indices = find(posEx > intervals(nLayer) & ...
         posEx <= intervals(nLayer+1));
    
        for ii = indices
            x = posEx(ii);
            hh2he = matrixHH2HE(omega, ks(nLayer), x, epsr(nLayer));
            HE = hh2he*H(:,nLayer);
            outStruct.Ex(ii) = -HE(1)*ky/omega/epsr(nLayer);
        end
    end
    if ~isempty(posExy)
        indices = find(posExy > intervals(nLayer) & ...
         posExy <= intervals(nLayer+1));
    
        for ii = indices
            x = posExy(ii);
            hh2he = matrixHH2HE(omega, ks(nLayer), x, epsr(nLayer));
            HE = hh2he*ddy*H(:,nLayer);
            outStruct.Exy(ii) = -HE(1)*ky/omega/epsr(nLayer);
        end
    end
    if ~isempty(posExx)
        indices = find(posExx > intervals(nLayer) & ...
         posExx <= intervals(nLayer+1));
    
        for ii = indices
            x = posExx(ii);
            hh2he = matrixHH2HE(omega, ks(nLayer), x, epsr(nLayer));
            HE = hh2he*ddx*H(:,nLayer);
            outStruct.Exx(ii) = -HE(1)*ky/omega/epsr(nLayer);
        end
    end
    
    if ~isempty(posEy)
        indices = find(posEy > intervals(nLayer) & ...
         posEy <= intervals(nLayer+1));
     
        for ii = indices
            x = posEy(ii);
            hh2he = matrixHH2HE(omega, ks(nLayer), x, epsr(nLayer));
            HE = hh2he*H(:,nLayer);
            outStruct.Ey(ii) = HE(2);
        end
    end
    if ~isempty(posEyy)
        indices = find(posEyy > intervals(nLayer) & ...
         posEyy <= intervals(nLayer+1));
     
        for ii = indices
            x = posEyy(ii);
            hh2he = matrixHH2HE(omega, ks(nLayer), x, epsr(nLayer));
            HE = hh2he*ddy*H(:,nLayer);
            outStruct.Eyy(ii) = HE(2);
        end
    end
    if ~isempty(posEyx)
        indices = find(posEyx > intervals(nLayer) & ...
         posEyx <= intervals(nLayer+1));
     
        for ii = indices
            x = posEyx(ii);
            hh2he = matrixHH2HE(omega, ks(nLayer), x, epsr(nLayer));
            HE = hh2he*ddx*H(:,nLayer);
            outStruct.Eyx(ii) = HE(2);
        end
    end
    
end

outStruct.H = H;