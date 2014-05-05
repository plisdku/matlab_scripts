function [dfdp, dfdv, dvdp] = calculateAdjointSensitivity(node, parameters,...
    outDir)

import t6.*
import t6.adjoint.*

fwdFileNames = {'boundary_dexx', 'boundary_deyy', 'boundary_dezz'};
adjFileNames = {'boundary_adjoint_exx', 'boundary_adjoint_eyy', ...
    'boundary_adjoint_ezz'};
adjRevFileNames = {'boundary_adjoint_exx.rev', 'boundary_adjoint_eyy.rev',...
    'boundary_adjoint_ezz.rev'};

if ~exist('outDir')
    outDir = 'output';
end
out = @(str) sprintf('%s/%s', outDir, str);

[coeffs, verts] = readDeps(out('Depsilon'));

bigJacobian = node.jacobian(parameters);
vertices = node.vertices(parameters);
vertJacobian = sparse(size(vertices, 1), 1);

%%

for fieldXYZ = 1:3

fprintf('Reversing file, field %i of 3\n', fieldXYZ);

fwdDE = OutputFile(out(fwdFileNames{fieldXYZ}));
adjoint.reverseFile(out(adjFileNames{fieldXYZ}));
adjDE = OutputFile(out(adjRevFileNames{fieldXYZ}));

assert(fwdDE.numFramesAvailable() == adjDE.numFramesAvailable());
numT = fwdDE.numFramesAvailable();

readAheadBytes = 100e6;
readAheadValues = readAheadBytes / 8;
chunkFrames = min(numT, ...
    floor(readAheadValues / (fwdDE.FrameSize  + adjDE.FrameSize)));

%%
numLags = 10; % one lag means static dielectric
adjBuffer = zeros(adjDE.FrameSize, numLags + chunkFrames);
fwdBuffer = zeros(fwdDE.FrameSize/2, 2, chunkFrames);

%fwdEBuffer = zeros(fwdDE.FrameSize/2, chunkFrames);
%fwdDBuffer = zeros(fwdDE.FrameSize/2, chunkFrames);

bufferIndex = @(t) mod( (t-1), numLags + chunkFrames ) + 1;

fwdDE.open();
adjDE.open();

%%

% Pre-load the adjoint buffer (numLags frames)
tAdjFirst = 1;
tAdjLast = min(numT, tAdjFirst + numLags - 1);

adjBuffer(:,bufferIndex(tAdjFirst:tAdjLast)) = adjDE.readFrames(...
    'NumFrames', tAdjLast - tAdjFirst + 1, 'Regions', 'Together');

%%

%warning('Only handling one lag at a time!');

tBeginChunk = 1;

while tBeginChunk <= numT
    fprintf('Processing chunks %i to %i\n', tBeginChunk, ...
        tBeginChunk+chunkFrames-1);
    tEndChunk = min(numT, tBeginChunk + chunkFrames - 1);
    
    % First we load the next frame into the buffer.
    tFwdFirst = tBeginChunk;
    tFwdLast = min(numT, tFwdFirst + chunkFrames - 1);
    loadFwdFrames = tFwdLast - tFwdFirst + 1;
    
    tAdjFirst = tFwdFirst + numLags;
    tAdjLast = min(numT, tAdjFirst + chunkFrames - 1);
    loadAdjFrames = numel(tAdjFirst:tAdjLast);
    
    fwdBuffer(:,:,1:loadFwdFrames) = fwdDE.readFrames(...
        'NumFrames', loadFwdFrames, 'Regions', 'Together');
    fwdEBuffer = squish(fwdBuffer(:,1,:), 2); % squeeze out the {E,D} dimension
    fwdDBuffer = squish(fwdBuffer(:,2,:), 2);
    
    %fprintf('Fwd: load %i to %i\n', tFwdFirst, tFwdLast);
    
    if tAdjFirst <= tAdjLast
        adjBuffer(:,bufferIndex(tAdjFirst:tAdjLast)) = adjDE.readFrames(...
            'NumFrames', loadAdjFrames, 'Regions', 'Together');
        %fprintf('Adj: load %i to %i\n\n', tAdjFirst, tAdjLast);
    end
    
    %lagsProvided = length(coeffs{movableVert,freeDir}.tensor{fieldXYZ, fieldXYZ}.DB);
    %lagsProvided = 1;
    lagsProvided = 10;
    for lag = 1:lagsProvided
        
        tAdj = [tFwdFirst + lag - 1, ...
            min(numT, tFwdFirst + lag + chunkFrames - 2)];
        tFwd = [tFwdFirst, tFwdFirst + tAdj(2) - tAdj(1)];
        
        adjE = adjBuffer(:,bufferIndex(tAdj(1):tAdj(2)));
        fwdD = fwdDBuffer(:, (tFwd(1):tFwd(2))-tFwd(1)+1);
        fwdE = fwdEBuffer(:, (tFwd(1):tFwd(2))-tFwd(1)+1);
        
        ED = sum(adjE.*fwdD, 2);
        EE = sum(adjE.*fwdE, 2);
    
        for movableVert = verts
        for freeDir = 1:3
        if isstruct(coeffs{movableVert, freeDir})
        if isstruct(coeffs{movableVert,freeDir}.tensor{1,1})
        if length(coeffs{movableVert,freeDir}.tensor) >= fieldXYZ
        if isstruct(coeffs{movableVert,freeDir}.tensor{fieldXYZ,fieldXYZ})
        if length(coeffs{movableVert,freeDir}.tensor{fieldXYZ, fieldXYZ}.DB) >= lag

            sumSensitivity = 0;
            
            coeffD = coeffs{movableVert,freeDir}.tensor{fieldXYZ, fieldXYZ}.DB{lag}.coefficients;
            indexD = coeffs{movableVert,freeDir}.tensor{fieldXYZ, fieldXYZ}.DB{lag}.indices;
            dotProdDE = -sum(coeffD.*ED(indexD));
            sumSensitivity = sumSensitivity + dotProdDE;
            
            coeffE = coeffs{movableVert,freeDir}.tensor{fieldXYZ, fieldXYZ}.EH{lag}.coefficients;
            indexE = coeffs{movableVert,freeDir}.tensor{fieldXYZ, fieldXYZ}.EH{lag}.indices;
            dotProdEE = sum(coeffE.*EE(indexE));
            sumSensitivity = sumSensitivity + dotProdEE;
            
            vertJacobian(freeDir + 3*(movableVert-1), 1) = ...
                vertJacobian(freeDir + 3*(movableVert-1), 1) ...
                + sumSensitivity;
            
        end
        end
        end
        end
        end
        end
        end
    end
    
    tBeginChunk = tBeginChunk + chunkFrames;
end

fwdDE.close();
adjDE.close();

end % for fieldXYZ = 1:3

%%

ts3 = [vertJacobian(1:3:end), vertJacobian(2:3:end), vertJacobian(3:3:end)];
vert3 = [vertices(1:3:end), vertices(2:3:end), vertices(3:3:end)];

dfdp = vertJacobian' * bigJacobian;

dfdv = vertJacobian';
dvdp = bigJacobian;
