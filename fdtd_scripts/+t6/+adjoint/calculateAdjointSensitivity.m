function [dfdp, dfdv, dvdp] = calculateAdjointSensitivity(node, parameters)

import t6.*
import t6.adjoint.*

fwdFileNames = {'boundary_dexx', 'boundary_deyy', 'boundary_dezz'};
adjFileNames = {'boundary_adjoint_exx', 'boundary_adjoint_eyy', ...
    'boundary_adjoint_ezz'};
adjRevFileNames = {'boundary_adjoint_exx.rev', 'boundary_adjoint_eyy.rev',...
    'boundary_adjoint_ezz.rev'};

[coeffs, verts] = readDeps('output/Depsilon');

bigJacobian = node.jacobian(parameters);
vertices = node.vertices(parameters);
vertJacobian = sparse(size(vertices, 1), 1);

%%

for fieldXYZ = 1:3

fprintf('Reversing file, field %i of 3\n', fieldXYZ);

fwdDE = OutputFile(['output/', fwdFileNames{fieldXYZ}]);
adjoint.reverseFile(['output/', adjFileNames{fieldXYZ}]);
adjDE = OutputFile(['output/', adjRevFileNames{fieldXYZ}]);

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
    
    
    
    for movableVert = verts
    for freeDir = 1:3
    if isstruct(coeffs{movableVert, freeDir})
    if isstruct(coeffs{movableVert,freeDir}.tensor{1,1})
    if length(coeffs{movableVert,freeDir}.tensor) >= fieldXYZ
    if isstruct(coeffs{movableVert,freeDir}.tensor{fieldXYZ,fieldXYZ})
        
        sumSensitivity = 0;
        
        lagsProvided = length(coeffs{movableVert,freeDir}.tensor{fieldXYZ, fieldXYZ}.DB);
        for lag = 1:lagsProvided
            
            tAdj = [tFwdFirst + lag - 1, ...
                min(numT, tFwdFirst + lag + chunkFrames - 2)];
            tFwd = [tFwdFirst, tFwdFirst + tAdj(2) - tAdj(1)];
            
            coeffD = coeffs{movableVert,freeDir}.tensor{fieldXYZ, fieldXYZ}.DB{lag}.coefficients;
            indexD = coeffs{movableVert,freeDir}.tensor{fieldXYZ, fieldXYZ}.DB{lag}.indices;
            
            adjE = adjBuffer(indexD,bufferIndex(tAdj(1):tAdj(2)));
            fwdD = fwdDBuffer(indexD, (tFwd(1):tFwd(2))-tFwd(1)+1);
            
            dotProdDE = -sum(coeffD.*sum(adjE.*fwdD, 2));
            absDotDE = sum(abs(coeffD).*sum(abs(adjE).*abs(fwdD), 2));
            
            sumSensitivity = sumSensitivity + dotProdDE;
            
            %{
            if norm(adjE(:)) ~= 0
                figure(1); clf
                plot(adjE); hold on; plot(fwdD)
                title(sprintf('E-E %s along %s, v%i', char('w'+fieldXYZ), char('w'+freeDir), movableVert));
                fprintf('Dot DE = %2.4g, abs dot %2.4f\n', dotProdDE, absDotDE);
                fprintf('\tcoeffD = %2.4f\n', coeffD);
                pause
            end
            %}
            
        end
        
        lagsProvided = length(coeffs{movableVert,freeDir}.tensor{fieldXYZ, fieldXYZ}.EH);
        for lag = 1:lagsProvided
            
            tAdj = [tFwdFirst + lag - 1, ...
                min(numT, tFwdFirst + lag + chunkFrames - 2)];
            tFwd = [tFwdFirst, tFwdFirst + tAdj(2) - tAdj(1)];
            
            coeffE = coeffs{movableVert,freeDir}.tensor{fieldXYZ, fieldXYZ}.EH{lag}.coefficients;
            indexE = coeffs{movableVert,freeDir}.tensor{fieldXYZ, fieldXYZ}.EH{lag}.indices;
            
            adjE = adjBuffer(indexE,bufferIndex(tAdj(1):tAdj(2)));
            fwdE = fwdEBuffer(indexE, (tFwd(1):tFwd(2))-tFwd(1)+1);
            
            dotProdEE = sum(coeffE.*sum(adjE.*fwdE, 2));
            absDotEE = sum(abs(coeffE).*sum(abs(adjE).*abs(fwdE), 2));
            
            sumSensitivity = sumSensitivity + dotProdEE;
            
            %{
            if norm(adjE(:)) ~= 0
                figure(2); clf
                plot(adjE); hold on; plot(fwdE)
                title(sprintf('E-E %s along %s, v%i', char('w'+fieldXYZ), char('w'+freeDir), movableVert));
                fprintf('Dot EE = %2.4g, abs dot %2.4f\n', dotProdEE, absDotEE);
                fprintf('\tcoeffE = %2.4f\n', coeffE);
                pause
            end
            %}
        end
        
        vertJacobian(freeDir + 3*(movableVert-1), 1) = ...
            vertJacobian(freeDir + 3*(movableVert-1), 1) ...
            + sumSensitivity;
    
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
