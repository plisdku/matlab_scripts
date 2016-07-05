%% Read DF data from COMSOL
function [F, dFdp, integrals] = evaluateFDF(movableMesh, params, filename, minFaceArea)

if nargin < 3
    filename = 'DF_on_surfaces.txt';
end

if nargin < 4
    minFaceArea = 1e-12;
end

%%
fid = fopen(filename);
% maxwell version (old version), without normal vectors
%AA = cell2mat(textscan(fid, '%n%n%n%n', 'CommentStyle', '%'));
AA = cell2mat(textscan(fid, '%n%n%n%n%n%n%n', 'CommentStyle', '%'));
fclose(fid);

x = AA(:,1);
y = AA(:,2);
z = AA(:,3);
DF = AA(:,4);

%% Get the triangles from the movable ball dude.

%params = [0 0 0 0 0 0]';
numParams = numel(params);

m = movableMesh.meshes(params);

faces = m{1}.faces;
verts = m{1}.patchVertices;

jac = m{1}.jacobian;

%%

%figure(1); clf
%plot3(x, y, z, 'o');
%hold on
%flatPatch('Faces', faces, 'Vertices', verts, 'FaceColor', 'r')
%axis image

%% Now some integrals!!

linearInterp = TriScatteredInterp(x,y,z,DF, 'linear');
nearestInterp = TriScatteredInterp(x,y,z,DF, 'nearest');
maxDF = max(abs(DF(:)));

omitCount = 0;

numFaces = size(faces, 1);
integrals = zeros(numParams, numFaces);
for ff = 1:numFaces
    %fprintf('Face %i of %i\n', ff, numFaces);
    % Obtain the three corners of the triangle
    vx = verts(faces(ff,:),1);
    vy = verts(faces(ff,:),2);
    vz = verts(faces(ff,:),3);
    
    dv1 = verts(faces(ff,2),:) - verts(faces(ff,1),:);
    dv2 = verts(faces(ff,3),:) - verts(faces(ff,1),:);
    faceArea = 0.5*norm(cross(dv1, dv2));
    
    if faceArea < minFaceArea
    %if abs(det([vx vy vz])) < 1e-5 % wrong.  grr.
        %warning('Null triangle!!  So tiny!  So cute!');
        omitCount = omitCount+1;
        continue;
    end
    
    jacobianIndices = @(vertNum) (1:3) + 3*(vertNum-1);
    
    jacobianChunk = jac([jacobianIndices(faces(ff,1)), ...
        jacobianIndices(faces(ff,2)), ...
        jacobianIndices(faces(ff,3))],:);
    
    % Get the quadrature points and weights in 3D
    N = 3;
    try
        [xx, ww, nv, bc] = simplexQuad3d(N, [vx vy vz]);
    catch exc
        keyboard
    end
    
    %{
    figure(2); clf
    flatPatch('Faces', faces, 'Vertices', verts, ...
        'FaceColor', 'g', 'FaceAlpha', 0.1)
    hold on
    flatPatch('Faces', faces(ff,:), 'Vertices', verts, ...
        'FaceColor', 'r', 'FaceAlpha', 0.5);
    plot3(xx(1,:), xx(2,:), xx(3,:), 'o')
    axis image vis3d
    view(3)
    keyboard
    %}
    
    
    
    % Now get DF at those points with griddata.
    pointwiseDF = saferTriScatteredInterp(linearInterp, nearestInterp, ...
        xx, maxDF );
    
    if any(isnan(pointwiseDF(:)))
        warning('We have NaN.\n');
        keyboard
    end
    
    %Vq = Vq * nv(3); % DF*nz
    %Vq = abs(nv(3))*ones(size(Vq));
    
    % Dot the normal vector with the vertex jacobians.
    nvDotJ = [jacobianChunk(1:3,:)'*nv, ...
        jacobianChunk(4:6,:)'*nv, ...
        jacobianChunk(7:9,:)'*nv];
    
    % size is [numParams 3].
    
    %nvDotJ = [dot(nv, jacobianChunk(1:3,:)), ...
    %    dot(nv, jacobianChunk(4:6,:)), ...
    %    dot(nv, jacobianChunk(7:9,:))];
    
    % Interpolate nvDotJ to each integration position.  Use barycentric
    % interpolation.
    verticalPerturbation = nvDotJ * bc; % size is [numParams numQuadPoints]
    
    % Now perform the integral, using the perturbation.
    
    addends = bsxfun(@times, ww.*pointwiseDF, verticalPerturbation);
    
    if any(isnan(addends(:)))
        warning('We have NaN.\n');
        keyboard
    end
    
    integrals(:,ff) = sum(addends, 2);
    %fprintf('Integral = %s\n', num2str(integrals(ff)));
    
    if abs(integrals(:,ff)) > 1e26
        keyboard
    end
    
end

%%

% maxwell uses this 1e-27 factor
%dFdp = transpose(sum(integrals*1e-27, 2)); % there are three m->nm conversions.
dFdp = transpose(sum(integrals, 2));

%% F?

% maxwell uses column 2, I think column 1 is frequency
F = ll.evaluateF('Column', 1);


    
 
function vals = saferTriScatteredInterp(linearInterp, nearestInterp, xx, maxF)

vals = linearInterp(xx(1,:), xx(2,:), xx(3,:));

maskNaN = isnan(vals);
maskOutOfBounds = abs(vals) > maxF;

%{
if any(maskNaN)
    warning(sprintf('There were %i NaNs, using nearest-neighbor interpolation', ...
        sum(maskNaN(:))));
end

if any(maskOutOfBounds)
    warning(sprintf('There were %i out of bounds values, using nearest-neighbor interpolation', ...
        sum(maskOutOfBounds(:))));
end
%}

maskBad = maskNaN | maskOutOfBounds;

if any(maskBad)
    vals(maskBad) = nearestInterp( ...
        xx(1,maskBad), xx(2,maskBad), xx(3,maskBad));
end




