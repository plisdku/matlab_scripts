function flux = fluxIntegral(vectorField, positions, normal)

sz = size(vectorField);
if sum(sz(1:3) > 1) == 3
    error('Cannot perform flux integral in 3D domain');
end

if ndims(vectorField) < 4
    error('Flux integrand must be a 4-D vector field');
end

normalComponent = vectorField(:,:,:,1,:)*normal(1) + ...
    vectorField(:,:,:,2,:)*normal(2) + ...
    vectorField(:,:,:,3,:)*normal(3);

for xyz = 3:-1:1 % might help it not lose indices
if sz(xyz) > 1
    if normal(xyz) ~= 0
        warning('Surface normal is not perpendicular to integration plane');
    end
    
    if length(positions{xyz}) ~= size(normalComponent, xyz)
        error('Positions must have compatible dimensions with normal component');
    end
    
    normalComponent = trapz(positions{xyz}, normalComponent, xyz);
end
end

flux = squeeze(normalComponent);
