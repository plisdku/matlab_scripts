function Eout = transformField(T, dx, x, y, z, E, xx, yy, zz)
%
% Apply coordinate transformation T followed by a translation dx to fields in E
%
% transformField(T, translation, x, y, z, E, xx, yy, zz)
%
% 

import modeInjection.*

E = reshape(E, numel(x), numel(y), numel(z), 3);
Eout = zeros(numel(xx), numel(yy), numel(zz), 3);

[xIn yIn zIn] = ndgrid(x, y, z);
[xOut yOut zOut] = ndgrid(xx - dx(1), yy - dx(2), zz - dx(3));

xyzIn = cat(4, xIn, yIn, zIn);
xyzOut = cat(4, xOut, yOut, zOut);

xyzInterp = t6.multTensor(xyzOut, inv(T), 4);

Erot = t6.multTensor(E, T, 4);

for dim = 1:3
    Eout(:,:,:,dim) = interpn(xIn, yIn, zIn, Erot(:,:,:,dim), ...
        xyzInterp(:,:,:,1), xyzInterp(:,:,:,2), xyzInterp(:,:,:,3));
    if any(isnan(Eout(:)))
        keyboard
    end
end





