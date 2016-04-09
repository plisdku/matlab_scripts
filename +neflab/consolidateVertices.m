function vertices = consolidateVertices(originalVertices, vertices, consolidateDistance)
% consolidateVertices   Correct the exact position of vertices from a
% geometry operation
%
% v = consolidateVertices(vOld,vNew) will find rows of vNew which differ by
%   less than 1e-13 from rows of vOld, and set them exactly equal to the 
%   nearby values from vOld.
%
% v = consolidateVertices(vOld, vNew, epsilon) also sets the distance under
%   which vertices are identified as equal.
%
% Purpose: when performing Boolean operations on solids with
% finite-precision vertex coordinates, sometimes a vertex that should not
% have been changed at all by the operation will end up perturbed by a bit
% or two out around the machine epsilon level.  This can cause later
% Boolean operations to behave badly or result in duplicated vertices later
% on in Matlab.  (All this badness comes about from translating Matlab's
% floating-point coordinates into CGAL's exact arithmetic, and back.)
%
% A better technical solution would be to improve the NefLab C++ program to
% guarantee that if an output point is supposed to match an input point
% that it agrees down to the last bit.  This is doable but I'm being lazy.

if nargin < 3
    consolidateDistance = 1e-13;
end

numOriginal = size(originalVertices,1);
numNew = size(vertices,1);

for iVert = 1:numNew
    for iOrig = 1:numOriginal
        delta = vertices(iVert,:) - originalVertices(iOrig,:);
        
        if norm(delta) < consolidateDistance
            vertices(iVert, :) = originalVertices(iOrig,:);
            break;
        end
        
    end
end

