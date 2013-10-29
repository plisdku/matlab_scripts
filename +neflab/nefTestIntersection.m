function doIntersect = nefTestIntersection(v1, f1, v2, f2)
% nefTestIntersection    Determine whether polyhedra intersect
%

[vv ff] = neflab.nefIntersection(v1, f1, v2, f2);

doIntersect = numel(ff) > 0;

