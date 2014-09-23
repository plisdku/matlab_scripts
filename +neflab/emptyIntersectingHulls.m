function areDisjoint = emptyIntersectingHulls(v1, v2)
% Return true if the intersection of the bounding boxes of v1 and v2 has
% zero volume.  The bounding boxes may touch but not intersect.

min1 = min(v1);
max1 = max(v1);

min2 = min(v2);
max2 = max(v2);

areDisjoint = any( (min1 >= max2) | (min2 >= max1) );