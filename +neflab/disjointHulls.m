function areDisjoint = disjointHulls(v1, v2)

min1 = min(v1);
max1 = max(v1);

min2 = min(v2);
max2 = max(v2);

areDisjoint = any( (min1 > max2) | (min2 > max1) );