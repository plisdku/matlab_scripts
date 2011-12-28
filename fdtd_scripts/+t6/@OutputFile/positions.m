function [xx, yy, zz] = positions(obj, varargin)

xx = cell(obj.numRegions(), length(obj.Fields));
yy = xx;
zz = xx;

[ii jj kk] = obj.yeeCells();

if ~iscell(ii)
    ii = {ii};
    jj = {jj};
    kk = {kk};
end

for rr = 1:obj.numRegions()
for ff = 1:length(obj.Fields)
    offset = obj.Fields{ff}.Offset;
    xx{rr,ff} = obj.Dxyz(1)*(ii{rr} + offset(1)) + obj.Origin(1);
    yy{rr,ff} = obj.Dxyz(2)*(jj{rr} + offset(2)) + obj.Origin(2);
    zz{rr,ff} = obj.Dxyz(3)*(kk{rr} + offset(3)) + obj.Origin(3);
end
end

if length(xx) == 1
    xx = xx{1};
    yy = yy{1};
    zz = zz{1};
end
