function [ii, jj, kk] = yeeCells(obj, varargin)

ii = cell(size(obj.Regions));
jj = cell(size(obj.Regions));
kk = cell(size(obj.Regions));

for rr = 1:length(obj.Regions)
        ii{rr} = obj.Regions{rr}.YeeCells(1):obj.Regions{rr}.Stride(1):obj.Regions{rr}.YeeCells(4);
        jj{rr} = obj.Regions{rr}.YeeCells(2):obj.Regions{rr}.Stride(2):obj.Regions{rr}.YeeCells(5);
        kk{rr} = obj.Regions{rr}.YeeCells(3):obj.Regions{rr}.Stride(3):obj.Regions{rr}.YeeCells(6);
end

if length(obj.Regions) == 1
    ii = ii{1};
    jj = jj{1};
    kk = kk{1};
end
