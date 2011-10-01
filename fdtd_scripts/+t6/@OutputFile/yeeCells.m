function [ii, jj, kk] = yeeCells(obj, varargin)

X.Regions = 'Separate'; % Separate or Together
X = parseargs(X, varargin{:});

if strcmp(X.Regions, 'Separate')
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
elseif strcmp(X.Regions, 'Together')
    ii = [];
    jj = [];
    kk = [];
    
    for rr = 1:length(obj.Regions)
        [ii1, jj1, kk1] = unrollRegion(obj.Regions{rr}.YeeCells);
        ii = [ii, ii1];
        jj = [jj, jj1];
        kk = [kk, kk1];
    end
else
    error('Invalid Regions option');
end


function [ii, jj, kk] = unrollRegion(yeeCells)

ii = [];
jj = [];
kk = [];

allX = yeeCells(1):yeeCells(4);

for zz = yeeCells(3):yeeCells(6)
for yy = yeeCells(2):yeeCells(5)
    ii = [ii, allX];
    jj = [jj, yy*ones(size(allX))];
    kk = [kk, zz*ones(size(allX))];
end
end
