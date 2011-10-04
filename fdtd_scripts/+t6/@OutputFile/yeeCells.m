function [ii, jj, kk] = yeeCells(obj, varargin)

X.Regions = 'Separate'; % Separate or Together
X = parseargs(X, varargin{:});

if strcmp(X.Regions, 'Separate')
    ii = cell(obj.numRegions(), 1);
    jj = cell(obj.numRegions(), 1);
    kk = cell(obj.numRegions(), 1);
    
    for rr = 1:length(obj.Regions)
        ii{rr} = obj.Regions.YeeCells(rr,1):obj.Regions.Stride(rr,1):obj.Regions.YeeCells(rr,4);
        jj{rr} = obj.Regions.YeeCells(rr,2):obj.Regions.Stride(rr,2):obj.Regions.YeeCells(rr,5);
        kk{rr} = obj.Regions.YeeCells(rr,3):obj.Regions.Stride(rr,3):obj.Regions.YeeCells(rr,6);
    end
    
    if obj.numRegions() == 1
        ii = ii{1};
        jj = jj{1};
        kk = kk{1};
    end
elseif strcmp(X.Regions, 'Together')
    ii = [];
    jj = [];
    kk = [];
    
    for rr = 1:obj.numRegions()
        [ii1, jj1, kk1] = unrollRegion(obj.Regions.YeeCells(rr,:));
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
