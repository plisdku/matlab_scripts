function yees = yeeCells(obj, varargin)

X.Regions = 'Separate'; % Separate or Together
X.Region = [];
X = parseargs(X, varargin{:});

validateArguments(obj, X);

if isempty(X.Region)
    X.Region = 1;
end

yees = cell(1,3);
if strcmpi(X.Regions, 'Separate')
    for xyz = 1:3
        yees{xyz} = obj.Regions.YeeCells(X.Region,xyz) : ...
            obj.Regions.Stride(X.Region,xyz) : ...
            obj.Regions.YeeCells(X.Region, xyz+3);
    end
elseif strcmpi(X.Regions, 'Together')
    
    for rr = 1:obj.numRegions()
        [ii, jj, kk] = unrollRegion(obj.Regions.YeeCells(rr,:));
        yees{1} = [yees{1}, ii];
        yees{2} = [yees{2}, jj];
        yees{3} = [yees{3}, kk];
    end
end

%{
if strcmp(X.Regions, 'Separate')
    ii = cell(obj.numRegions(), 1);
    jj = cell(obj.numRegions(), 1);
    kk = cell(obj.numRegions(), 1);
    
    for rr = 1:obj.numRegions
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
%}

function validateArguments(obj, X)

if ~strcmpi(X.Regions, 'Separate') && ~strcmpi(X.Regions, 'Together')
    error('Regions must be Separate or Together');
end

if strcmpi(X.Regions, 'Separate')
    if obj.numRegions() > 1 && isempty(X.Region)
        error('More than one region is present in this file.  Please select one using the Region argument.');
    end
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
