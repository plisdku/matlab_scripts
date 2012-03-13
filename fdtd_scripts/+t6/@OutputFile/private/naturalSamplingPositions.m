function pos = naturalSamplingPositions(obj, bounds, region, numSamples)

pos = cell(1,3);

if isempty(bounds)
    bounds = obj.Regions.Bounds(region,:);
end

if isempty(numSamples)
    numSamples = (bounds(4:6) - bounds(1:3)) ./ obj.Dxyz + 1;
end

numSamples = ceil(numSamples);
numSamples(numSamples < 1) = 1;
numSamples(obj.Regions.YeeCells(region,4:6) == obj.Regions.YeeCells(region,1:3)) = 1;

for xyz = 1:3
    if numSamples(xyz) > 1
        pos{xyz} = linspace(bounds(xyz), bounds(xyz+3), numSamples(xyz));
    else
        pos{xyz} = bounds(xyz);
    end
end    
