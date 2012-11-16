function tt = naturalSamplingTimes(obj, duration, numSamples)

if isempty(duration)
    duration = obj.Durations{1}.Duration;
end

if isempty(numSamples)
    numSamples = ceil((duration(2)-duration(1)) ./ obj.Dt + 1);
end

if duration(2) == duration(1)
    numSamples = 1;
end

assert(numSamples >= 1);

tt = linspace(duration(1), duration(2), numSamples);
