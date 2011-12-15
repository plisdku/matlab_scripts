function sourcePhasors = makeSourcePhasors(obj)

source = obj.source;
%{
if coarsening > 1
    numSourceTimesteps = size(obj.source,2);
    numSamples = length(1:coarsening:numSourceTimesteps+coarsening-1);
    source = zeros(4, numSamples);
    sampledSource = downsample(obj.source', coarsening)';
    source(:,1:size(sampledSource,2)) = sampledSource;
    assert(length(source) == numSamples);
end
%}

sourcePhasors = ifft(source, [], 2);

%plot(abs(sourcePhasors(1,:)));
%pause