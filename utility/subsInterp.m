function B = subsInterp(A, samplePoints, interpDims)
% B = subsInterp(A, subscriptArray)
% B = subsInterp(A, subscriptArray, interpDims)
% 
% Quick array interpolation.  If size(A) = [10 10] then this function will
% interpolate to find values A(x,y) for any x, y between 1 and 10.

if nargin == 2
    interpDims = 1:ndims(A);
    keepDims = [];
else
    % setdiff is so bloody slow!
    %keepDims = setdiff(1:ndims(A), interpDims);
    keepDims = true(1, ndims(A));
    keepDims(interpDims) = false;
    keepDims = find(keepDims);
end

nInterp = numel(interpDims);
nKeep = numel(keepDims);
nDims = nInterp + nKeep;
nSamples = size(samplePoints, 1);

szA = size(A);

szInterp = szA(interpDims);
szKeep = szA(keepDims);

for dd = 1:size(samplePoints,2)
    assert(all(samplePoints(:,dd) >= 1 & ...
        samplePoints(:,dd) <= szA(interpDims(dd))));
end

% Index combinations
lowerSubscript = floor(samplePoints);
upperSubscript = ceil(samplePoints);
subscripts = permute(cat(3, lowerSubscript, upperSubscript), [1 3 2]);
%assert(equalSize(size(subscripts), [nSamples 2 nInterp]));

% Weight combinations
upperWeight = permute(samplePoints - lowerSubscript, [1 3 2]);
weights = cat(2, 1 - upperWeight, upperWeight);
%assert(equalSize(size(weights), [nSamples 2 nInterp]));

% What I have now:
%
% subscripts(samplePoints, twoIndices, dimensions)
% weights(samplePoints, twoIndices, dimensions)

% For a single sample point, subscripts(n, :, :) is
%   [x0 y0 z0 w0
%    x1 y1 z1 w1]
% in 4D in this case.  To generate all combinations of these coords,
%   (x0 y0 z0 w0), (x0 y0 z0 w1), (x0 y0 z1 w0), ...
% use an array of indices,
%   [ 1 1 1 1
%     1 1 1 2
%     1 1 2 1
%     1 1 2 2
%     . . .   ].
% This can be generated with allWords().  But to get my desired array
%  [x0 y0 z0 w0
%   x0 y0 z0 w1
%   x0 y0 z1 w0 ...]
% I need to reach into the "indices" array (above) and grab elements
%  [1 3 5 7
%   1 3 5 8
%   1 3 6 7 ...]
% which can be had from sub2ind:

words = allWords(2, nInterp) + 1;
nWords = size(words,1);
%assert(equalSize(size(words), [2^nInterp, nInterp]));

ll = sub2ind([2 nInterp], words, repmat(1:nInterp, nWords, 1));
%assert(equalSize(size(ll), size(words)));

% Now ll are linear indices into the arrays of 
% 1) index combinations, [x0 y0; x1 y1], and
% 2) weight combinations, [wx0 wy0; wx1 wy1].
%
% The value of each sample in turn is
%
% A(x0,y0)*wx0*wy0 + A(x1,y0)*wx1*wy0 + A(x0,y1)*wx0*wy1 + A(x1,y1)*wx1*wy1
%

sampleSubscripts = reshape(subscripts(:,ll), nSamples, [], nInterp);
%assert(equalSize(size(sampleSubscripts), [nSamples, nWords, nInterp]));
% for each sample:
% sampleSubscripts(n,:,:) = [x0 y0
%                            x0 y1
%                            x1 y0
%                            x1 y1]

sampleCoordWeights = reshape(weights(:,ll), nSamples, [], nInterp);
%assert(equalSize(size(sampleCoordWeights), [nSamples, nWords, nInterp]));
% for each sample:
% sampleCoordWeights(n,:,:) = [wx0 wy0
%                              wx0 wy1
%                              wx1 wy0
%                              wx1 wy1]

ind = vecSub2Ind(szInterp, reshape(sampleSubscripts, [], nInterp));
sampleIndices = reshape(ind, nSamples, []);
%assert(equalSize(size(sampleIndices), [nSamples nWords]));
% for each sample:
% sampleIndices(n,:) = [ indexOf(x0,y0) indexOf(x0, y1) ... ]
% where indexOf() is from sub2ind

sampleWeights = reshape(prod(sampleCoordWeights,3), ...
    [ones(1,nKeep), nSamples, nWords]);
% the ones(1,nKeep) is so I don't need to shiftdim it later.
%assert(equalSize(size(sampleWeights), [nSamples nWords]));
% for each sample:
% sampleWeights(n,:) = [ wx0*wy0 wx0*wy1 wx1*wy0 wx1*wy1 ]

% The operation now is essentially B_out = sum(B(pointIndices).*weights).
% pointIndices is the right size now, and we need to get the corresponding
% values of B.

%% Create a sparse matrix to do this operation.

numInterpElements = prod(szInterp);
numKeepElements = prod(szKeep);

rows = repmat((1:1000)', 1, 8);
interpMatrix = sparse(rows(:), sampleIndices(:), sampleWeights(:), ...
    nSamples, prod(szInterp));

C = permute(A, [interpDims keepDims nDims+1]);
C = reshape(C, numInterpElements, []);

B = interpMatrix*C;
B = reshape(B, [nSamples szKeep 1]);

%{
B = permute(A, [keepDims interpDims nDims+1]);

%indicesB = repmat({':'}, 1, numel(keepDims)+1);
%indicesB{end} = sampleIndices;
indicesB = {':',':',':',':',':',':',':',':',':',':'};
indicesB{nKeep+1} = sampleIndices;

Bsubs = B(indicesB{1:nKeep+1});
B = reshape(Bsubs, [szKeep, size(sampleIndices)]);
%B = reshape(B(indicesB{1:nKeep+1}), [szKeep, size(sampleIndices)]);

% Size of B is now [szKeep nSamples nWords].  Sum over the nWords.
szB = size(B);
if nKeep > 0
    assert(isequal(szB(1:nKeep), szKeep));
end
%assert(isequal(szB(nKeep+1), nSamples));
%assert(isequal(szB(nKeep+2), size(sampleIndices,2)));

B_addends = bsxfun(@times, B, sampleWeights); %shiftdim(sampleWeights, -nKeep));

B = sum(B_addends, nKeep+2);
B = permute(B, [ndims(B), 1:ndims(B)-1]); % put interpolation at the front
%}

end



function truth = equalSize(szA, szB)

N = min(length(szA), length(szB));

truth = isequal(szA(1:N), szB(1:N)) && ...
    all(szA(N+1:end) == 1) && ...
    all(szB(N+1:end) == 1);

end












