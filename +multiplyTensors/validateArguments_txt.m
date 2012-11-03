function validateArguments_txt(A, szA, B, szB, inA, inB, replaceDims)

% Make sure that the user is not lying about the sizes of A and B.
% We can compare the sizes up through the last non-singleton dimension.
% Scalars are size [1 1] in Matlab and size [1] (or []?) in tensor-land, so
% sometimes there are no indices to compare at all...

if ~isequal(stripSingletons(szA), stripSingletons(size(A)))
	errString = ['szA is not compatible with size(A):\n', ...
    	sprintf('\tszA = %s\nsize(A) = %s', num2str(szA), num2str(size(A)))];
    error(errString);
end

if ~isequal(stripSingletons(szB), stripSingletons(size(B)))
	errString = ['szA is not compatible with size(A):\n', ...
    	sprintf('\tszA = %s\nsize(A) = %s', num2str(szA), num2str(size(A)))];
    error(errString);
end

if ~areUnique(inA)
    error('Some elements of inA are repeating:\n\tinA = %s', ...
        num2str(inA));
end

if ~areUnique(inB)
    error('Some elements of inB are repeating:\n\tinB = %s', ...
        num2str(inB));
end

if length(inA) ~= length(inB)
    errString = ['Number of contracting dimensions is inconsistent:\n', ...
        sprintf('\tnumel(inA) = %i\n\tnumel(inB) = %i', ...
            length(inA), length(inB))];
    error(errString);
end

if ~isequal(szA(inA), szB(inB))
    errString = ['Contracting dimensions are not all of equal size: \n', ...
        sprintf('\tszA(inA) = %s\n\tszB(inB) = %s', num2str(szA(inA)), ...
        num2str(szB(inB)))];
    error(errString);
end

if ~isempty(replaceDims)
    if length(replaceDims) ~= length(inB)
        errString = ['If replaceDims is specified, each contracting'...,
            ' index of B must have a replacement dimension:\n', ...
            sprintf('\tlength(replaceDims) = %i\n\tlength(inB) = %i', ...
                length(replaceDims), length(inB))];
        error(errString);
    end
    
end

if ~isempty(replaceDims)
    if numel(replaceDims) ~= numel(inB)
        error('numel(inB) must equal numel(replaceDims)');
    end
    
    if ~areUnique(replaceDims)
        error('Some elements of replaceDims are repeated:\n\t%s',...
            replaceDims);
    end
    
    if ~areUnique([inB replaceDims])
        error(['Some elements of replaceDims overlap inB:\n', ...
            '\tinB = %s\n\treplaceDims = %s'], inB, replaceDims);
    end
end


function sz = stripSingletons(szIn)

[~, ii] = find(szIn > 1, 1, 'last');
sz = szIn(1:ii);


function b = areUnique(A)

b = numel(unique(A)) == numel(A);
