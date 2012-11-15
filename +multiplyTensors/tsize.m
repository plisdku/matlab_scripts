function sz = tsize(A, numDims)

if numDims == 1
    if size(A,2) ~= 1
        error('size(A) (%s) is incompatible with numDims (%s)', ...
            size(A), numDims);
    end
    
    sz = numel(A);
else
    
    if numDims < ndims(A)
        error('size(A) (%s) is incompatible with numDims(%s)', ...
            size(A), numDims);
    end
    
    sz = [size(A), ones(1, numDims - ndims(A))];
end