%% Test the tensor product function!

import multiplyTensors.*

isClose = @(A, B) norm(A(:) - B(:))/(0.5*norm(A(:)) + 0.5*norm(B(:))) < 1e-6;

%% Tensor size

A = rand(1);
assert(isequal(tsize(A, 1), 1));
assert(isequal(tsize(A, 2), [1 1]));

A = rand(3,3);
assert(isequal(tsize(A, 2), [3 3]));
assert(isequal(tsize(A, 3), [3 3 1]));

%% Tensor dimensions

sz = tensorOrder(1, 1, [], [], []);
assert(isequal(sz, [1 1]));

sz = tensorOrder([2 3], [2 2], 1, 1);
assert(isequal(sz, [3 2]));

sz = tensorOrder([2 3], [2 2], 1, 1, 2);
assert(isequal(sz, [2 3]));

sz = tensorOrder([10 1 1], [8 1 1]);
assert(isequal(sz, [10 1 1 8 1 1]));

sz = tensorOrder([10 1 1], [8 1 1], 2:3, 2:3);
assert(isequal(sz, [10 8]));

fprintf('Tensor size PASSED\n');

%% Pure outer product of two 2x2 matrices

A = rand(3,3);
B = rand(2,2);
[C szC] = txt(A, ndims(A), B, ndims(B));

assert(isequal(size(C), [size(A) size(B)]));

for ii = 1:size(A, 1)
    for jj = 1:size(A, 2)
        for kk = 1:size(B, 1)
            for ll = 1:size(B, 2)
                assert(C(ii,jj,kk,ll) == A(ii,jj)*B(kk,ll));
            end
        end
    end
end

fprintf('Outer product test PASSED\n');


%% Frobenius norm

A = rand(3,3,3,3);

A_frobenius = txt(A, ndims(A), A, ndims(A), 1:4, 1:4);

assert(isClose(A_frobenius, sum(A(:).^2)));

fprintf('Frobenius norm test PASSED\n');

%% Inner product of tensors

A = rand(3,3,2);
B = rand(3,3,2);

AB = txt(A, ndims(A), B, ndims(B), 1:3);

assert(isClose(AB, dot(A(:), B(:))));
fprintf('Inner product test PASSED\n');

%% Mixed inner-outer product of tensors: 2 outer indices, 1 inner index

A = rand(3,2,4);
B = rand(4,1,3);

AB = txt(A, ndims(A), B, ndims(B), 1, 3);

assert(isequal(size(AB), [2 4 4]));

% Do it by brute force

AoxB = txt(A, ndims(A), B, ndims(B)); % outer product

AB2 = zeros(2,4,4);

for ii = 1:size(A,2)
    for jj = 1:size(A,3)
        for kk = 1:size(B,1)
            for ll = 1:size(B,2)
                subA = A(:,ii,jj);
                subB = B(kk,ll,:);
                AB2(ii,jj,kk,ll) = dot(subA(:), subB(:));
            end
        end
    end
end

assert(isequal(AB2, AB));
fprintf('Inner-1 outer-2 product test PASSED\n');

%% Mixed inner-outer product of tensors, 2 inner indices

A = rand(3,2,4);
B = rand(1,4,3,2);

AB = txt(A, ndims(A), B, ndims(B), 1:2, 3:4);

assert(isequal(size(AB), [4 1 4]));

% Do it by brute force

AB2 = zeros(4,1,4);
for ii = 1:size(A,3)
    for jj = 1:size(B,1)
        for kk = 1:size(B,2)
            subA = A(:,:,ii);
            subB = B(jj,kk,:,:);
            AB2(ii,jj,kk) = dot(subA(:), subB(:));
        end
    end
end

assert(isequal(AB, AB2));
fprintf('Inner-2 product test PASSED\n');

%% Big scalar times matrix with txt (regression test)

A = 1;
dimsA = 5;
inA = 4;

B = rand(2, 1);
inB = 2;
replace = 1;

[C szC] = txt(A, dimsA, B, ndims(B), inA, inB, replace);

assert(isequal(szC, [1 1 1 2 1]));
assert(isequal(A*B(:), C(:)));

fprintf('Big scalar test PASSED\n');


%% Matrix-vector product (a special case worth testing)

A = rand(3,3);
B = rand(3,1);

% Matrix times column vector
AB = txt(A, ndims(A), B, 2, 2, 1);
assert(isequal(size(AB), [3 1]));
assert(isequal(AB, A*B));
fprintf('Matrix-vector product test PASSED\n');

% Matrix times column vector, worked backwards (why not)
% This is not the same as B'*A.  I'm just re-ordering the output indices
% from A*B.  Fun, huh.
BA = txt(B, 2, A, ndims(A), 1, 2);
assert(isequal(size(BA), [1 3]));
assert(isClose(BA, (A*B)'));
fprintf('Vector-matrix product test PASSED\n');

%% Matrix times vector with txt (regression test)

A = eye(3);
B = rand(3,1);

AB = txt(A, 2, B, 1, 2, 1);

assert(isClose(AB, A*B));
fprintf('Matrix-vector product test with truly 1D vector PASSED\n');

%% Tensor-matrix product

T = rand(4,4,3,4);
B = rand(5,3);

[C szC] = txa(T, ndims(T), B, 3);

assert(isequal(szC, [4 4 5 4]));
assert(isequal(size(C), [4 4 5 4]));

for ii = 1:size(T,1)
    for jj = 1:size(T,2)
        for kk = 1:size(T,4)
            bt = B*reshape(T(ii,jj,:,kk), [], 1);
            assert(isequal(bt, squish(C(ii,jj,:,kk))));
        end
    end
end

fprintf('Tensor-matrix product test PASSED\n');

%% Tensor-matrix product, matrix is Nx1, sparse (regression test)
% The point is that the sparse matrix has a trailing singleton dim.

T = rand([1 6 6 1 2]);
B = sparse(3, 1, 1);

[C szC] = txa(T, 5, B, 1, 2);

assert(isequal(szC, size(C)));

%% Tensor field times tensor field, no replacement

A = rand(10,7,3,2,5);   % [x y i j z]    free j
B = rand(5,3,7,10,4);   % [z i y x k]    k replaces i

C = tfxtf(A, ndims(A), [1 2 5], B, ndims(B), [4 3 1], ...
    3, 2); % 5

assert(isequal(size(C), [10 7 2 5 4]));

for ii = 1:size(A,1)
    for jj = 1:size(A,2)
        for kk = 1:size(A,5)
            a = squish(A(ii,jj,:,:,kk));
            b = squish(B(kk,:,jj,ii,:));
        
            c = txt(a, ndims(a), b, ndims(b), 1, 1); % 2
            
            assert(isequal(c, squish(C(ii, jj, :, kk, :))));
        end
    end
end

fprintf('Tensor field product, no replacement PASSED\n');

%%
C = tfxtf(A, ndims(A), [1 2 5], B, ndims(B), [4 3 1], ...
    3, 2, 5);

assert(isequal(size(C), [10 7 4 2 5]));

for ii = 1:size(A,1)
    for jj = 1:size(A,2)
        for kk = 1:size(A,5)
            a = squish(A(ii,jj,:,:,kk));
            b = squish(B(kk,:,jj,ii,:));
        
            c = txt(a, ndims(a), b, ndims(b), 1, 1, 2);
            
            assert(isequal(c, squish(C(ii, jj, :, :, kk))));
        end
    end
end

fprintf('Tensor field product with replacement PASSED\n');

%% Another test?  I was having an issue with interleaving dimensions properly.

A = rand(10, 11, 3);
B = rand(5, 10, 3);
C = tfxtf(A, ndims(A), 3, B, ndims(B), 3, 1, 2, 1);

assert(isequal(size(C), [5 11 3]));

fprintf('Tensor field product size test PASSED\n');

%% Tensor fields with leading singular dims

A = randi(9, 1, 3, 10);
B = randi(9, 1, 3, 10);

C = tfxtf(A, ndims(A), [1 3], B, ndims(B), [1 3]);

for xx = 1:size(A,1)
    for yy = 1:size(A,3);
        axy = squish(A(xx,:,yy));
        bxy = squish(B(xx,:,yy));
        AB = outerProduct(axy, 1, bxy, 1);
        
        cxy = squish(C(xx,:,yy,:));
        
        assert(isClose(AB, cxy));
    end
end

fprintf('Tensor field product with no inner products PASSED\n');

%% Tensor field times cell array of matrices
% This provides a way to use sparse matrices simply.

A = rand(100, 100, 3, 2);

e = ones(100, 1);
b = spdiags([e -2*e e], -1:1, 100, 100); % Matlab example: 2nd difference
B = repmat({b}, [1 3]);

C = txca(A, ndims(A), 3, B, 2, 2);

for nn = 1:3
    c = txa(A(:,:,nn,:), ndims(A), B{nn}, 2, 2);
    assert(isequal(C(:,:,nn,:), c));
end

fprintf('Tensor times cell array test PASSED\n');

%% TXCA regression test

A = rand([1 26 26 4 2]);
B = repmat({rand(1,2)}, [1 4]);
[C szC] = txca(A, 5, 4, B, 5, 2);

assert(isequal(size(C), [1 26 26 4]));
assert(isequal(szC, [1 26 26 4 1]));

%% Another test, about sizes...

A = rand(8, 30);
B = rand(31, 30);

% We consider A to be size [8 30 1], so ndims == 3.
[C szC] = txa(A, 3, B, 2, 2);
assert(isequal(szC, [8 31 1]));


%% Speed comparison
% Compare full tensor field multiplication to cell-sparse method.

Nx = 80;
Ny = 3000;
Nyy = 2990;

A = rand(Nx, Ny, 3);

% Apply a second-difference operation to the Y coordinate of tensor field
% A.  Furthermore resize the rows a little bit, from Ny to Nyy elements.
e = ones(Ny, 1);
b = spdiags([e -2*e e], -1:1, Nyy, Ny);

B_mat = repmat(full(b), [1 1 3]);  % [100 100 3]
B_cell = repmat({b}, [1 3]);

% Method 1: tfxtf
tic
C1 = tfxtf(A, ndims(A), 3, B_mat, ndims(B_mat), 3, 2, 2, 1);
tFull = toc;

%% Method 2: cells
tic
C2 = txca(A, ndims(A), 3, B_cell, 2, 2);
tSparse = toc;

assert(isClose(C1, C2));

fprintf('Transform size [%i %i %i] tensor to [%i %i %i] tensor:\n', ...
    Nx, Ny, 3, Nx, Nyy, 3);
fprintf('\t%2.3e seconds (full matrix method)\n', tFull);
fprintf('\t%2.3e seconds (cell array of sparse matrices)\n', tSparse);

%% All done.

fprintf('All tests done.\n');





