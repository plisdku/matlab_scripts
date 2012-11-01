%% Test the tensor product function!

import multiplyTensors.*

isClose = @(A, B) norm(A(:) - B(:))/(0.5*norm(A(:)) + 0.5*norm(B(:))) < 1e-6;

%% Pure outer product of two 2x2 matrices

A = rand(3,3);
B = rand(2,2);

C = txt(A, B);

assert(isequal(size(C), [size(A) size(B)]))

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

A_frobenius = txt(A, A, 1:4);

assert(isequal(A_frobenius, sum(A(:).^2)));

fprintf('Frobenius norm test PASSED\n');

%% Inner product of tensors

A = rand(3,3,2);
B = rand(3,3,2);

AB = txt(A, B, 1:3);

assert(isequal(AB, dot(A(:), B(:))));
fprintf('Inner product test PASSED\n');

%% Mixed inner-outer product of tensors: 2 outer indices, 1 inner index

A = rand(3,2,4);
B = rand(4,1,3);

AB = txt(A, B, 1, 3);

assert(isequal(size(AB), [2 4 4]));

% Do it by brute force

AoxB = txt(A, B); % outer product

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

AB = txt(A, B, 1:2, 3:4);

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

%% Matrix-vector product (a special case worth testing)

A = rand(3,3);
B = rand(3,1);

% Matrix times column vector
AB = txt(A, B, 2, 1);
assert(isequal(size(AB), [3 1]));
assert(isequal(AB, A*B));
fprintf('Matrix-vector product test PASSED\n');

% Matrix times column vector, worked backwards (why not)
% This is not the same as B'*A.  I'm just re-ordering the output indices
% from A*B.  Fun, huh.
BA = txt(B, A, 1, 2);
assert(isequal(size(BA), [1 3]));
assert(isClose(BA, (A*B)'));
fprintf('Vector-matrix product test PASSED\n');

%% Tensor-matrix product

T = rand(4,4,3,4);
B = rand(5,3);

[C szC] = txa(T, B, 3);

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

%% Tensor field times tensor field

A = rand(10,11,3,2,2);   % [x y i j z]    free j
B = rand(2,3,11,10,4);   % [z i y x k]    k replaces i

C = tfxtf(A,B,[1 2 5],[4,3,1],3,2,5);

assert(isequal(size(C), [10 11 4 2 2]));

for ii = 1:size(A,1)
    for jj = 1:size(A,2)
        for kk = 1:size(A,5)
            a = squish(A(ii,jj,:,:,kk));
            b = squish(B(kk,:,jj,ii,:));
        
            c = txt(a, b, 1, 1, 2);
            
            assert(isequal(c, squish(C(ii, jj, :, :, kk))));
        end
    end
end

fprintf('Tensor field product test PASSED\n');

%% Another test?  I was having an issue with interleaving dimensions properly.

A = rand(10, 11, 3);
B = rand(5, 10, 3);
C = tfxtf(A, B, 3, 3, 1, 2, 1);

assert(isequal(size(C), [5 11 3]));

fprintf('\tTensor field product size test PASSED\n');

%% Tensor field times cell array of matrices
% This provides a way to use sparse matrices simply.

A = rand(100, 100, 3, 2);

e = ones(100, 1);
b = spdiags([e -2*e e], -1:1, 100, 100); % Matlab example: 2nd difference
B = repmat({b}, [1 3]);

C = txca(A, B, 2, 3);

for nn = 1:3
    c = txa(A(:,:,nn,:), B{nn}, 2);
    assert(isequal(C(:,:,nn,:), c));
end

fprintf('Tensor times cell array test PASSED\n');

        
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
C1 = tfxtf(A, B_mat, 3, 3, 2, 2, 1);
tFull = toc;

% Method 2: cells
tic
C2 = txca(A, B_cell, 2, 3);
tSparse = toc;

assert(isClose(C1, C2));

fprintf('Transform size [%i %i %i] tensor to [%i %i %i] tensor:\n', ...
    Nx, Ny, 3, Nx, Nyy, 3);
fprintf('\t%2.3e seconds (full matrix method)\n', tFull);
fprintf('\t%2.3e seconds (cell array of sparse matrices)\n', tSparse);