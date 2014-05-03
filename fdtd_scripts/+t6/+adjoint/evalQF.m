function [f Df] = evalQF(data, filters, kernel, xs, ys, zs, ts)
% [f Df] = evalQF(data, filters, kernel, xs, ys, zs, ts)
%

import multiplyTensors.*
import t6.adjoint.*

% Evaluate f

d = data;

for mm = 1:numel(filters)
    
    %if isa(filters{mm}.Data, 'function_handle')
    %    U = applyCoordinates(filters{mm}.Data, filters{mm}.dim, ...
    %        xs, ys, zs, ts);
    %else
    %    U = filters{mm}.Data;
    %end
    
    if any(strcmpi(filters{mm}.Operation, {'Matrix', 'MatrixArray'}))
        
        U = applyCoordinates(filters{mm}.Data, filters{mm}.dim, ...
            xs, ys, zs, ts);
        
        if iscell(U)        % Cell array of matrices
            d = txca(d, 5, 4, U, filters{mm}.dim, 2);
        else                % Multiply everything by one matrix
            d = txa(d, 5, U, filters{mm}.dim);
        end
    
    elseif strcmpi(filters{mm}.Operation, 'Pointwise')
        
        U = applyCoords_Multiply(filters{mm}.Data, filters{mm}.dim, ...
            xs, ys, zs, ts);
        
        %d = bsxfun(@times, shape(U, size(d), filters{mm}.dim), d);
        d = bsxfun(@times, U, d);
        
    end
    
    %if mm == 4
    %    assignin('base', 'qfData', d);
    %end
    
end

%assignin('base', 'qfKet', U);
%assignin('base', 'qfd', d);
%assignin('base', 'qfts', ts);

% At this point any remaining nonsingular dimensions of d must be eaten by
% a quadratic form, as d' * K * d, to obtain f.  The kernel needs to do
% this.

% Build up this kernel term in stages, then?
% I also will need transpose(k)*d, so build that one too.

kd = d;
ktd = d;

if ~iscell(kernel) && ~isempty(kernel)
    kernel = {kernel};
end
for mm = 1:numel(kernel)
    
    %if isa(kernel{mm}.Data, 'function_handle')
    %    U = applyCoordinates(kernel{mm}.Data, kernel{mm}.dim, ...
    %        xs, ys, zs, ts);
    %else
    %    U = kernel{mm}.Data;
    %end
    
    if any(strcmpi(kernel{mm}.Operation, {'Matrix', 'MatrixArray'}))
        U = applyCoordinates(kernel{mm}.Data, kernel{mm}.dim, ...
            xs, ys, zs, ts);
        
        if iscell(U)
            error('Not supporting this now, what would it mean?');
            %kd = txca(kd, 5, 4, U, kernel{mm}.dim, 2);
            %ktd = txca(ktd, 5, 4, U, kernel{mm}.dim, 2);
        else
            %ktkd = txa(ktkd, 5, U+U', kernel{mm}.dim);
            kd = txa(kd, 5, U, kernel{mm}.dim);
            ktd = txa(ktd, 5, U', kernel{mm}.dim);
        end
    elseif strcmpi(kernel{mm}.Operation, 'Pointwise')
        
        U = applyCoords_Multiply(kernel{mm}.Data, kernel{mm}.dim, ...
            xs, ys, zs, ts);
        %kd = bsxfun(@times, shape(U, size(kd), filters{mm}.dim), kd);
        %ktd = bsxfun(@times, shape(U, size(ktd), filters{mm}.dim), ktd);
        kd = bsxfun(@times, U, kd);
        ktd = bsxfun(@times, U, ktd);
    end
    
end

% ktkd ?? (K' + K) * d
ktkd = kd + ktd;

% And now the last part: sum it up...
% This used to be f = colVec(d)'*kernel*colVec(d)
colVec = @(x) reshape(x, [], 1);
f = colVec(d)'*colVec(kd);

% Evaluate Df.
% 
% Df = Dx'*U'*K*U*x + x'*U'*K*U*Dx
%
% Presently, d = U*x.  I will let Dx be the identity operator here.  So

dd = ktkd;

%dd = txa(d, 5, kernel, 4, 1); % kernel*colVec(d);
%fprintf('dd is this big: \n');
%disp(size(dd));

for mm = numel(filters):-1:1
    
    %fprintf('Filter:\n');
    %disp(filters{mm});
    
    %if isa(filters{mm}.Data, 'function_handle')
    %    U = applyCoordinates(filters{mm}.Data, filters{mm}.dim, ...
    %        xs, ys, zs, ts);
    %else
    %    U = filters{mm}.Data;
    %end
    
    if any(strcmpi(filters{mm}.Operation, {'Matrix', 'MatrixArray'})) 
        
        U = applyCoordinates(filters{mm}.Data, filters{mm}.dim, ...
            xs, ys, zs, ts);
        
        if iscell(U)        % Cell array of matrices
            
            % Matlab's cellfun() does not support cell arrays of sparse
            % matrices.  WTF Mathworks, ??????.
            txU = cell(size(U));
            for nn = 1:length(U)
                txU{nn} = U{nn}';
            end
            
            dd = txca(dd, 5, 4, txU, filters{mm}.dim, 2);
        else                % Multiply everything by one matrix
            dd = txa(dd, 5, U', filters{mm}.dim);
        end
    
    elseif strcmpi(filters{mm}.Operation, 'Pointwise')
        
        U = applyCoords_Multiply(filters{mm}.Data, filters{mm}.dim, ...
            xs, ys, zs, ts);
        %dd = bsxfun(@times, shape(conj(U), size(dd), filters{mm}.dim), dd);
        dd = bsxfun(@times, conj(U), dd);
        
    end
    
end

% Now dd = U'*(K'+K)*U*x.  Perfect.  Right?

Df = dd;

% Now dd = U'*K*U*x.  Twice its real part is the sensitivity, I guess.
% Df = 2*real(dd);


function y = shape(x, sz, dims)
% make x 

allCoords = 1:numel(sz);

ySize = sz;
ySize(setdiff(allCoords, dims)) = 1;

xSize = size(x);


ySize(dims) = size(x);

y = reshape(x, ySize);


