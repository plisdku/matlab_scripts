function [f Df] = evalQF(data, filters, xs, ys, zs, ts)
% [f Df] = evalQuadraticForm(data, filters, coords)
%

import multiplyTensors.*

% Evaluate f

d = data;

for mm = 1:numel(filters)
    
    %fprintf('Filter:\n');
    %disp(filters{mm});
    
    if isa(filters{mm}.Data, 'function_handle')
        U = applyCoordinates(filters{mm}.Data, filters{mm}.dim, ...
            xs, ys, zs, ts);
    else
        U = filters{mm}.Data;
    end
    
    
    if any(strcmpi(filters{mm}.Operation, {'Matrix', 'MatrixArray'})) 
        
        if iscell(U)        % Cell array of matrices
            d = txca(d, 5, 4, U, filters{mm}.dim, 2);
        else                % Multiply everything by one matrix
            d = txa(d, 5, U, filters{mm}.dim);
        end
    
    elseif strcmpi(filters{mm}.Operation, 'Pointwise')
        
        d = bsxfun(@times, U, d);
        
    end
    
end

% At this point any remaining nonsingular dimensions of d must be eaten by
% a quadratic form, as d' * K * d, to obtain f.

% Now the quadratic form: it should operate over FIELDS to achieve
% something like a Poynting vector optimization!  That's dimension 4.
%kernel = speye([1 1]*size(d, 4));
kernel = 1;

colVec = @(x) reshape(x, [], 1);

f = colVec(d)'*kernel*colVec(d);

%f = txt(conj(d), 5, txa(d, 5, kernel, 4, 1), 5, 4);

%fprintf('f is:\n');
%disp(f);

% and integrate if that's missing:
if numel(f) > 1
    error('Still need an integral.');
end

%fprintf('f = %2.5e\n', f);

% Evaluate Df.
% 
% Df = Dx'*U'*K*U*x + x'*U'*K*U*Dx
%
% Presently, d = U*x.  I will let Dx be the identity operator here.  So

dd = txa(d, 5, kernel, 4, 1); % kernel*colVec(d);
%fprintf('dd is this big: \n');
%disp(size(dd));

for mm = numel(filters):-1:1
    
    %fprintf('Filter:\n');
    %disp(filters{mm});
    
    if isa(filters{mm}.Data, 'function_handle')
        U = applyCoordinates(filters{mm}.Data, filters{mm}.dim, ...
            xs, ys, zs, ts);
    else
        U = filters{mm}.Data;
    end
    
    
    if any(strcmpi(filters{mm}.Operation, {'Matrix', 'MatrixArray'})) 
        
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
        
        dd = bsxfun(@times, conj(U), dd);
        
    end
    
    %{
    if  isfield(filters{mm}, 'field')
        ff = filters{mm}.field;
    else
        ff = [];
    end
    
    if isfield(filters{mm}, 'U')
    
        if isa(filters{mm}.U, 'function_handle')
            U = filters{mm}.U(coordinates{filters{mm}.dim});
        else
            U = filters{mm}.U;
        end
        
        if iscell(U)
            % Matlab's cellfun() does not support cell arrays of sparse
            % matrices.  WTF Mathworks, ??????.
            txU = cell(size(U));
            for nn = 1:length(U)
                txU{nn} = U{nn}';
            end
            
            dd = txca(dd, 5, 4, txU, filters{mm}.dim, 2);
        elseif ~isempty(ff)
            dd(:,:,:,ff,:) = txa(dd(:,:,:,ff,:), 5, U', filters{mm}.dim);
        else
            dd = txa(dd, 5, U', filters{mm}.dim);
        end
    
    elseif isfield(filters{mm}, 'y')
        
        if isa(filters{mm}.y, 'function_handle')
            y = filters{mm}.y(coordinates{filters{mm}.dim});
        else
            y = filters{mm}.y;
        end
        
        if ~isempty(ff)
            dd(:,:,:,ff,:) = bsxfun(@times, conj(y), dd(:,:,:,ff,:));
        else
            dd = bsxfun(@times, conj(y), dd);
        end
    end
    %}
    
end

% Now dd = U'*K*U*x.  Twice its real part is the sensitivity, I guess.

Df = 2*real(dd);

