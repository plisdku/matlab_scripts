function [f, Df, f_w, freqs, filteredData] = evalQuadForm(data, pos, timeVals, objFunStruct)
% evalQuadForm  Evaluate objective function and its gradient
%
% [f Df] = evalQuadForm(data, pos, timeVals, objFunStruct)
%
% [f Df] = evalQuadForm(file, objFunStruct)
%
% [f Df] = evalQuadForm(data) provides the sum of squares of the
% field values in data and the corresponding adjoint source Df.
% 
% The order of operations is first to do spatial and temporal filtering,
% then do operations on the fields (e.g. ExH and such).  The spatial
% operations are:
% 
% filterT * filterXYZ * filterZ * filterY * filterX * fields
%
% filterT * filterXYZ * filterZ * filterY * filterX * interpZ * interpY * interpX
% 
% Math background: let the objective function be x'*A*x.  Then the adjoint
% source should be (A' + A)*x.  I think.  I'm trying to do that.

import t6.*

if nargout > 2
    [data_w, freqs] = mySpectrum(data, timeVals);
end

numFields = size(data, 4);
if isempty(objFunStruct.Kernel)
    objFunStruct.Kernel = eye(numFields);
end

for ff = 1:numFields
    oneFieldData = data(:,:,:,ff,:);
    
    % Spatial interpolation
    if ~isempty(objFunStruct.InterpX)
        oneFieldData = multTensor(oneFieldData, objFunStruct.InterpX{ff}, 1);
    end
    
    if ~isempty(objFunStruct.InterpY)
        oneFieldData = multTensor(oneFieldData, objFunStruct.InterpY{ff}, 2);
    end
    
    if ~isempty(objFunStruct.InterpZ)
        oneFieldData = multTensor(oneFieldData, objFunStruct.InterpZ{ff}, 3);
    end
    
    % User integrations
    if ~isempty(objFunStruct.X)
        oneFieldData = multTensor(oneFieldData, ...
            evalFilterMatrix(objFunStruct.X, pos, ff, 1), 1);
    end
    
    if ~isempty(objFunStruct.Y)
        oneFieldData = multTensor(oneFieldData, ...
            evalFilterMatrix(objFunStruct.Y, pos, ff, 2), 2);
    end
    
    if ~isempty(objFunStruct.Z)
        oneFieldData = multTensor(oneFieldData, ...
            evalFilterMatrix(objFunStruct.Z, pos, ff, 3), 3);
    end
    
    if ~isempty(objFunStruct.T)
        oneFieldData = multTensor(oneFieldData, ...
            evalTimeFilter(objFunStruct.T, timeVals{ff}, ff), 5);
        
%         figure(9); clf
%         %image(complex2rgb(transpose(squish(oneFieldData))));
%         imagesc(abs(transpose(squish(oneFieldData))));
%         axis xy image;
%         prettify('Times', 14)
%         title(sprintf('Measured magnitude %i', ff)); colorbar
%         pause
    end
    
    if ~isempty(objFunStruct.XYZ)
        
        filterArray = evalFilterArray(objFunStruct.XYZ, pos, ff);
        
        figure(2); clf
        image(complex2rgb(transpose(squish(filterArray))));
        axis xy image; title(sprintf('Filter %i', ff));
        pause
        
        % I'm adding ones to the end of the size array because reshape
        % will barf if I ever happen to reshape to a
        % less-than-two-dimensional form.
        szDat = [size(oneFieldData) 1 1 1 1];
        
        oneFieldData = reshape(oneFieldData, [prod(szDat(1:3)), szDat(4:end)]);
        filterArray = reshape(filterArray, 1, []);
        
        oneFieldData = multTensor(oneFieldData, filterArray, 1);
        oneFieldData = reshape(oneFieldData, [1 1 1 szDat(4:end)]);
        
        fprintf('Field %i tensor product: %2.2f + %2.2fi \n', ...
            ff, real(oneFieldData), imag(oneFieldData));
    end
    
%    if ~isempty(objFunStruct.T)
%        oneFieldData = multTensor(oneFieldData, ...
%            evalTimeFilter(objFunStruct.T, timeVals{ff}, ff), 5);
%    end
    
    if ff == 1, filteredData = oneFieldData;
    else filteredData(:,:,:,ff,:) = oneFieldData;
    end
end

% Evaluate the objective function itself
% This should be f = x'*A*x.

quadraticAddend = conj(filteredData) .* ...
    multTensor(filteredData, objFunStruct.Kernel, 4);
f = objFunStruct.Factor * sum(quadraticAddend(:));

% Obtain Df
% it should be Df = (A' + A)*x, which we evaluate carefully.
% If K is the kernel, the other operations can be U.  Then
%
% A = U'*K*U
% 
% Df = 2*real(x'*U'*K*U).
%
% The U*x part has been done already (it's data).
% So x'*U' = filteredData' (INCLUDING complex conjugates here, careful!).
%
% y' := x'*U'
% z' := y'*K   (I'll call this variable z' := zz in the code)
% Df = 2*real(z'*U)
%

assert(isequal(objFunStruct.Kernel, objFunStruct.Kernel'));

zz = objFunStruct.Factor * ...
    multTensor(conj(filteredData), objFunStruct.Kernel, 4);

for ff = 1:numFields
    oneFieldData = zz(:,:,:,ff,:);
    
    if ~isempty(objFunStruct.T)
        oneFieldData = multTensor(oneFieldData, ...
            transpose(evalTimeFilter(objFunStruct.T, timeVals{ff}, ff)), 5);
    end
    
    if ~isempty(objFunStruct.XYZ)
        filterArray = evalFilterArray(objFunStruct.XYZ, pos, ff);
        
        % I'm adding a one to the end of the size array because reshape
        % will barf if I ever happen to reshape to a
        % less-than-two-dimensional form.
        szFilter = size(filterArray);
        szFilter(end+1:3) = 1;
        szDat = [size(oneFieldData) 1];
        
        %oneFieldData = reshape(oneFieldData, [prod(szDat(1:3)), szDat(4:end)]);
        filterArray = reshape(filterArray, 1, []);
        
        oneFieldData = multTensor(oneFieldData, transpose(filterArray), 1);
        oneFieldData = reshape(oneFieldData, [szFilter szDat(4:end)]);
        
        % this was never right.
        %oneFieldData = bsxfun(@times, oneFieldData, ...
        %    (evalFilterArray(objFunStruct.XYZ, pos, ff)));
    end
    
    if ~isempty(objFunStruct.Z)
        oneFieldData = multTensor(oneFieldData, ...
            transpose(evalFilterMatrix(objFunStruct.Z, pos, ff, 3)), 3);
    end
    
    if ~isempty(objFunStruct.Y)
        oneFieldData = multTensor(oneFieldData, ...
            transpose(evalFilterMatrix(objFunStruct.Y, pos, ff, 2)), 2);
    end
    
    if ~isempty(objFunStruct.X)
        oneFieldData = multTensor(oneFieldData, ...
            transpose(evalFilterMatrix(objFunStruct.X, pos, ff, 1)), 1);
    end
    
    % Spatial interpolation
    if ~isempty(objFunStruct.InterpZ)
        oneFieldData = multTensor(oneFieldData, ...
            transpose(objFunStruct.InterpZ{ff}), 3);
    end
    
    if ~isempty(objFunStruct.InterpY)
        oneFieldData = multTensor(oneFieldData, ...
            transpose(objFunStruct.InterpY{ff}), 2);
    end
    
    if ~isempty(objFunStruct.InterpX)
        oneFieldData = multTensor(oneFieldData, ...
            transpose(objFunStruct.InterpX{ff}), 1);
    end
    
    if ff == 1; Df = 2*real(oneFieldData);
    else Df(:,:,:,ff,:) = 2*real(oneFieldData);
    end

end

% Obtain f_w if asked-for

if nargout > 2
    for ff = 1:numFields
        oneFieldData_w = data_w(:,:,:,ff,:);
        
        % Spatial interpolation
        if ~isempty(objFunStruct.InterpX)
            oneFieldData_w = multTensor(oneFieldData_w, ...
                objFunStruct.InterpX{ff}, 1);
        end
        
        if ~isempty(objFunStruct.InterpY)
            oneFieldData_w = multTensor(oneFieldData_w, ...
                objFunStruct.InterpY{ff}, 2);
        end
        
        if ~isempty(objFunStruct.InterpZ)
            oneFieldData_w = multTensor(oneFieldData_w, ...
                objFunStruct.InterpZ{ff}, 3);
        end
        
        % User integration
        if ~isempty(objFunStruct.X)
            oneFieldData_w = multTensor(oneFieldData_w, ...
                evalFilterMatrix(objFunStruct.X, pos, ff, 1), 1);
        end

        if ~isempty(objFunStruct.Y)
            oneFieldData_w = multTensor(oneFieldData_w, ...
                evalFilterMatrix(objFunStruct.Y, pos, ff, 2), 2);
        end

        if ~isempty(objFunStruct.Z)
            oneFieldData_w = multTensor(oneFieldData_w, ...
                evalFilterMatrix(objFunStruct.Z, pos, ff, 3), 3);
        end


        if ~isempty(objFunStruct.XYZ)

            filterArray = evalFilterArray(objFunStruct.XYZ, pos, ff);

            % I'm adding a one to the end of the size array because reshape
            % will barf if I ever happen to reshape to a
            % less-than-two-dimensional form.
            szDat = [size(oneFieldData) 1 1 1 1];

            oneFieldData = reshape(oneFieldData, [prod(szDat(1:3)), szDat(4:end)]);
            filterArray = reshape(filterArray, 1, []);

            oneFieldData = multTensor(oneFieldData, filterArray, 1);
            oneFieldData = reshape(oneFieldData, [1 1 1 szDat(4:end)]);

            %oneFieldData = sum(sum(sum(...
            %    bsxfun(@times, oneFieldData, ...
            %    evalFilterArray(objFunStruct.XYZ, pos, ff)), 3), 2), 1);
        end

        % never correct, this...
%        if ~isempty(objFunStruct.XYZ)
%            oneFieldData_w = bsxfun(@times, oneFieldData_w, ...
%                evalFilterArray(objFunStruct.XYZ, pos, ff));
%        end
        
        if ~isempty(objFunStruct.T)
            warning('Cannot do a time filter in freq domain yet');
        end

        if ff == 1, filteredData_w = oneFieldData_w;
        else filteredData_w(:,:,:,ff,:) = oneFieldData_w;
        end
    end

    f_w = objFunStruct.Factor * ...
        squeeze(sum(sum(sum(sum(conj(filteredData_w) .* ...
        multTensor(filteredData_w, objFunStruct.Kernel, 4), ...
        4), 3), 2), 1));
end

end


function vals = evalFilterArray(filter, pos, whichField)
% filter is a cell array of function handles, or a function handle

[xx yy zz] = ndgrid(pos{1}, pos{2}, pos{3});

dxyz = [1 1 1];
for xyz = 1:3
    if length(pos{xyz}) > 1
        dxyz(xyz) = pos{xyz}(2) - pos{xyz}(1);
    end
end

if iscell(filter)
    vals = filter{whichField}(xx, yy, zz) * prod(dxyz);
else
    vals = filter(xx, yy, zz) * prod(dxyz);
end

end


function vals = evalFilterMatrix(filter, pos, whichField, xyz)

if length(pos{xyz}) > 1
    dx = pos{xyz}(2) - pos{xyz}(1);
else
    dx = 1;
end

if iscell(filter)
    if isa(filter{whichField}, 'function_handle')
        vals = filter{whichField}(pos{xyz})*dx;
        
        if length(vals) > 1
            vals([1 end]) = vals([1 end])/2; % implement trapezoidal integration
        end
    else
        vals = filter{whichField};
    end
else
    if isa(filter, 'function_handle')
        vals = filter(pos{xyz})*dx;
        
        if length(vals) > 1
            vals([1 end]) = vals([1 end])/2; % implement trapezoidal integration
        end
    else
        vals = filter;
    end
end

end

function vals = evalTimeFilter(filter, timeVals, whichField)

if length(timeVals) > 2 % I use 2 because of frequency filter?what a kludge!
    dt = timeVals(2) - timeVals(1);
else
    dt = 1;
end

isVector = @(v) ndims(v) == 2 && any(size(v) == 1);

if iscell(filter)
    if isa(filter{whichField}, 'function_handle')
        vals = filter{whichField}(timeVals);
        
        if isVector(vals)
            vals = vals*dt;
            
            if length(vals) > 2
                vals([1 end]) = vals([1 end])/2; % implement trapezoidal integration
            end
        else
            % vals must be a filter matrix
        end
        
    else
        vals = filter{whichField};
    end
else
    if isa(filter, 'function_handle')
        vals = filter(timeVals);
        
        if isVector(vals)
            vals = vals*dt;
            
            if length(vals) > 2
                vals([1 end]) = vals([1 end])/2; % implement trapezoidal integration
            end
        else
            % vals must be a filter matrix
        end
    else
        vals = filter;
    end
end

end


function [data_w, freqs] = mySpectrum(data, timeVals)

data_w = data;
nFields = size(data, 4);

for ff = 1:nFields
    
    if iscell(timeVals)
        [data_w(:,:,:,ff,:), freqs] = t6.analysis.spectrum(...
            data(:,:,:,ff,:), 'Time', timeVals{ff});
    else
        [data_w(:,:,:,ff,:), freqs] = t6.analysis.spectrum(...
            data(:,:,:,ff,:), 'Time', timeVals);
    end
end

end

%% Make sparse diagonal matrix

function B = mySparseDiag(d)

validateattributes(d, {'numeric'}, {'vector'});
B = sparse(1:length(d), 1:length(d), d(:));

end


%% Multiply Tensors

function B = multTensor(T, A, dim)
% mult(T, A, dim) multiplies tensor T by matrix A over dimension dim
% (A*T)

%if isscalar(T) && dim ~= 1
%    warning('Tensor T is a scalar so multiplication is being shifted to first dimension');
%    dim = 1;
%end

if isscalar(A)
    B = A*T;
    return
end

if dim < 1
    error('dim = 0');
elseif dim > ndims(T)
    % Singleton dimensions get stripped off the end of arrays in Matlab.
    % That's the only good reason why this sort of problem should crop up.
    % Fix by adding a singleton dimension to the front of T, doing the A*T
    % operation, and then pushing that singleton dimension to the back of
    % T.  We can't expect to give T any dimensions beyond dim.  Stupid
    % Matlab.
    tPerm = reshape(T, [1 size(T)]);
    
    szPermuted = [1 size(T)];
    
    prodTensor = reshape(A*tPerm(:,:), [size(A,1), size(T)]);
    B = ipermute(prodTensor, [dim 1:dim-1]);
else
    % This is the usual case.  Put desired dimension first.
    permutedDims = [dim, 1:dim-1, dim+1:ndims(T)];
    sz = size(T);
    szPermuted = sz(permutedDims);
    tPerm = permute(T, permutedDims); % Put the desired dimension first
    
    prodTensor = reshape(A*tPerm(:,:), [size(A,1), szPermuted(2:end)]);
    B = ipermute(prodTensor, permutedDims);
end

%if dim < 1 || dim > ndims(T)
%    if size(A, 2) > 1
%        error('dim must be valid dimension of T');
%    else
%        %dim = 1;
%        dim 
%    end
%end

end



