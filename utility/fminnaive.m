function [x fval iter] = fminnaive(fn, x0, minX, maxX, bigDelta, smallDelta, ...
    shrinkFactor, growFactor, goalF)
% [x fval iter] = fminnaive(fn, x0, minX, maxX, bigDelta, smallDelta, ...
%   shrinkFactor, growFactor, goalF)

if nargin < 5
    bigDelta = 1;
    smallDelta = 1;
    shrinkFactor = 0.3;
    growFactor = 1.8;
    goalF = realmax;
end

x = x0;

done = 0;

fHist = [];
DfHist = [];
xPrev = x0;
xBest = x0;
fBest = Inf;

stepSizes = 0.03 * (maxX - minX);

maxStepSizes = 0.2 * (maxX - minX);
%minStepSizes = 0*ones(size(x0));

%shrinkFactor = 0.5;
%growFactor = 1.25;
%bigDelta = 0.7;
%smallDelta = 1.0;

iter = 0;
while ~done
    fprintf('# %i\n', iter);
    iter = iter + 1;
    [f, Df] = fn(x);
    
    if (-f) > goalF || iter > 100
        fval = f;
        return;
    end
    
    fHist(end+1) = f;
    DfHist(end+1,:) = Df;
    
    activeConstraints = (Df' < 0 & x == minX) | ...
        (Df' > 0 & x == maxX);
    ii = ~activeConstraints';
    
    if any(activeConstraints)
        fprintf('\tconstraints: ');
        fprintf('%i ', find(activeConstraints));
        fprintf('\n');
    end
    
    if f < fBest
        xBest = x;
        fBest = f;
    end
    
    if iter >= 2
        isBetter = (fHist(end) < fHist(end-1));
        
        if ~isBetter % step halfway across and average the gradients.  uh yeah.
            stepSizes = stepSizes * shrinkFactor;
            x = 0.5*(xBest + x);
            %DfHist(end,:) = DfHist(end-1,:);
            
            Df = 0.5*(DfHist(end-1,:) + Df);
            
            fprintf('\tregressed and averaging\n');
            
        elseif iter >= 3
            
            relChange = @(A, B) abs(A-B) ./ abs(B);
            
            bigChanges = (relChange(DfHist(end,:), DfHist(end-1,:)) > bigDelta) ...
                & (relChange(DfHist(end-1,:), DfHist(end-2,:)) > bigDelta);
            
            smallChanges = (relChange(DfHist(end,:), DfHist(end-1,:)) < smallDelta) ...
                & (relChange(DfHist(end-1,:), DfHist(end-2,:)) < smallDelta);
            
            %signSwitches = (sign(DfHist(end,:)) == sign(DfHist(end-2,:))) & ...
            %    (sign(DfHist(end,:)) ~= sign(DfHist(end-1,:)));
            %sameSigns = (sign(DfHist(end,:)) == sign(DfHist(end-2,:))) & ...
            %    (sign(DfHist(end,:)) == sign(DfHist(end-1,:)));
            
            if any(bigChanges)
                fprintf('\tbig changes: ');
                fprintf('%i ', find(bigChanges));
                fprintf('\n');
            end
            
            if any(smallChanges)
                fprintf('\tsmall changes: ');
                fprintf('%i ', find(smallChanges));
                fprintf('\n');
            end

            stepSizes(bigChanges & ii) = stepSizes(bigChanges & ii) * shrinkFactor;
            stepSizes(smallChanges & ii) = stepSizes(smallChanges & ii) * growFactor;
            
            
            isBetter = all(diff(fHist(end-2:end)) < 0);
            
            if isBetter
                stepSizes = stepSizes * growFactor;
            end
        end
    end
    
    stepSizes(stepSizes > maxStepSizes) = ...
        maxStepSizes(stepSizes > maxStepSizes);
    
    unitGrad = Df / norm(Df(ii));
    unitGrad(isnan(unitGrad) | isinf(unitGrad)) = 0;
    
    step = -stepSizes .* unitGrad';
    
    xPrev = x;
    x(:) = x(:) + step(:);
    
    x(x < minX) = minX(x < minX);
    x(x > maxX) = maxX(x > maxX);
end
    
    