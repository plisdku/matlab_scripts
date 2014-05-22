function [x, fval, iter, xHist, fHist, DfHist] = extremize(fn, x0, varargin)
% extremize    minimize or maximize a function
%
% [x, fval, iter] = extremize(fun, x0, named parameters)
%
% Minimize (default) or maximize a function, and try really hard to succeed
% even when the function is choppy as hell.  Lots of juicy,
% hard-to-understand parameters.
%
% Full output:
%
% [x fval iter xHist fHist DfHist] = extremize(fn, x0, varargin)
% 
%
% Named parameters:
%
% Direction     'Minimize' or 'Maximize'
% Bounds        N x 2 array [minX, maxX] (search bounds)
% Deltas        [smallDelta bigDelta], relative change considered small/big
% GrowthFactors [decel accel], step size factors
% GoalF         value of F at which extremize will terminate
% MaxIter       maximum number of iterations
% MaxStep       Nx1 or scalar, maximum step size
% Silent        true or false, governs what gets written to stdout
% TolF          absolute change in F for termination
% RelTolF       relative change in F for termination
% Callback      a function f(xHistory, fHistory, DfHistory)
%

    X.Direction = 'Minimize';
    X.Bounds = [];
    X.Deltas = [];
    X.GrowthFactors = [];
    X.GoalF = [];
    X.MaxIter = 100;
    X.MaxStep = [];
    X.Silent = false;
    X.TolF = 1e-10;
    X.RelTolF = 1e-6;
    X.Callback = @(xHist, fHist, DfHist) pause(0);

    X = parseargs(X, varargin{:});
    
    validate_x0_and_Bounds(x0, X.Bounds);
    
    setDefaults(x0);
    validateInput(x0, X);
    
    if strcmpi(X.Direction, 'Minimize')
        fSign = 1;
    else
        fSign = -1;
    end
    
    % Extract some parameters from X; give them intuitive names.
    bigDelta = X.Deltas(2);
    smallDelta = X.Deltas(1);
    shrinkFactor = X.GrowthFactors(1);
    growFactor = X.GrowthFactors(2);
    
    x = x0;
    done = 0;
    fHist = [];
    DfHist = [];
    xHist = [];
    xBest = x0;
    fBest = Inf;
    stepSizes = 0.1*X.MaxStep;
    
    iter = 0;
    while ~done
        report('# %i\n', iter);
        
        iter = iter + 1;
        
        % EVALUATE THE FUNCTION
        % ... and flip the sign if we're maximizing instead of minimizing.
        [f_eval, Df_eval] = fn(x);
        
        fHist(end+1) = f_eval;
        DfHist(end+1,:) = Df_eval;
        xHist(:,end+1) = x;
        
        f = fSign*f_eval;
        Df = fSign*Df_eval;
        
        X.Callback(xHist, fHist, DfHist);
        
        % Check if we've reached the goal value.
        if (fSign < 0 && f_eval > X.GoalF) || ...
            (fSign > 0 && f_eval < X.GoalF)
            
            report('Goal value of F reached.\n');
            fval = f_eval;
            return;
        end
        
        if iter >= X.MaxIter
            report('Maximum number of iterations reached.\n');
            fval = f_eval;
            return;
        end
        
        % Check if change in F is below TolF in last N iterations
        N_tolF = 3; % number of small changes in F that we require.
        if iter >= N_tolF + 1
            if all(abs(diff(fHist(end-N_tolF:end))) < X.TolF)
                report('%i changes in f smaller than TolF = %g\n', ...
                    N_tolF, X.TolF);
                fval = f_eval;
                return;
            end
            
            if all(abs(diff(fHist(end-N_tolF:end))/f) < X.RelTolF)
                report('%i relative changes in f smaller than RelTolF = %g\n', ...
                    N_tolF, X.RelTolF);
                fval = f_eval;
                return
            end;
        end
        
        
        activeConstraints = (Df' < 0 & x == X.Bounds(:,1)) | ...
            (Df' > 0 & x == X.Bounds(:,2));
        ii = ~activeConstraints';
    
        if any(activeConstraints)
            report('\tconstraints: ');
            report('%i ', find(activeConstraints));
            report('\n');
        end
    
        if f < fBest
            xBest = x;
            fBest = f;
        end
    
        if iter >= 2
            isBetter = (fSign*fHist(end) < fSign*fHist(end-1));
        
            if ~isBetter % step halfway across and average the gradients.  uh yeah.
                stepSizes = stepSizes * shrinkFactor;
                x = 0.5*(xBest + x);
                %DfHist(end,:) = DfHist(end-1,:);
            
                Df = 0.5*(fSign*DfHist(end-1,:) + Df);
            
                report('\tregressed and averaging\n');
            
            elseif iter >= 3
            
                relChange = @(A, B) abs(A-B) ./ abs(B);
            
                bigChanges = (relChange(DfHist(end,:), DfHist(end-1,:)) > bigDelta) ...
                    & (relChange(DfHist(end-1,:), DfHist(end-2,:)) > bigDelta);
            
                smallChanges = (relChange(DfHist(end,:), DfHist(end-1,:)) < smallDelta) ...
                    & (relChange(DfHist(end-1,:), DfHist(end-2,:)) < smallDelta);

                if any(bigChanges)
                    report('\tbig changes: ');
                    report('%i ', find(bigChanges));
                    report('\n');
                end

                if any(smallChanges)
                    report('\tsmall changes: ');
                    report('%i ', find(smallChanges));
                    report('\n');
                end

                stepSizes(bigChanges & ii) = stepSizes(bigChanges & ii) * shrinkFactor;
                stepSizes(smallChanges & ii) = stepSizes(smallChanges & ii) * growFactor;
                
                isBetter = all(diff(fHist(end-2:end)) < 0);
            
                if isBetter
                    stepSizes = stepSizes * growFactor;
                end
            end
        end
    
        iStepTooBig = stepSizes > X.MaxStep;
        stepSizes(iStepTooBig) = X.MaxStep(iStepTooBig);
    
        unitGrad = Df / norm(Df(ii));
        unitGrad(isnan(unitGrad) | isinf(unitGrad)) = 0;
    
        step = -stepSizes .* unitGrad';
        step(activeConstraints) = 0;
        
        xPrev = x;
        x(:) = x(:) + step(:);
    
        iBelowMin = x < X.Bounds(:,1);
        iAboveMax = x > X.Bounds(:,2);
        x(iBelowMin) = X.Bounds(iBelowMin,1);
        x(iAboveMax) = X.Bounds(iAboveMax,2);
    end
    
    
    
    
    
    
    
    
    return;
    
    % Nested functions below.
    
    
    
    
    
    function report(varargin)
        if ~X.Silent
            fprintf(varargin{:});
        end
    end
    
    % Fill in options that the user did not set
    function setDefaults(x0)
        if isempty(X.Bounds)
            X.Bounds = [-Inf(numel(x0),1), Inf(numel(x0),1)];
        end
        
        if isempty(X.Deltas)
            X.Deltas = [1 1];
        end
        
        if isempty(X.GrowthFactors)
            X.GrowthFactors = [0.3 1.8];
        end
        
        if isempty(X.GoalF)
            if strcmpi(X.Direction, 'Minimize')
                X.GoalF = -realmax;
            else
                X.GoalF = realmax;
            end
        end
        
        if isempty(X.MaxStep)
            X.MaxStep = 0.2*(X.Bounds(:,2) - X.Bounds(:,1));
        end
        
        % Adjust max step if needed
        if numel(X.MaxStep) == 1 && numel(x0) > 1
            X.MaxStep = repmat(X.MaxStep, size(x0));
        end
            
    end
    

end


function validate_x0_and_Bounds(x0, bounds)
    if ~iscolumn(x0)
        error('x0 must be a column vector.  Sorry for the inconvenience, have a nice day!');
    end
    
    if isempty(bounds)
        error('Please provide ''Bounds'' to search in, I''m not clairvoyant!');
    end
    
    if ~isequal(size(bounds), size(x0).*[1 2])
        error('Bounds must be [minX maxX] where minX, maxX are size of x0');
    end
end

function validateInput(x0, X)
    
    if ~strcmpi(X.Direction, 'Minimize') && ~strcmpi(X.Direction, 'Maximize')
        error('Direction must be Maximize or Minimize');
    end
    
    iBad = find(X.Bounds(:,2) < X.Bounds(:,1));
    if ~isempty(iBad)
        error('Some max bounds are less than min bounds (first one is #%i)',...
            iBad(1));
    end
    
    iBad = find(x0 < X.Bounds(:,1) | x0 > X.Bounds(:,2));
    if ~isempty(iBad)
        error('x0 is out of bounds (first error at element #%i)', ...
            iBad(1));
    end
    
    if numel(X.GrowthFactors) ~= 2
        error('Growth factors must be a two-element array [shrink grow].');
    end
    
    if ~isscalar(X.GoalF)
        error('GoalF must be a scalar.');
    end
    
    if ~isscalar(X.MaxIter)
        error('MaxIter must be a scalar.');
    end
    
    if ~isequal(size(X.MaxStep), size(x0))
        error('MaxStep should be a scalar or the same size as x0.');
    end
    
    if X.GrowthFactors(1) > 1
        error('First growth factor must be < 1');
    end
    
    if X.GrowthFactors(2) < 1
        error('Second growth factor must be > 1.');
    end
    
    if ~isscalar(X.TolF) || X.TolF < 0
        error('TolF must be a nonnegative scalar.');
    end
    
    if ~isscalar(X.RelTolF) || X.RelTolF < 0
        error('RelTolF must be a nonnegative scalar.');
    end
end