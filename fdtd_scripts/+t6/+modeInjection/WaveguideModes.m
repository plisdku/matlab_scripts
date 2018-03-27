classdef WaveguideModes
% WaveguideModes  Class handling 1D and 2D waveguide mode sourcing
%
% Usage: WaveguideModes(dimensionString, named parameters)
%
% dimensionString is either '1D' or '2D'
%
% 1D named parameters:
%   Wavelength
%   IndexBounds         [minIndex, maxIndex]: range of mode indices to search
%   Permittivities      length-N array of relative permittivities
%   Permeabilities      length-N array of relative permeabilities
%   Mode                'tey' or 'tmy'
%   BoundariesX         length-N-1 array of boundary positions
%   X                   length-M ascending array of x positions
%
%
% 2D named parameters:
%   Wavelength
%   GuessIndex          Approximate mode index to search near
%   NumModes            How many modes to solve for
%   X                   length-M ascending array of x positions
%   Y                   length-N ascending array of y positions
%   Permittivity        size MxN array of relative permittivities
%   BoundaryConditions  a string, e.g. '0000', compatible with modesolver package
%
    
    properties
        E
        H
        lambda
        xE
        yE
        xH
        yH
        nEff
    end
    
    methods
        function obj = WaveguideModes(dimensionString, varargin)
            
            if strcmpi(dimensionString, '2D')        
                obj = obj.init2D(varargin{:});
            elseif strcmpi(dimensionString, '1D')
                obj = obj.init1D(varargin{:});
            else
                error('Please specify 1D or 2D');
            end
            
        end
        
        function obj = init1D(obj, varargin)
            X.Wavelength = [];
            X.IndexBounds = [];
            X.Permittivities = [];
            X.Permeabilities = [];
            X.Mode = {'tmy', 'tey'};
            X.BoundariesX = [];
            X.X = [];
            X = parseargs(X, varargin{:});
            
            k0 = 2*pi/X.Wavelength;
            omega = k0;
            
            if strcmpi(X.Mode, 'tmy')
                mode = 'tm';
            else
                mode = 'te';
            end
            
            k = tmm.modes(X.BoundariesX, X.Permittivities, X.Permeabilities, ...
                omega, k0*X.IndexBounds(1), k0*X.IndexBounds(2), mode);
            
            obj.nEff = k/k0;
            obj.lambda = X.Wavelength;
            obj.xE = X.X;
            obj.xH = X.X;
            obj.yE = 0;
            obj.yH = 0;
            obj.E = zeros(numel(X.X), 1, 3, numel(k));
            obj.H = obj.E;
                        
            % We can figure out whether or not a mode is bound of course:
            indices = real(sqrt(X.Permittivities.*X.Permeabilities));
            isBound = @(modeNum) any(obj.nEff > indices(1)) | ...
                any(obj.nEff) > indices(end);
            
            if strcmpi(mode, 'te')
                for mm = 1:numel(k)
                    [ex, hy, hz] = tmm.solveTE(X.BoundariesX, X.Permittivities,...
                        X.Permeabilities, omega, k(mm), X.X, isBound(mm));
                    
                    % Forward permutation: ex becomes ey and so on.
                    
                    obj.E(:,1,2,mm) = ex;
                    obj.H(:,1,3,mm) = hy;
                    obj.H(:,1,1,mm) = hz;
                end
            elseif strcmpi(mode, 'tm')
                for mm = 1:numel(k)
                    [hx, ey, ez] = tmm.solveTM(X.BoundariesX, X.Permittivities,...
                        X.Permeabilities, omega, k(mm), X.X, isBound(mm));
                    
                    obj.H(:,:,2,mm) = hx; % forward permute to our coordinates.
                    obj.E(:,:,3,mm) = ey;
                    obj.E(:,:,1,mm) = ez;
                end
            else
                error('what?');
            end
            
        end
        
        function obj = init2D(obj, varargin)
            X.Wavelength = [];
            X.GuessIndex = [];
            X.NumModes = 1;
            X.X = [];
            X.Y = [];
            X.Permittivity = [];
            X.BoundaryConditions = '0000';
            
            X = parseargs(X, varargin{:});
            
            obj.lambda = X.Wavelength;
            
            rowVec = @(A) reshape(A, 1, []);
            
            obj.xH = [0, rowVec(0.5*(X.X(1:end-1) + X.X(2:end))), 0];
            obj.yH = [0, rowVec(0.5*(X.Y(1:end-1) + X.Y(2:end))), 0];
            obj.xH([1 end]) = 2*obj.xH([2 end-1]) - obj.xH([3 end-2]);
            obj.yH([1 end]) = 2*obj.yH([2 end-1]) - obj.yH([3 end-2]);

            obj.xE = X.X;
            obj.yE = X.Y;
            
            [ex, ey, ez, hx, hy, hz, obj.nEff] = modesolver.fvmodes(X.Wavelength, ...
                X.GuessIndex, X.NumModes, ...
                diff(obj.xH), diff(obj.yH), X.Permittivity, ...
                X.BoundaryConditions);
            
            totalEnergy = obj.complexPower(ex, ey, ez, hx, hy, hz, ...
                obj.xE, obj.yE, obj.xH, obj.yH);
            sqrtEnergy = sqrt(totalEnergy);

            ex = bsxfun(@times, ex, 1./sqrtEnergy);
            ey = bsxfun(@times, ey, 1./sqrtEnergy);
            ez = bsxfun(@times, ez, 1./sqrtEnergy);
            hx = bsxfun(@times, hx, 1./sqrtEnergy);
            hy = bsxfun(@times, hy, 1./sqrtEnergy);
            hz = bsxfun(@times, hz, 1./sqrtEnergy);
            
            newEnergy = obj.complexPower(ex, ey, ez, hx, hy, hz, ...
                obj.xE, obj.yE, obj.xH, obj.yH);
            
            obj.E = permute(cat(4, ex, ey, ez), [1 2 4 3]);
            obj.H = permute(cat(4, hx, hy, hz), [1 2 4 3]);
        end
    
        function totalEnergy = calcEnergy(obj, ex, ey, ez, hx, hy, hz, ...
            xE, yE, xH, yH)

            exhy = ex.*gridInterp(xH, yH, hy, xE, yE);
            eyhx = ey.*gridInterp(xH, yH, hx, xE, yE);
            
            totalEnergy = 0.5*real(trapzn({xE, yE}, exhy - eyhx));
        end
        
        function totalEnergy = complexPower(obj, ex, ey, ez, hx, hy, hz, ...
            xE, yE, xH, yH)

            exhy = ex.*gridInterp(xH, yH, conj(hy), xE, yE);
            eyhx = ey.*gridInterp(xH, yH, conj(hx), xE, yE);
            
            totalEnergy = 0.5*trapzn({xE, yE}, exhy - eyhx);
        end
        
        
        function [Eout, Hout] = interpolate(obj, varargin)
            import t6.modeInjection.*
            X.X = [];
            X.Y = [];
            X.Z = [];
            X.Rotation = eye(3);
            X.Translation = [0 0 0]';
            X = parseargs(X, varargin{:});
            
            szE = size(obj.E);
            szH = size(obj.H);
            
            numModes = size(obj.E, 4);
            
            % We will need to add some fictitious dimensions (one or two) to 
            % make use of Matlab's interpolation functions for rotating the
            % fields.  The z dimension must be duplicated once over, and the
            % y dimension may need this treatment as well.
            
            repeats = [1 1 2 1 1];
            if numel(obj.yE) == 1
                repeats(2) = 2;
            end
            
            E_repeated = repmat(reshape(obj.E, [szE(1:2) 1 szE(3) numModes]), ...
                repeats);
            H_repeated = repmat(reshape(obj.H, [szH(1:2) 1 szH(3) numModes]), ...
                repeats);
            
            % Spatial locations of these duplicated samples: really far away!
            distantPositions = realmax*[-1 1];
            
            Eout = zeros(numel(X.X), numel(X.Y), numel(X.Z), 3, numModes);
            Hout = Eout;
            
            for mm = 1:numModes
                if numel(obj.yE) == 1 % happens when I use transfer matrices
                    Eout(:,:,:,:,mm) = transformField(X.Rotation, X.Translation, ...
                        obj.xE, distantPositions, distantPositions, ...
                        E_repeated(:,:,:,:,mm), X.X, X.Y, X.Z);
                    Hout(:,:,:,:,mm) = transformField(X.Rotation, X.Translation, ...
                        obj.xH, distantPositions, distantPositions, ...
                        H_repeated(:,:,:,:,mm), X.X, X.Y, X.Z);
                else % happens when I use the 2D mode solver
                    Eout(:,:,:,:,mm) = transformField(X.Rotation, X.Translation, ...
                        obj.xE, obj.yE, distantPositions, ...
                        E_repeated(:,:,:,:,mm), X.X, X.Y, X.Z);
                    Hout(:,:,:,:,mm) = transformField(X.Rotation, X.Translation, ...
                        obj.xH, obj.yH, distantPositions, ...
                        H_repeated(:,:,:,:,mm), X.X, X.Y, X.Z);
                end
            end
            
        end
        
        
    end % methods
    
end