classdef WaveguideModes
    
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
        function obj = WaveguideModes(varargin)
            
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
            
            [ex ey ez hx hy hz obj.nEff] = modesolver.fvmodes(X.Wavelength, ...
                X.GuessIndex, X.NumModes, ...
                diff(obj.xH), diff(obj.yH), X.Permittivity, ...
                X.BoundaryConditions);

            totalEnergy = obj.calcEnergy(ex, ey, ez, hx, hy, hz, ...
                obj.xE, obj.yE, obj.xH, obj.yH);
            sqrtEnergy = sqrt(totalEnergy);

            ex = bsxfun(@times, ex, 1./sqrtEnergy);
            ey = bsxfun(@times, ey, 1./sqrtEnergy);
            ez = bsxfun(@times, ez, 1./sqrtEnergy);
            hx = bsxfun(@times, hx, 1./sqrtEnergy);
            hy = bsxfun(@times, hy, 1./sqrtEnergy);
            hz = bsxfun(@times, hz, 1./sqrtEnergy);
            
            obj.E = permute(cat(4, ex, ey, ez), [1 2 4 3]);
            obj.H = permute(cat(4, hx, hy, hz), [1 2 4 3]);
        end
    
        function totalEnergy = calcEnergy(obj, ex, ey, ez, hx, hy, hz, ...
            xE, yE, xH, yH)

            exhy = ex.*gridInterp(xH, yH, hy, xE, yE);
            eyhx = ey.*gridInterp(xH, yH, hx, xE, yE);

            totalEnergy = 0.5*real(trapz(yE, trapz(xE, exhy - eyhx, 1), 2));
        end
        
        
        function [Eout Hout] = interpolate(obj, varargin)
            import modeInjection.*
            X.X = [];
            X.Y = [];
            X.Z = [];
            X.Rotation = eye(3);
            X.Translation = [0 0 0]';
            X = parseargs(X, varargin{:});
            
            szE = size(obj.E);
            szH = size(obj.H);
            
            numModes = size(obj.E, 4);
            
            % Put in a fictitious Z dimension
            E_with_z = repmat(reshape(obj.E, [szE(1:2) 1 szE(3) numModes]), ...
                [1 1 2 1 1]);
            H_with_z = repmat(reshape(obj.H, [szH(1:2) 1 szH(3) numModes]), ...
                [1 1 2 1 1]);
            
            zFake = realmax*[-1 1];

            Eout = zeros(numel(X.X), numel(X.Y), numel(X.Z), 3, numModes);
            Hout = Eout;
            
            for mm = 1:numModes
                Eout(:,:,:,:,mm) = transformField(X.Rotation, X.Translation, ...
                    obj.xE, obj.yE, zFake, E_with_z(:,:,:,:,mm), ...
                    X.X, X.Y, X.Z);
                Hout(:,:,:,:,mm) = transformField(X.Rotation, X.Translation, ...
                    obj.xH, obj.yH, zFake, H_with_z(:,:,:,:,mm), ...
                    X.X, X.Y, X.Z);
            end
            
        end
            
            
        
        
    end % methods
    
end