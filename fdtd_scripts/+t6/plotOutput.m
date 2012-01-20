function plotOutput(fileName, varargin)
% plotOutput(fileName, period)
%   Plot 1D and 2D output data from FDTD program.  1D data will be
%   plotted vs time, and 2D output will be presented frame by frame
%   without pause.  Use getOutputFrames or openOutputFile and
%   getOutputFrame to work with the data.
%

import t6.*

X.Period = 1;
X.YLim = [-1 1];
X.CLim = [];
X.UnitString = [];
X = parseargs(X, varargin{:});

period = X.Period;

imagesc_args = {};
if ~isempty(X.CLim)
    imagesc_args = {X.CLim};
end

file = OutputFile(fileName);

fieldNames = {};
for nn = 1:file.numFields()
    fieldNames{nn} = file.Fields{nn}.Name;
end

if file.numRegions() > 1
    error('Function does not yet work with multi-region outputs.');
end

if file.numFields() > 1
    xyzPos = file.positions('Field', 1);
    dim = [cellfun(@numel, xyzPos), file.numFields()];
    %dim = [file.Regions.Size(1,:), file.numFields()];
else
    %dim = file.Regions.Size;
    xyzPos = file.positions();
    dim = cellfun(@numel, xyzPos);
end

numFrames = file.numFramesAvailable;

% Dimension cases
%
% 1 1 1         One field component at a point: time trace
% 1 1 1 3       Three components at a point: threefold time trace
% 1 1 N, perms  1D output: plot 1D movie
% 1 1 N 3       1D output: plot threefold 1D movie
% M N 1, perms  2D output: plot 2D movie
% M N 1 3       2D output: plot 2D movie with three parts

if (prod(dim(1:3)) == 1)   % All 0D cases
    file.open();
    data = file.readFrames();
    file.close();
    if (length(dim) == 3)
        plot(squeeze(data));
    elseif (length(dim) == 4)  % assuming dim4 = 3
        plot(squeeze(data)');
        legend(fieldNames{:})
    end
elseif (nnz(dim(1:3) == 1) == 2)   % All 1D cases
    if (length(dim) == 3)
        frameNum = 1;
        file.open
        data = file.readFrames('NumFrames', 1);
        while frameNum <= numFrames
            if (mod(frameNum, period) == 0)
                plot(squeeze(data));
                ylim(X.YLim);
                legend(fieldNames{:})
                title(sprintf('Frame %i', frameNum));
                pause(0.01);
            end
            if frameNum < numFrames
                data = file.readFrames('NumFrames', 1);
            end
            frameNum = frameNum + 1;
        end
        file.close
    elseif (length(dim) == 4)
        frameNum = 1;
        file.open
        data = file.readFrames('NumFrames', 1);
        while frameNum <= numFrames
            if (mod(frameNum, period) == 0)
                plot(squeeze(data(:,:,:,:)));
                ylim(X.YLim);
                legend(fieldNames{:})
                title(sprintf('Frame %i', frameNum));
                pause(0.01);
            end
            if frameNum < numFrames
                data = file.readFrames('NumFrames', 1);
            end
            frameNum = frameNum + 1;
        end
        file.close
    end
elseif (nnz(dim(1:3) == 1) == 1)    % All 2D cases
    
    if file.Regions.Size(1,1) == 1
        row = 2; col = 3;
    elseif file.Regions.Size(1,2) == 1
        row = 1; col = 3;
    else
        row = 1; col = 2;
    end
    
    if file.numFields() == 1
        frameNum = 1;
        file.open();
        data = file.readFrames('NumFrames', 1);
        
        xyzPos = file.positions(); % size {3}
        
        while frameNum <= numFrames
            if (mod(frameNum, period) == 0)
                imagesc_centered(xyzPos{row}, xyzPos{col}, ...
                    transpose(squeeze(data)), imagesc_args{:});
                axis image
                set(gca, 'YDir', 'Normal');
                colorbar;
                title(sprintf('Frame %i', frameNum));
                xlabel(sprintf('%s%s', char('w'+row), X.UnitString));
                ylabel(sprintf('%s%s', char('w'+col), X.UnitString));
                pause(0.01);
            end
            if frameNum < numFrames
                data = file.readFrames('NumFrames', 1);
            end
            frameNum = frameNum + 1;
        end
    elseif file.numFields() == 3
        frameNum = 1;
        file.open
        data = file.readFrames('NumFrames', 1);
        
        
        while frameNum <= numFrames
            if (mod(frameNum, period) == 0)
                datsize = size(squeeze(data(:,:,:,1)));
                clf
                if diff(xyzPos{row}{2}([1 end])) < ...
                    diff(xyzPos{col}{2}([1 end]))
                    subplot(131);
                    imagesc_centered(xyzPos{row}{1}, xyzPos{col}{1}, ...
                        transpose(squeeze(data(:,:,:,1))));
                    axis image
                    set(gca, 'YDir', 'Normal');
                    colorbar;
                    title(file.Fields{1}.Name);
                    xlabel(sprintf('%s%s', char('w'+row), X.UnitString));
                    ylabel(sprintf('%s%s', char('w'+col), X.UnitString));
                    
                    subplot(132)
                    imagesc_centered(xyzPos{row}{2}, xyzPos{col}{2}, ...
                        transpose(squeeze(data(:,:,:,2))));
                    axis image
                    set(gca, 'YDir', 'Normal');
                    colorbar;
                    title(file.Fields{2}.Name);
                    xlabel(sprintf('%s%s', char('w'+row), X.UnitString));
                    ylabel(sprintf('%s%s', char('w'+col), X.UnitString));
                    
                    subplot(133)
                    imagesc_centered(xyzPos{row}{3}, xyzPos{col}{3}, ...
                        transpose(squeeze(data(:,:,:,3))));
                    axis image
                    set(gca, 'YDir', 'Normal');
                    colorbar;
                    title(file.Fields{3}.Name);
                    xlabel(sprintf('%s%s', char('w'+row), X.UnitString));
                    ylabel(sprintf('%s%s', char('w'+col), X.UnitString));
                else
                    subplot(311)
                    imagesc_centered(xyzPos{row}{1}, xyzPos{col}{1}, ...
                        transpose(squeeze(data(:,:,:,1))));
                    axis image
                    set(gca, 'YDir', 'Normal');
                    colorbar;
                    title(file.Fields{1}.Name);
                    xlabel(sprintf('%s%s', char('w'+row), X.UnitString));
                    ylabel(sprintf('%s%s', char('w'+col), X.UnitString));
                    
                    subplot(312)
                    imagesc_centered(xyzPos{row}{2}, xyzPos{col}{2}, ...
                        transpose(squeeze(data(:,:,:,2))));
                    axis image
                    set(gca, 'YDir', 'Normal');
                    colorbar;
                    title(file.Fields{2}.Name);
                    xlabel(sprintf('%s%s', char('w'+row), X.UnitString));
                    ylabel(sprintf('%s%s', char('w'+col), X.UnitString));
                    
                    subplot(313)
                    imagesc_centered(xyzPos{row}{3}, xyzPos{col}{3}, ...
                        transpose(squeeze(data(:,:,:,3))));
                    axis image
                    set(gca, 'YDir', 'Normal');
                    colorbar;
                    title(file.Fields{3}.Name);
                    xlabel(sprintf('%s%s', char('w'+row), X.UnitString));
                    ylabel(sprintf('%s%s', char('w'+col), X.UnitString));
                end
                pause(0.01);
            end
            if frameNum < numFrames
                data = file.readFrames('NumFrames', 1);
            end
            frameNum = frameNum + 1;
        end
        file.close
    end
elseif (nnz(dim(1:3) == 1) == 0)   % All 3D cases
    error('Cannot plot 3D data with plotOutput (yet).');
end

