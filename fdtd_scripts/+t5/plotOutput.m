function plotOutput(fileName, period)
% plotOutput(fileName, period)
%   Plot 1D and 2D output data from FDTD program.  1D data will be
%   plotted vs time, and 2D output will be presented frame by frame
%   without pause.  Use getOutputFrames or openOutputFile and
%   getOutputFrame to work with the data.
%
%   Provided for back-compatibility with Trogdor 4.

import t5.*

if (nargin < 2)
    period = 1;
end

file = OutputFile(fileName);

fieldNames = {};
for nn = 1:length(file.Fields)
    fieldNames{nn} = file.Fields{nn}.Name;
end

if length(file.Regions) > 1
    error('Function does not yet work with multi-region outputs.');
end

if length(file.Fields) > 1
    dim = [file.Regions{1}.Size, length(file.Fields)];
else
    dim = file.Regions{1}.Size;
end

numFrames = file.numFramesAvailable;

% Dimension cases
%
% 1 1 1         One field componen at a point: time trace
% 1 1 1 3       Three components at a point: threefold time trace
% 1 1 N, perms  1D output: plot 1D movie
% 1 1 N 3       1D output: plot threefold 1D movie
% M N 1, perms  2D output: plot 2D movie
% M N 1 3       2D output: plot 2D movie with three parts

if (prod(dim(1:3)) == 1)   % All 0D cases
    data = file.read;
    if (length(dim) == 3)
        plot(squeeze(data));
    elseif (length(dim) == 4)  % assuming dim4 = 3
        plot(squeeze(data)');
        legend(fieldNames{:})
        %{
        plot(squeeze(data(1,1,1,1,:)), 'r');
        hold on
        plot(squeeze(data(1,1,1,2,:)), 'g');
        plot(squeeze(data(1,1,1,3,:)), 'b');
        hold off
        legend('Component 1', 'Component 2', 'Component 3');
        %}
    end
elseif (nnz(dim(1:3) == 1) == 2)   % All 1D cases
    if (length(dim) == 3)
        frameNum = 1;
        file.open
        data = file.readFrames(1);
        while frameNum <= numFrames
            if (mod(frameNum, period) == 0)
                plot(squeeze(data));
                ylim([-2, 2]);
                legend(fieldNames{:})
                title(sprintf('Frame %i', frameNum));
                pause(0.01);
            end
            if frameNum < numFrames
                data = file.readFrames(1);
            end
            frameNum = frameNum + 1;
        end
        file.close
    elseif (length(dim) == 4)
        frameNum = 1;
        file.open
        data = file.readFrames(1);
        while frameNum <= numFrames
            if (mod(frameNum, period) == 0)
                plot(squeeze(data(:,:,:,:)));
                %{
                hold on
                plot(squeeze(data(:,:,:,2)), 'g');
                plot(squeeze(data(:,:,:,3)), 'b');
                hold off
                %}
                ylim([-2, 2]);
                title(sprintf('Frame %i', frameNum));
                pause(0.01);
                %pause  
            end
            if frameNum < numFrames
                data = file.readFrames(1);
            end
            frameNum = frameNum + 1;
        end
        file.close
    end
elseif (nnz(dim(1:3) == 1) == 1)    % All 2D cases
    
    coords = 'xyz';
    [xPos yPos zPos] = file.positions;
    xyzPos = {xPos yPos zPos};
    
    if file.Regions{1}.Size(1) == 1
        row = 2; col = 3;
    elseif file.Regions{1}.Size(2) == 1
        row = 1; col = 3;
    else
        row = 1; col = 2;
    end
    
    if (length(dim) == 3)
        frameNum = 1;
        file.open
        data = file.readFrames(1);
        while frameNum <= numFrames
            if (mod(frameNum, period) == 0)
                imagesc_centered(xyzPos{row}, xyzPos{col}, ...
                    transpose(squeeze(data)));
                axis image
                set(gca, 'YDir', 'Normal');
                colorbar;
                title(sprintf('Frame %i', frameNum));
                xlabel(sprintf('%s (m)', coords(row)));
                ylabel(sprintf('%s (m)', coords(col)));
                pause(0.01);
            end
            if frameNum < numFrames
                data = file.readFrames(1);
            end
            frameNum = frameNum + 1;
        end
    elseif (length(dim) == 4)
        frameNum = 1;
        file.open
        data = file.readFrames(1);
        
        
        while frameNum <= numFrames
            if (mod(frameNum, period) == 0)
                datsize = size(squeeze(data(:,:,:,1)));
                clf
                if datsize(1) < datsize(2)   % tall skinny
                    subplot(131);
                    imagesc_centered(xyzPos{row}{1}, xyzPos{col}{1}, ...
                        transpose(squeeze(data(:,:,:,1))));
                    axis image
                    set(gca, 'YDir', 'Normal');
                    colorbar;
                    title(file.Fields{1}.Name);
                    xlabel(sprintf('%s (m)', coords(row)));
                    ylabel(sprintf('%s (m)', coords(col)));
                    
                    subplot(132)
                    imagesc_centered(xyzPos{row}{2}, xyzPos{col}{2}, ...
                        transpose(squeeze(data(:,:,:,2))));
                    axis image
                    colorbar;
                    title(file.Fields{2}.Name);
                    xlabel(sprintf('%s (m)', coords(row)));
                    ylabel(sprintf('%s (m)', coords(col)));
                    
                    subplot(133)
                    imagesc_centered(xyzPos{row}{3}, xyzPos{col}{3}, ...
                        transpose(squeeze(data(:,:,:,3))));
                    axis image
                    colorbar;
                    title(file.Fields{3}.Name);
                    xlabel(sprintf('%s (m)', coords(row)));
                    ylabel(sprintf('%s (m)', coords(col)));
                    %suptitle(sprintf('Frame %i', frameNum));
                else
                    subplot(311)
                    imagesc_centered(xyzPos{row}{1}, xyzPos{col}{1}, ...
                        transpose(squeeze(data(:,:,:,1))));
                    axis image
                    colorbar;
                    title(file.Fields{1}.Name);
                    xlabel(sprintf('%s (m)', coords(row)));
                    ylabel(sprintf('%s (m)', coords(col)));
                    
                    subplot(312)
                    imagesc_centered(xyzPos{row}{2}, xyzPos{col}{2}, ...
                        transpose(squeeze(data(:,:,:,2))));
                    axis image
                    colorbar;
                    title(file.Fields{2}.Name);
                    xlabel(sprintf('%s (m)', coords(row)));
                    ylabel(sprintf('%s (m)', coords(col)));
                    
                    subplot(313)
                    imagesc_centered(xyzPos{row}{3}, xyzPos{col}{3}, ...
                        transpose(squeeze(data(:,:,:,3))));
                    axis image
                    colorbar;
                    title(file.Fields{3}.Name);
                    xlabel(sprintf('%s (m)', coords(row)));
                    ylabel(sprintf('%s (m)', coords(col)));
                    %suptitle(sprintf('Frame %i', frameNum));
                end
                pause(0.01);
            end
            if frameNum < numFrames
                data = file.readFrames(1);
            end
            frameNum = frameNum + 1;
        end
        file.close
    end
elseif (nnz(dim(1:3) == 1) == 0)   % All 3D cases
    error('Cannot plot 3D data with plotOutput (yet).');
end

