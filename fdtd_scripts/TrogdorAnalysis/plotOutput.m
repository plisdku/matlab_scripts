function plotOutput(filePrefix, period)
% plotOutput(filePrefix, period)
%   Plot 1D and 2D output data from FDTD program.  1D data will be
%   plotted vs time, and 2D output will be presented frame by frame
%   without pause.  Use getOutputFrames or openOutputFile and
%   getOutputFrame to work with the data.
%
%   filePrefix  prefix for .dat and .txt FDTD output files
%   
%   period      set to integer n >= 1 to skip frames in time sequences
%
%   version 4.5
%   July 29, 2008


if (nargin < 2)
    period = 1;
end

[fid, dim] = openOutputFile(filePrefix);

    
disp(sprintf('Opened file %s with id %i', [filePrefix, '.dat'], fid));

% Dimension cases
%
% 1 1 1         One field componen at a point: time trace
% 1 1 1 3       Three components at a point: threefold time trace
% 1 1 N, perms  1D output: plot 1D movie
% 1 1 N 3       1D output: plot threefold 1D movie
% M N 1, perms  2D output: plot 2D movie
% M N 1 3       2D output: plot 2D movie with three parts

if (prod(dim(1:3)) == 1)   % All 0D cases
        [data, count] = getOutputFrames(fid, dim);
    if (length(dim) == 3)
        plot(squeeze(data));
    elseif (length(dim) == 4)  % assuming dim4 = 3
        plot(squeeze(data(1,1,1,1,:)), 'r');
        hold on
        plot(squeeze(data(1,1,1,2,:)), 'g');
        plot(squeeze(data(1,1,1,3,:)), 'b');
        hold off
        legend('Component 1', 'Component 2', 'Component 3');
    end
elseif (nnz(dim(1:3) == 1) == 2)   % All 1D cases
    if (length(dim) == 3)
        frameNum = 1;
        [data, count] = getOutputFrame(fid, dim);
        while (count ~= 0)
            if (mod(frameNum, period) == 0)
                plot(squeeze(data));
                ylim([-2, 2]);
                title(sprintf('Frame %i', frameNum));
                pause(0.1);
            end
            [data, count] = getOutputFrame(fid, dim);
            frameNum = frameNum + 1;
        end
    elseif (length(dim) == 4)
        frameNum = 1;
        [data, count] = getOutputFrame(fid, dim);
        while (count ~= 0)
            if (mod(frameNum, period) == 0)
                plot(squeeze(data(:,:,:,1)), 'r');
                hold on
                plot(squeeze(data(:,:,:,2)), 'g');
                plot(squeeze(data(:,:,:,3)), 'b');
                hold off
                ylim([-2, 2]);
                title(sprintf('Frame %i', frameNum));
                legend('Component 1', 'Component 2', 'Component 3');
                pause(0.1);
                %pause  
            end
            [data, count] = getOutputFrame(fid, dim);
            frameNum = frameNum + 1;
        end
    end
elseif (nnz(dim(1:3) == 1) == 1)    % All 2D cases
    if (length(dim) == 3)
        frameNum = 1;
        [data, count] = getOutputFrame(fid, dim);
        
        while (count ~= 0 && feof(fid) == 0)
            if (mod(frameNum, period) == 0)
                imagesc_centered(squeeze(data));
                axis image
                colorbar;
                title(sprintf('Frame %i', frameNum));
                pause(0.1);
            end
            [data, count] = getOutputFrame(fid, dim);
            frameNum = frameNum + 1;
        end
    elseif (length(dim) == 4)
        frameNum = 1;
        [data, count] = getOutputFrame(fid, dim);
        
        
        while (count ~= 0 && feof(fid) == 0)
            if (mod(frameNum, period) == 0)
                datsize = size(squeeze(data(:,:,:,1)));
                
                %if 1
                if datsize(1) > datsize(2)   % tall skinny
                    subplot('position', [0.05 0.05 0.3 0.9]);
                    imagesc_centered(squeeze(data(:,:,:,1)));
                    axis image
                    colorbar;
                    subplot('position', [0.36 0.05 0.3 0.9]);
                    imagesc_centered(squeeze(data(:,:,:,2)));
                    set(gca, 'ytick', []);
                    axis image
                    colorbar;
                    title(sprintf('Frame %i', frameNum));
                    subplot('position', [0.67 0.05 0.3 0.9]);
                    imagesc_centered(squeeze(data(:,:,:,3)));
                    set(gca, 'ytick', []);
                    axis image
                    colorbar;
                else
                    subplot('position', [0.05 0.05 0.9 0.3]);
                    imagesc_centered(squeeze(data(:,:,:,1)));
                    axis image
                    colorbar;
                    title(sprintf('Frame %i', frameNum));
                    subplot('position', [0.05 0.36 0.9 0.3]);
                    imagesc_centered(squeeze(data(:,:,:,2)));
                    set(gca, 'xtick', []);
                    axis image
                    colorbar;
                    subplot('position', [0.05 0.67 0.9 0.3]);
                    imagesc_centered(squeeze(data(:,:,:,3)));
                    set(gca, 'xtick', []);
                    axis image
                    colorbar;
                end
                %{
                if datsize(1) > datsize(2)
                    plotnums = [131 132 133];
                    cbarpos = 'SouthOutside';
                else
                    plotnums = [311 312 313];
                    cbarpos = 'EastOutside';
                end
                
                subplot(plotnums(1))
                imagesc_centered(squeeze(data(:,:,:,1)));
                colorbar(cbarpos);
                title(sprintf('Frame %i', frameNum));
                subplot(plotnums(2))
                imagesc_centered(squeeze(data(:,:,:,2)));
                colorbar(cbarpos);
                subplot(plotnums(3));
                imagesc_centered(squeeze(data(:,:,:,3)));
                colorbar(cbarpos)
                %}
                pause(0.1);
            end
            [data, count] = getOutputFrame(fid, dim);
            frameNum = frameNum + 1;
        end
    end
elseif (nnz(dim(1:3) == 1) == 0)   % All 3D cases
    error('Cannot plot 3D data with plotOutput (yet).');
end

%%%%
fclose(fid);

