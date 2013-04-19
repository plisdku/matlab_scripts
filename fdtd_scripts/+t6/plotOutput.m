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
X.Colormap = [];
X.UnitString = [];
X.FilePattern = '';
X.Subplots = [];
X.Times = [];
X.Callback = [];
X = parseargs(X, varargin{:});

if isempty(X.Colormap) && exist('orangecrush', 'file')
    colormap orangecrush;
elseif ~isempty(X.Colormap)
    colormap(X.Colormap);
end

imagesc_args = {};
if ~isempty(X.CLim)
    imagesc_args = {X.CLim};
end

file = OutputFile(fileName);

if file.numRegions() > 1
    error('Function does not yet work with multi-region outputs.');
end

if file.numFields() > 1
    xyzPos = file.positions('Field', 1);
    dim = [cellfun(@numel, xyzPos); file.numFields()];
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
    handle0D(file, X);
elseif (nnz(dim(1:3) == 1) == 2)   % All 1D cases
    handle1D(file, X);
elseif (nnz(dim(1:3) == 1) == 1)    % All 2D cases
    handle2D(file, X, imagesc_args);
elseif (nnz(dim(1:3) == 1) == 0)   % All 3D cases
    error('Cannot plot 3D data with plotOutput (yet).');
end

end


function handle0D(file, X)

fieldNames = {};
for nn = 1:file.numFields()
    fieldNames{nn} = file.Fields{nn}.Name;
end

file.open();
data = file.readFrames('Times', X.Times);
file.close();

fileTimes = {};
for ff = 1:file.numFields
    fileTimes{ff} = file.times('Field', ff);

    ft(:,ff) = fileTimes{ff}(1:size(data,5));
end

if ~isempty(X.Times)
    plot(X.Times, squish(data));
else
    %for ff = 1:file.numFields
        %plot(fileTimes{ff}, squish(data(:,:,:,ff,:)));
        plot(ft, squish(data)');
    %    hold on
    %end
end

xlabel('Time')
ylabel('Field')
legend(fieldNames{:});
%lineColors('jet')
    
end

function handle1D(file, X)

fieldNames = {};
for nn = 1:file.numFields()
    fieldNames{nn} = file.Fields{nn}.Name;
end

xyz = file.positions();
numFrames = file.numFramesAvailable;

figureHandle = gcf;

frameNum = 1;
file.open
data = file.readFrames('NumFrames', 1);
while frameNum <= numFrames
    
    if ~ishandle(figureHandle)
        return;
    end
    
    if (mod(frameNum, X.Period) == 0)
        
        plot(squeeze(data(:,:,:,:)));
        ylim(X.YLim);
        legend(fieldNames{:})
        title(sprintf('Frame %i of %i', frameNum, numFrames));
        if ~isempty(X.Callback)
            X.Callback();
        end
        pause(0.01);
        
        if frameNum < numFrames
            data = file.readFrames('NumFrames', 1);
        end
        frameNum = frameNum + 1;
        
    end
end
file.close;


end



function nxny = bestSubplots(Lx, Ly, numFields)

% Try to keep the total plot extent small.

desiredAspectRatio = 1.5;
desiredAngle = atan(1/desiredAspectRatio);

smallestArea = inf;
bestAngle = inf;

bestRowCol = [1 1];

for numRows = 1:numFields
for numCols = 1:numFields
if numRows*numCols >= numFields    % this is such a dumb way to do this :-D
    
    xSize = Lx*numCols;
    ySize = Ly*numRows;
    
    area = xSize*ySize;
    
    thisAngle = atan(ySize/xSize);
    
    if area < smallestArea
        bestRowCol = [numRows numCols];
        bestAngle = thisAngle;
        smallestArea = area;
    elseif area == smallestArea
        
        if abs(thisAngle-desiredAngle) < abs(bestAngle-desiredAngle)
            bestRowCol = [numRows numCols];
            bestAngle = thisAngle;
            smallestArea = area;
        end
        
    end
    
end
end
end

nxny = bestRowCol;

end





function handle2D(file, X, imagesc_args)

numFrames = file.numFramesAvailable();

% Determine which axes we're using
bounds = file.Regions.Bounds(1,:);

if bounds(1) == bounds(4)
    row = 2; col = 3;
elseif bounds(2) == bounds(5)
    row = 1; col = 3;
else
    row = 1; col = 2;
end


xyzPos = file.positions();

if isempty(X.Subplots)
    Lx = bounds(row+3) - bounds(row);
    Ly = bounds(col+3) - bounds(col);
    
    nxny = bestSubplots(Lx, Ly, file.numFields);
else
    nxny = X.Subplots;
end

movieFrame = 0;

figureHandle = gcf;

file.open;
for frame = 0:numFrames-1
    
    if ~ishandle(figureHandle)
        return;
    end
    
    data = file.readFrames('NumFrames', 1);
    
    if mod(frame, X.Period) == 0
        for ff = 1:file.numFields
            
            if prod(nxny) > 1
                subplot(nxny(2), nxny(1), ff, 'align');
            end
            
            myRow = floor(1 + (ff-1)/nxny(1));
            myCol = 1 + mod(ff-1, nxny(1));
            
            imagesc_centered(xyzPos{row,1}, xyzPos{col,1}, ...
                transpose(squeeze(data(:,:,:,ff))), imagesc_args{:});
            axis xy image
            colorbar
            title(file.Fields{ff}.Name);
            
            if myRow < nxny(2)
                set(gca, 'xtick', []);
            else
                xlabel(sprintf('%s %s', char('w'+row), X.UnitString));
            end
            
            if myCol > 1
                set(gca, 'ytick', []);
            else
                ylabel(sprintf('%s %s', char('w'+col), X.UnitString));
            end
            
        end
        
        if ~isempty(X.FilePattern)
            saveas(gcf, sprintf(X.FilePattern, movieFrame));
            movieFrame = movieFrame + 1;
        end
        
        if ~isempty(X.Callback)
            X.Callback();
        end
        
        pause(0.01);
    end
end
file.close

end

