function writeSourceSpec(sim, source, varargin)

X.AutoMaskFile = '';
X.AutoTimeFile = '';
X.AutoSpaceTimeFile = '';
X = parseargs(X, varargin{:});


% The purpose of the spec file is mostly to provide the DIMENSIONS of data.
% This will be equivalent to source.yeeCells for space-varying data, but will
% not be equivalent for time data, for instance.

%if isfield(source, 'timeData') && ~isempty(source.timeData)
if ~isempty(X.AutoTimeFile)
    fh = fopen([X.AutoTimeFile, '.txt'], 'wt');
    fprintf(fh, 'trogdor5data\n');
    fprintf(fh, 'precision %s\n', sim.Precision);
    fprintf(fh, 'date %s\n', date());
    fprintf(fh, 'dxyz [%g, %g, %g]\n', ...
        sim.Dxyz(1), sim.Dxyz(2), sim.Dxyz(3));
    fprintf(fh, 'dt %g\n', sim.Dt);
    fprintf(fh, 'origin [%i, %i, %i]\n', sim.Grid.Origin);
    fprintf(fh, 'specfile %s\n', [X.AutoTimeFile, '.txt']);
    
    % CustomTFSFSource must provide all E and H fields and has no "field".
    if isfield(source, 'field')
        fieldList = source.field;
    else
        fieldList = {'ex', 'ey', 'ez', 'hx', 'hy', 'hz'};
    end
    
    for ff = 1:length(fieldList)
        offset = t6.fieldOffset(fieldList{ff});
        fprintf(fh, 'field %s [%g, %g, %g] %g\n', fieldList{ff}, offset);
    end
    
    for dd = 1:size(source.timesteps, 1)
        fprintf(fh, 'duration from %i to %i period 1\n', source.timesteps(dd,:));
    end
    
    if isfield(source, 'bounds') && ~isempty(source.bounds)
        for bb = 1:size(source.bounds, 1)
            fprintf(fh, 'region [[%i, %i, %i], [%i, %i, %i]] stride [1, 1, 1]',...
                source.yeeCells(bb,:));
            fprintf(fh, ' bounds [[%i, %i, %i], [%i, %i, %i]]\n', ...
                source.bounds(bb,:))
        end
    else
        for bb = 1:size(source.yeeCells, 1)
            fprintf(fh, 'region [[%i, %i, %i], [%i, %i, %i]] stride [1, 1, 1]\n',...
                source.yeeCells(bb,:));
        end
    end
    
    fclose(fh);
end

%if isfield(source, 'spaceTimeData') && ~isempty(source.spaceTimeData)
if ~isempty(X.AutoSpaceTimeFile)
    fh = fopen([X.AutoSpaceTimeFile, '.txt'], 'wt');
    fprintf(fh, 'trogdor5data\n');
    fprintf(fh, 'precision %s\n', sim.Precision);
    fprintf(fh, 'date %s\n', date());
    fprintf(fh, 'dxyz [%g, %g, %g]\n', ...
        sim.Dxyz(1), sim.Dxyz(2), sim.Dxyz(3));
    fprintf(fh, 'dt %g\n', sim.Dt);
    fprintf(fh, 'origin [%i, %i, %i]\n', sim.Grid.Origin);
    fprintf(fh, 'specfile %s\n', [X.AutoTimeFile, '.txt']);
    %fprintf(fh, 'unitVector0 (1, 0, 0)\n');
    %fprintf(fh, 'unitVector1 (0, 1, 0)\n');
    %fprintf(fh, 'unitVector2 (0, 0, 1)\n');
    
    % CustomTFSFSource must provide all E and H fields and has no "field".
    if isfield(source, 'field')
        fieldList = source.field;
    else
        fieldList = {'ex', 'ey', 'ez', 'hx', 'hy', 'hz'};
    end
    
    for ff = 1:length(fieldList)
        offset = t6.fieldOffset(fieldList{ff});
        fprintf(fh, 'field %s [%g, %g, %g] %g\n', fieldList{ff}, offset);
    end
    
    for dd = 1:size(source.timesteps, 1)
        fprintf(fh, 'duration from %i to %i period 1\n', source.timesteps(dd,:));
    end
    
    if ~isempty(source.bounds)
        for bb = 1:size(source.bounds, 1)
            fprintf(fh, 'region [[%i, %i, %i], [%i, %i, %i]] stride [1, 1, 1]',...
                source.yeeCells(bb,:));
            fprintf(fh, ' bounds [[%i, %i, %i], [%i, %i, %i]]\n', ...
                source.bounds(bb,:));
        end
    else
        for bb = 1:size(source.yeeCells, 1)
            fprintf(fh, 'region [[%i, %i, %i], [%i, %i, %i]] stride [1, 1, 1]\n',...
                source.yeeCells(bb,:));
        end
    end
    
    %for yy = 1:size(source.yeeCells, 1)
    %    fprintf(fh, 'region [[%i, %i, %i], [%i, %i, %i]] stride [1, 1, 1]\n',...
    %        source.yeeCells(yy,:));
    %end
    
    fclose(fh);
end

if isfield(source, 'maskFile')
    %warning('Not writing specfile for maskFile (do with data request?)');
end

if isfield(source, 'spaceTimeFile')
    %warning('Not writing specfile for spaceTimeFile (do with data request?)');
end
