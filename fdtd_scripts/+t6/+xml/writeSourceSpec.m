function writeSourceSpec(source, varargin)
sim = t6.TrogdorSimulation.instance; % to get dxyz, dt, etc.

X.AutoMaskFile = '';
X.AutoTimeFile = '';
X.AutoSpaceTimeFile = '';
X = parseargs(X, varargin{:});

precisionString = t6.TrogdorSimulation.instance().Precision;

% The purpose of the spec file is mostly to provide the DIMENSIONS of data.
% This will be equivalent to source.yeeCells for space-varying data, but will
% not be equivalent for time data, for instance.

if isfield(source, 'timeData')
    fh = fopen([X.AutoTimeFile, '.txt'], 'wt');
    fprintf(fh, 'trogdor5data\n');
    fprintf(fh, 'precision %s\n', precisionString);
    fprintf(fh, 'date %s\n', date());
    fprintf(fh, 'dxyz [%g, %g, %g]\n', ...
        sim.Dxyz(1), sim.Dxyz(2), sim.Dxyz(3));
    fprintf(fh, 'dt %g\n', sim.Dt);
    fprintf(fh, 'specfile %s\n', [X.AutoTimeFile, '.txt']);
%    fprintf(fh, 'unitVector0 (1, 0, 0)\n');
%    fprintf(fh, 'unitVector1 (0, 1, 0)\n');
%    fprintf(fh, 'unitVector2 (0, 0, 1)\n');
    
    % CustomTFSFSource must provide all E and H fields and has no "field".
    if isfield(source, 'field')
        fieldList = source.field;
    else
        fieldList = {'ex', 'ey', 'ez', 'hx', 'hy', 'hz'};
    end
    
    for ff = 1:length(fieldList)
        offset = t6.xml.fieldOffset(fieldList{ff});
        fprintf(fh, 'field %s [%g, %g, %g] %g\n', fieldList{ff}, offset);
    end
    
    for dd = 1:size(source.duration, 1)
        fprintf(fh, 'duration from %i to %i period 1\n', source.duration(dd,:));
    end
    
    fprintf(fh, 'region [[0, 0, 0], [0, 0, 0]] stride [1, 1, 1]\n');
    
    fclose(fh);
end

if isfield(source, 'spaceTimeData')
    fh = fopen([X.AutoSpaceTimeFile, '.txt'], 'wt');
    fprintf(fh, 'trogdor5data\n');
    fprintf(fh, 'precision %s\n', precisionString);
    fprintf(fh, 'date %s\n', date());
    fprintf(fh, 'dxyz [%g, %g, %g]\n', ...
        sim.Dxyz(1), sim.Dxyz(2), sim.Dxyz(3));
    fprintf(fh, 'dt %g\n', sim.Dt);
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
        offset = t6.xml.fieldOffset(fieldList{ff});
        fprintf(fh, 'field %s [%g, %g, %g] %g\n', fieldList{ff}, offset);
    end
    
    for dd = 1:size(source.duration, 1)
        fprintf(fh, 'duration from %i to %i period 1\n', source.duration(dd,:));
    end
    
    for yy = 1:size(source.yeeCells, 1)
        fprintf(fh, 'region [[%i, %i, %i], [%i, %i, %i]] stride [1, 1, 1]\n',...
            source.yeeCells(yy,:));
    end
    
    fclose(fh);
end

if isfield(source, 'maskData')
    fh = fopen([X.AutoMaskFile, '.txt'], 'wt');
    fprintf(fh, 'trogdor5data\n');
    fprintf(fh, 'precision %s\n', precisionString);
    fprintf(fh, 'date %s\n', date());
    fprintf(fh, 'dxyz [%g, %g, %g]\n', ...
        sim.Dxyz(1), sim.Dxyz(2), sim.Dxyz(3));
    fprintf(fh, 'dt %g\n', sim.Dt);
    fprintf(fh, 'specfile %s\n', [X.AutoMaskFile, '.txt']);
    fprintf(fh, 'unitVector0 (1, 0, 0)\n');
    fprintf(fh, 'unitVector1 (0, 1, 0)\n');
    fprintf(fh, 'unitVector2 (0, 0, 1)\n');
    
    % CustomTFSFSource must provide all E and H fields and has no "field".
    if isfield(source, 'field')
        fieldList = source.field;
    else
        fieldList = {'ex', 'ey', 'ez', 'hx', 'hy', 'hz'};
    end
    
    for ff = 1:length(fieldList)
        offset = t6.xml.fieldOffset(fieldList{ff});
        fprintf(fh, 'field %s [%g, %g, %g] %g\n', fieldList{ff}, offset);
    end
    
    fprintf(fh, 'duration from 1 to 1 period 1\n');
    
    for yy = 1:size(source.yeeCells, 1)
        fprintf(fh, 'region [[%i, %i, %i], [%i, %i, %i]] stride [1, 1, 1]\n',...
            source.yeeCells(yy,:));
    end
    fclose(fh);
end

if isfield(source, 'spaceTimeData')
    fh = fopen([X.AutoSpaceTimeFile, '.txt'], 'wt');
    fprintf(fh, 'trogdor5data\n');
    fprintf(fh, 'precision %s\n', precisionString);
    fprintf(fh, 'date %s\n', date());
    fprintf(fh, 'dxyz [%g, %g, %g]\n', ...
        sim.Dxyz(1), sim.Dxyz(2), sim.Dxyz(3));
    fprintf(fh, 'dt %g\n', sim.Dt);
    fprintf(fh, 'specfile %s\n', [X.AutoSpaceTimeFile, '.txt']);
%    fprintf(fh, 'unitVector0 (1, 0, 0)\n');
%    fprintf(fh, 'unitVector1 (0, 1, 0)\n');
%    fprintf(fh, 'unitVector2 (0, 0, 1)\n');
    
    % CustomTFSFSource must provide all E and H fields and has no "field".
    if isfield(source, 'field')
        fieldList = source.field;
    else
        fieldList = {'ex', 'ey', 'ez', 'hx', 'hy', 'hz'};
    end
    
    for ff = 1:length(fieldList)
        offset = t6.xml.fieldOffset(fieldList{ff});
        fprintf(fh, 'field %s [%g, %g, %g] %g\n', fieldList{ff}, offset);
    end
    
    for dd = 1:size(source.duration, 1)
        fprintf(fh, 'duration from %i to %i period 1\n', source.duration(dd,:));
    end
    
    for yy = 1:size(source.yeeCells, 1)
        fprintf(fh, 'region [[%i, %i, %i], [%i, %i, %i]] stride [1, 1, 1]\n',...
            source.yeeCells(yy,:));
    end
    fclose(fh);
end

if isfield(source, 'maskFile')
    %warning('Not writing specfile for maskFile (do with data request?)');
end

if isfield(source, 'spaceTimeFile')
    %warning('Not writing specfile for spaceTimeFile (do with data request?)');
end
