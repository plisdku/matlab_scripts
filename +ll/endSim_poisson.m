function disjointMeshes = endSim_poisson(varargin)
    % endSim
    
    X.MPH = 'fromMatlab.mph';
    X.StopEarly = false;
    X.SaveFields = false;
    X.Gradient = true;
    X.Solver = 'Direct';
    X = parseargs(X, varargin{:});

    global LL_MODEL;
    bounds = LL_MODEL.bounds;

    import com.comsol.model.*
    import com.comsol.model.util.*
    
    model = ModelUtil.create('Model');
    model.modelNode.create('mod1');
    geom = model.geom.create('geom1', 3);
    

    %% Figure out whether we can read things out of a file or not.
    
    [stepFile, movableDomainsFile, domainMaterialsFile] = ...
        generateFileNamesFromGeometry(LL_MODEL.meshes, bounds);

    % I'm not using the cache at present---though I should---but here is
    % the test for its existence.
    cacheExists = exist(stepFile, 'file') && exist(movableDomainsFile, 'file') ...
        && exist(domainMaterialsFile, 'file');
    cacheExists = false;
    warning('Assuming cache does not exist');
    
    %% Create the materials

    comsolMaterials(model, LL_MODEL.meshes, LL_MODEL.materials);
    numMaterials = numel(unique(cellfun(@(x) x.material, LL_MODEL.meshes)));
    
    %% Create output/source/measurement volumes, planes and points!
    % Outputs, sources and measurements can all be on points, volumes or planes.
    % They must be incorporated with the geometry somehow.
    
    sourceStructs = makeSourcesOrMeasurements(model, geom, 'src', LL_MODEL.sources);
    measStructs = makeSourcesOrMeasurements(model, geom, 'meas', LL_MODEL.measurements);
    
    %% Build the STEP file if necessary
    % It's needed unless it's cached and I'm using the cache.
    
    % chunks: array of structs with vertices, faces and a material tag.
    % they are mutually disjoint and their material tags may not be unique.
    [disjointMeshes, chunks] = processGeometry(LL_MODEL.meshes, [sourceStructs measStructs], stepFile);
    
    comsolSTEPImport(geom, stepFile);
    
    model.save([pwd filesep 'preRunGeometry.mph']);
    comsolRunGeometry(geom);
    model.save([pwd filesep 'postRunGeometry.mph']);
    
    %% Assign a material to each domain
    % Get the material for each chunk, then get the chunk for each domain.
    
    comsolAssignMaterials(model, geom, chunks, ...
        cacheExists, domainMaterialsFile, numMaterials)
    
    %% Get the domain and boundary indices for each input mesh
    % We'll use these to select electrode boundaries to apply voltages on
    % and to exclude domains from meshing as needed (e.g. for electrodes).
    
    %[meshDomains, meshOuterBoundaries] = findMeshEntities(model, geom, LL_MODEL.meshes);
    [meshDomains, meshOuterBoundaries, meshNewBoundaries] = findMeshEntities(model, geom, disjointMeshes);
    
    % Create selections named 'boundary_%i' for each mesh.
    % These are for *new* boundaries added by each input mesh, so each
    % boundary entity is present in one and only one selection.
    comsolBoundarySelections(model, meshNewBoundaries);

    %% Physics!
    
    elementOrderString = '3';
    warning('Element order %s', elementOrderString);
        
    comsolForwardPhysics(model, LL_MODEL.meshes, elementOrderString);
    
    comsolAdjointPhysics(model, LL_MODEL.meshes, LL_MODEL.measurements, ...
        measStructs, elementOrderString);
    
    %% View!

    model.view('view1').set('renderwireframe', false);
    
    %% Study!
    
    comsolStudy(model, X.Gradient);
    
    %% Surface outputs
    
    comsolMovableMeshSelection(model, LL_MODEL.meshes, cacheExists, ...
        movableDomainsFile);
    
    comsolMeasurements(model, LL_MODEL.measurements, X.Gradient);
    
    %% Variables for handling the adjoint calculation
    comsolVariables(model);
    
    %% Plots
    comsolPlots(X, model);
    
    %% Mesh!
    % I'm building the mesh here, right before the solve step, because
    % sometimes the solution fails due to meshing problems.  I want to
    % respond to solution failures by tweaking the mesh.
    
    % For every mesh attempt
    %  if it fails, tweak some parameters sorta randomly
    %
    % For ever solution failure
    %  make the mesh a bit smaller
    
    comsolMesh(model, meshDomains, LL_MODEL.meshes, ...
        LL_MODEL.measurements, ...
        LL_MODEL.hmax, LL_MODEL.hmin, LL_MODEL.hgrad);
    
    model.save([pwd filesep X.MPH]);
    
    %% Save the model before running
    
    model.save([pwd filesep X.MPH]);
    if X.StopEarly
        model.save([pwd filesep 'saveEarly.mph']);
        warning('quitting early')
        return
    end

    %% Run it all
    
    maxSolveAttempts = 5;
    %hmaxes = success.hmax*linspace(1.0, 0.75, maxSolveAttempts);
    
    succeeded = 0;
    tries = 1;
    while ~succeeded
        try
            %model.sol('sol1').runAll;
            
            % Forward: run Study Step 1 through Save Solutions 1
            model.sol('sol1').runFromTo('st1', 'su1');
            
            % INSERT CHARGED PARTICLE OPTICS HERE
            % Also... run the field export here as needed!
            
            % Dual: run Study Step 2 through Stationary 2
            model.sol('sol1').runFromTo('st2', 's2');
            
            succeeded = 1;
            if X.SaveFields
                model.save([pwd filesep 'fields_' X.MPH]);
            end
        catch crapception
            keyboard
            %{
            if tries < maxSolveAttempts
                warning('Solve failed.  Re-meshing and re-solving.');
                tries = tries + 1;
                
                attemptMeshing(model, theMesh, hmaxes(tries), success.hgrad, ...
                    maxMeshAttempts);
            else
                error('Solution failed with %i tries', tries);
            end
            %}
        end
    end

    if tries > 1
        fprintf('Solve succeeded with %i tries\n', tries);
    end
    
    if isa(LL_MODEL.forwardCallback, 'function_handle')
        LL_MODEL.forwardCallback()
    end
    
    %% Save the objective function and sensitivity information

    model.result.table('tblF').clearTableData;
    model.result.numerical('intF').set('table', 'tblF');
    model.result.numerical('intF').setResult;

    if X.Gradient
        model.result.export('exportSurfaceDF').run;
    end

    model.result.export('expTableF').run;
    
    %% Save other requested exports
    for nn = 1:length(LL_MODEL.exports)
    if strcmpi(LL_MODEL.exports{nn}.mode, 'Forward')
        doExport(model, LL_MODEL.exports{nn}, sprintf('export_%i',nn));
    end
    end

end

function doExport(model, exportStruct, export_name)
    
    pointFile = sprintf([pwd filesep '%s_points.txt'], export_name);
    assert(size(exportStruct.points,2) == 3);
    dlmwrite(pointFile, exportStruct.points, 'delimiter', '\t');
    
    export = model.result.export.create(export_name, 'Data');
    export.label(export_name);
    export.set('location', 'file');
    %export.set('descr', {'Electric potential'});
    export.set('filename', exportStruct.file);
    %export.set('unit', {'V'});
    export.set('coordfilename', [pwd filesep pointFile]);
    export.set('expr', {'V'});
%    model.save([pwd filesep 'mid_export.mph']);
    try
        export.run;
    catch exc
        warning('Got a prollem')
    end
end



function [stepFile, movableDomainsFile, domainMaterialsFile] = ...
        generateFileNamesFromGeometry(meshes, bounds)
% Calculate checksum from the geometry and generate filenames with it.

    checksum = geometryChecksum(meshes, bounds);
    assert(isa(checksum, 'uint16'));

    stepFile = sprintf('structure_%.4x.step', checksum);
    movableDomainsFile = sprintf('movableDomains_%.4x.txt', checksum);
    domainMaterialsFile = sprintf('domainMaterials_%.4x.txt', checksum);
end


function comsolMaterials(model, meshes, materials)

    tensorElems = @(T) arrayfun(@(a) sprintf('%2.8f+%2.8fi',real(a),imag(a)), T(:), ...
        'UniformOutput', false);

    numMats = numel(unique(cellfun(@(x) x.material, meshes)));

    for mm = 1:numMats

        matName = sprintf('mat%i', mm);
        mat = model.material.create(matName);
        mat.name(matName);
        mat.propertyGroup('def').set('relpermittivity', ...
            tensorElems(materials{mm}.epsr * eye(3)));
        mat.propertyGroup('def').set('relpermeability', ...
            tensorElems(materials{mm}.mur * eye(3)));
        mat.propertyGroup('def').set('electricconductivity', ...
            tensorElems(materials{mm}.sigma * eye(3)));

    end

end

function comsolSTEPImport(geom, stepFileName)

    stepImport = geom.feature.create('impSTEP', 'Import');
    stepImport.set('createselection', true);
    stepImport.set('type', 'cad');
    stepImport.set('filename', [pwd filesep stepFileName]);
	stepImport.set('importtol', '1.0E-9');
    %stepImport.set('unit', 'source');

    % The STEP file will be in millimeters.  Scale to meters.
    % Pre-scaling results in internal geometry errors since COMSOL sucks.
    scale = geom.feature.create('STEP_unit_correction', 'Scale');
    scale.set('isotropic', '1e3');
    scale.set('factor', '1e3');
    scale.selection('input').named('impSTEP');

end


function comsolRunGeometry(geom)
    try
        geom.run();
    catch exc
        warning('Failed initial attempt');

        % Now do pairwise unions of things from nonPMLChunks.
        % The hope is that some union of two chunks will prevent geometry
        % failures.

        pairs = combnk(1:3,2);

        succeeded = false;
        for pp = 1:size(pairs,1)
        if succeeded == false
            fprintf('Attempting with union of %i and %i\n', ...
                pairs(pp,1), pairs(pp,2));

            try
                un = geom.feature.create('saveMyButtUnion', 'Union');
                un.selection('input').set(...
                    { sprintf('sca1(%i)', pairs(pp,1)), ...
                      sprintf('sca1(%i)', pairs(pp,2)) } );
                model.save([pwd filesep 'preRunGeometry.mph']);
                geom.run();
                succeeded = true;
                fprintf('Succeded!\n');
            catch exc
                warning('Failed geom.run()!');
                geom.feature.remove('saveMyButtUnion');
            end
        end
        end

        if succeeded == false
            error('Failed every single time!');
        end
    end
end

function [meshDomains, meshBoundaries, contributedBoundaries] = ...
    findMeshEntities(model, geom, meshes)
    % The original structural meshs, which may have several chunks
    domainMeshes = smallestEnclosingChunks(model, geom, meshes); % used in loop
    meshDomains = cell(size(meshes)); % for meshing
    meshBoundaries = cell(size(meshes)); % comsolElectrodeSelections
    contributedBoundaries = cell(size(meshes));
    
    for cc = 1:length(meshes)
        meshDomains{cc} = find(domainMeshes == cc);
        meshBoundaries{cc} = ll.outerDomainBoundaryEntities(model, meshDomains{cc});
    end
    
    % Now I want something else... the boundary entities that were added
    % by the addition of each mesh in order.  Gotta work backwards to get
    % these.
    
    boundariesToExclude = [];
    for cc = length(meshes):-1:1
        contributedBoundaries{cc} = setdiff(meshBoundaries{cc},...
            boundariesToExclude);
        boundariesToExclude = [boundariesToExclude; ...
            contributedBoundaries{cc}];
    end
    
end

function comsolAssignMaterials(model, geom, chunks, ...
        cacheExists, domainMaterialsFile, numMaterials)
    % Seems like I don't define materials at all now.  How to handle this?
    if cacheExists
        domainMaterial = dlmread(domainMaterialsFile);
    else
        domainChunks = smallestEnclosingChunks(model, geom, chunks);
        %chunkBoundaries = cell(size(chunks)); % not used
        %for cc = 1:length(chunks)
        %    chunkBoundaries{cc} = ll.outerDomainBoundaryEntities(model, find(domainChunks == cc));
        %end

        allMaterials = cellfun(@(a) a.material, chunks);
        domainMaterial = allMaterials(domainChunks);
        dlmwrite(domainMaterialsFile, domainMaterial);
    end

    for mm = 1:numMaterials
        matTag = sprintf('mat%i', mm);
        selTag = sprintf('selMat%i', mm);
        sel = model.selection.create(selTag, 'Explicit');
        sel.set(find(domainMaterial == mm));
        model.material(matTag).selection.named(selTag);
    end
end

function comsolBoundarySelections(model, meshBoundaries)
    for mm = 1:numel(meshBoundaries)
        selName = sprintf('boundary_%i', mm);
        sel = model.selection.create(selName, 'Explicit');
        sel.name(selName);
        sel.geom('geom1', 2);
        sel.set(meshBoundaries{mm});
    end
end
% 
% function comsolElectrodeSelections(model, meshes, meshBoundaries)
%     for mm = 1:numel(meshes)
%     if ~isempty(meshes{mm}.voltage)
%         electrodeName = sprintf('electrode_%i', mm);
%         sel = model.selection.create(electrodeName, 'Explicit');
%         sel.name(electrodeName);
%         sel.geom('geom1', 2);
%         sel.set(meshBoundaries{mm});
%     end
%     end
% end
    
function comsolForwardPhysics(model, meshes, elementOrderString)

    model.save([pwd filesep 'prePhysics.mph']);
    
    model.physics.create('es', 'Electrostatics', 'geom1');
    model.physics('es').prop('ShapeProperty').set('order_electricpotential',...
        elementOrderString);

    % Forward electric potentials!
    % (Boundary conditions with surface charge don't need to create any
    % object analogous to this potential object... just write to rhoq.)
    for mm = 1:numel(meshes)
    if ~isempty(meshes{mm}.voltage)
        potentialName = sprintf('potential%i', mm);
        pot = model.physics('es').create(potentialName, ...
            'ElectricPotential', 2);
        pot.selection.named(sprintf('boundary_%i', mm));
        pot.set('V0', meshes{mm}.voltage);
        pot.name(sprintf('Potential %i', mm));
    elseif ~isempty(meshes{mm}.surfacecharge)
        chargeName = sprintf('surfaceCharge%i', mm);
        charge = model.physics('es').create(chargeName, 'SurfaceChargeDensity', 2);
        charge.selection.named(sprintf('boundary_%i', mm));
        charge.set('rhoqs', meshes{mm}.surfacecharge);
    end
    end
end


function comsolAdjointPhysics(model, meshes, measurements, measStructs, ...
    elementOrderString)
    
    %fprintf('Adjoint physics\n')
    model.physics.create('es2', 'Electrostatics', 'geom1');
    model.physics('es2').prop('ShapeProperty').set('order_electricpotential',...
        elementOrderString);

    % Adjoint electric potentials
    % Every surface must have V=0 on it.

    for mm = 1:numel(meshes)
    if ~isempty(meshes{mm}.voltage)
        potentialName = sprintf('adjointPotential%i', mm);
        pot = model.physics('es2').create(potentialName, ...
            'ElectricPotential', 2);
        pot.selection.named(sprintf('boundary_%i', mm));
        pot.set('V0', 0);
        pot.name(sprintf('Adjoint potential %i', mm));
    elseif ~isempty(meshes{mm}.surfacecharge)
        chargeName = sprintf('adjointSurfaceCharge%i', mm);
        charge = model.physics('es2').create(chargeName, ...
            'SurfaceChargeDensity', 2);
        charge.selection.named(sprintf('boundary_%i', mm));
        charge.set('rhoqs', 0);
    end
    end

    % Adjoint space charge!

    numMeasurements = numel(measurements);
    
    measDims = [];
    measSel = {};
    for ss = 1:numMeasurements
        
        measDims(ss) = measurements{ss}.dimensions;
        
        %bounds = measurements{ss}.bounds;
        %extent = bounds(4:6) - bounds(1:3);
        

        %if ~X.Gradient
        %    continue;
        %end

        currName = sprintf('adjSpaceCharge%i', ss);

        if measurements{ss}.dimensions == 0 % point
            error('not handling points yet')
        elseif measurements{ss}.dimensions == 2 % surface current
            error('not handling surfaces yet')
        elseif measurements{ss}.dimensions == 3 % volume current

            scd = model.physics('es2').create(currName, 'SpaceChargeDensity', 3);
            scd.selection.named(measStructs{ss}.selectionName);
            %scd.set('rhoq', LL_MODEL.measurements{ss}.rhoq);
            scd.set('rhoq', measurements{ss}.g); % should it be g or not?
            scd.name(sprintf('Adjoint charge %i', ss));

            measSel = {measSel{:} measStructs{ss}.selectionName};
        end
    end

    if numel(unique(measDims)) > 1
        error('Mixing measurement dimensions!');
    end

    if numMeasurements > 0
        measurementSel = model.selection.create('measSel', 'Union');
        measurementSel.geom('geom1', measDims(1));
        measurementSel.name('Measurement selection');
        measurementSel.set('input', measSel);
    end
end


function comsolMovableMeshSelection(model, meshes, cacheExists, movableDomainsFile)
    if cacheExists
        %fprintf('Using cached movable domains file.\n');
        try
            movableMeshDomains = dlmread(movableDomainsFile);
        catch exc
            warning('Could not read movable domains file, is it empty?');
            movableMeshDomains = [];
        end
    else
        movableMeshDomains = findMovableBoundaries(model, meshes);
        dlmwrite(movableDomainsFile, movableMeshDomains);
    end

    sel = model.selection.create('movableMeshes', 'Explicit');
    sel.name('Movable meshes');
    sel.geom('geom1', 2);
    sel.set(movableMeshDomains);
end

% Create the geometrical representation for a source or measurement
% manifold.  3D makes a triangulated boundary representation, and 2D makes
% a rect in a COMSOL work plane.  Also returns the number of dimensions
% (why not) and the COMSOL selection name.  The selections are actually
% created in this function, which seems out of step with how I do things
% elsewhere.
function srcMeasRecords = makeSourcesOrMeasurements(model, geom,...
        prefix, srcMeasInputs)
    
    N = numel(srcMeasInputs);
    srcMeasRecords = cell(size(srcMeasInputs));

    for nn = 1:N
        bounds = srcMeasInputs{nn}.bounds;
        extent = bounds(4:6) - bounds(1:3);


        if nnz(extent) == 0 % it's a point
            
            ptName = sprintf('pt_%s_%i', prefix, nn);
            pt = geom.feature.create(ptName, 'Point');
            pt.setIndex('p', num2str(bounds(1)), 0, 0);
            pt.setIndex('p', num2str(bounds(2)), 1, 0);
            pt.setIndex('p', num2str(bounds(3)), 2, 0);
            pt.set('createselection', true);
            
            selectionName = sprintf('geom1_%s_pnt', ptName);
            
        elseif nnz(extent) == 1 % it's a line
            error('I do not handle lines yet');
        elseif nnz(extent) == 2 % it's a plane
            if bounds(1) == bounds(4)
                plane = 'yz'; quickPlane = 'quickx';
                sz = bounds([5 6]) - bounds([2 3]);
                center = 0.5*(bounds([5 6]) + bounds([2 3]));
                quickPlaneCoord = bounds(1);
            elseif bounds(2) == bounds(5)
                plane = 'zx'; quickPlane = 'quicky';
                sz = bounds([6 4]) - bounds([3 1]);
                center = 0.5*(bounds([6 4]) + bounds([3 1]));
                quickPlaneCoord = bounds(2);
            elseif bounds(3) == bounds(6)
                plane = 'xy'; quickPlane = 'quickz';
                sz = bounds([4 5]) - bounds([1 2]);
                center = 0.5*(bounds([4 5]) + bounds([1 2]));
                quickPlaneCoord = bounds(3);
            end

            wpName = sprintf('wp_%s_%i', prefix, nn);
            wp = geom.feature.create(wpName, 'WorkPlane');
            wp.set('quickplane', plane);
            wp.set(quickPlane, quickPlaneCoord);
            wp.geom.feature.create('r1', 'Rectangle');
            wp.geom.feature('r1').set('base', 'center');
            wp.geom.feature('r1').set('size', {num2str(sz(1)), num2str(sz(2))});
            wp.geom.feature('r1').set('pos', {num2str(center(1)), num2str(center(2))});
            wp.set('createselection', true);

            selectionName = sprintf('geom1_%s_bnd', wpName);

        elseif nnz(extent) == 3 % it's a volume

            % The point here is just to mark off a chunk of space and force the
            % tetrahedra to line up with its surfaces.

            r = t6.model.Rect(@(p) bounds);
            m = r.meshes;
            srcMeasRecords{nn}.vertices = m{1}.patchVertices;
            srcMeasRecords{nn}.faces = m{1}.faces;

            boxName = sprintf('box_%s_%i', prefix, nn);
            
            % Create the selection of the inside elements.
            margin = 1e-6;
            box = model.selection.create(boxName, 'Box');
            box.set('xmin', num2str(bounds(1)-margin));
            box.set('ymin', num2str(bounds(2)-margin));
            box.set('zmin', num2str(bounds(3)-margin));
            box.set('xmax', num2str(bounds(4)+margin));
            box.set('ymax', num2str(bounds(5)+margin));
            box.set('zmax', num2str(bounds(6)+margin));
            box.set('condition', 'inside');
            box.name(boxName);

            selectionName = boxName;

        else
            error('wut.');
        end

        srcMeasRecords{nn}.dimensions = nnz(extent);
        srcMeasRecords{nn}.selectionName = selectionName;
    end

end


% For each mesh in meshes, return a mesh that is also disjoint to all the meshes
% after it:
%
% disjointMeshes{N} = meshes{N}
% disjointMeshes{N-1} = meshes{N-1} - meshes{N};
% disjointMeshes{N-2} = meshes{N-2} - meshes{N-1} - meshes{N};
% . . .
%
% An optional bounding box argument may be specified as well.  Any mesh that is
% contained in the bounding box will be excluded from these operations.  The
% use of this is to allow the huge number of PML intersections that are
% calculated later to be done without involving geometry that's never going to
% touch the PML anywhere.
function disjointMeshes = makeDisjointInputs(meshes, excludeBounds)
    
    myPatch = @(v,f,c) patch('Vertices', v, 'Faces', f, 'FaceColor', 'c',...
        'EdgeAlpha', 0.1, 'FaceAlpha', 0.3);
    
    bbox = @(vertArray) [min(vertArray) max(vertArray)];
    rectInRect = @(r1, r2) all(r1(1:3) >= r2(1:3)) && all(r1(4:6) <= r2(4:6));
        
    if nargin < 2
        % Make a bounding box that EVERYTHING reaches outside.
        nullBounds = [Inf Inf Inf -Inf -Inf -Inf];
        excludeBounds = nullBounds;
    end
    
    numMeshes = numel(meshes);
    disjointMeshes = cell(size(meshes));
    
    % now figure out the disjointeries
    for mm = 1:numMeshes
    if ~rectInRect(bbox(meshes{mm}.vertices), excludeBounds)
    
        v = meshes{mm}.vertices;
        f = meshes{mm}.faces;
        
        for nn = (mm+1):numMeshes
        if ~rectInRect(bbox(meshes{nn}.vertices), excludeBounds)
            
            v2 = meshes{nn}.vertices;
            f2 = meshes{nn}.faces;
            
            [v, f] = neflab.nefDifference(v, f, v2, f2);
        end
        end
        
        disjointMeshes{mm}.vertices = v;
        disjointMeshes{mm}.faces = f;
        disjointMeshes{mm}.material = meshes{mm}.material;
        
    end
    end
end

function meshes = uniteMaterials(inMeshes)
    
    meshes = {};
    
    numInputs = numel(inMeshes);
    for mm = 1:numInputs
    if ~isempty(inMeshes{mm})
        
        iMat = inMeshes{mm}.material;
        
        if iMat <= numel(meshes) && isfield(meshes{iMat}, 'vertices')
            % If the material has been unioned before, unite with previous.
            [meshes{iMat}.vertices, meshes{iMat}.faces] = neflab.nefUnion(...
                inMeshes{mm}.vertices, inMeshes{mm}.faces, ...
                meshes{iMat}.vertices, meshes{iMat}.faces);
        else
            % New material, not yet unioned.  Make a new output mesh.
            meshes{iMat}.vertices = inMeshes{mm}.vertices;
            meshes{iMat}.faces = inMeshes{mm}.faces;
            meshes{iMat}.material = iMat;
        end
        
    end
    end
    
end

% outChunks is basically the Venn diagram from overlaying subChunks on
% inChunks.  The materials in outChunks come from the materials in
% inChunks.
function outChunks = vennChunks(inChunks, subChunks)
    
    outChunks = inChunks;
    for iSub = 1:numel(subChunks)
    if subChunks{iSub}.dimensions == 3
        
        numInChunks = numel(outChunks);
        nNew = numInChunks + 1;
        
        for iIn = 1:numInChunks
            
            [vInter, fInter] = neflab.nefIntersection(...
                outChunks{iIn}.vertices, outChunks{iIn}.faces, ...
                subChunks{iSub}.vertices, subChunks{iSub}.faces);
            
            if ~isempty(vInter)
                
                [vDiff, fDiff] = neflab.nefDifference(...
                    outChunks{iIn}.vertices, outChunks{iIn}.faces, ...
                    subChunks{iSub}.vertices, subChunks{iSub}.faces);
                
                if ~isempty(vDiff)
                    outChunks{nNew}.vertices = vInter;
                    outChunks{nNew}.faces = fInter;
                    outChunks{nNew}.material = outChunks{iIn}.material;
                    nNew = nNew + 1;
                    
                    outChunks{iIn}.vertices = vDiff;
                    outChunks{iIn}.faces = fDiff;
                    
                end
            end
        end
        
    end
    end
end

function r = calcBoundingBox(A)
    r = [min(A) max(A)];
end

% I ought to use mphselectcoords to find the indices of all boundaries that
% are within some little tolerance of a movable vertex.  Perhaps.  Eh?  I
% think that's a better way.
function movableMeshDomains = findMovableBoundaries(model, meshes)
    movableMeshDomains = [];
    
    %calcBoundingBox = @(A) [min(A) max(A)];
    
    for mm = 1:numel(meshes)    
    if any(meshes{mm}.jacobian(:))
    
        bbox = calcBoundingBox(meshes{mm}.vertices) + [-1 -1 -1 1 1 1]*1e-6;
    
        comsolIndices = [1 4; 2 5; 3 6];
    
        boundaryDomains = mphselectbox(model, 'geom1', bbox(comsolIndices),...
            'boundary');
    
        movableMeshDomains = [movableMeshDomains boundaryDomains];
    
    end
    end
end

% function electrodeDomains = findElectrodes(model, meshes)
%     electrodeDomains = [];
%     
%     for mm = 1:numel(meshes)
%         if ~isempty(meshes{mm}.voltage)
%             electrodeDomains = [electrodeDomains, mm];
%         end
%     endf
% end

function enclosingChunks = smallestEnclosingChunks(model, geom, chunks)
    % Approximately find the smallest chunk enclosing each domain in the
    % COMSOL geometry.  This is probably the chunk of material that the
    % domain belongs to.
    %
    % Method: Each chunk in Matlab may enclose several COMSOL domains.
    % The domain most likely belongs to the smallest chunk which encloses
    % it.  So: for each domain, find its smallest enclosing chunk.  Then
    % assign the material of that chunk to the domain.
    %
    % Actually we don't find the domains enclosed by each chunk, we find
    % the domains enclosed by each chunk's bounding box.  This is WRONG but
    % tends to work for structures with lots of axis-aligned boxes.
    
    numChunks = numel(chunks);
    
    % 1. Find all bounding boxes and their volumes
    
    enclosures = sparse(numChunks, geom.getNDomains);
    
    volumes = zeros(1,numChunks);
    
    for ss = 1:numChunks
        bbox = calcBoundingBox(chunks{ss}.vertices);
        volumes(ss) = prod(bbox(4:6)-bbox(1:3));
        
        domains = selectBoxHelper(model, bbox + [-1 -1 -1 1 1 1]*1e-6);
        enclosures(ss,domains) = 1;
    end
    
    % Find the chunk with the smallest bounding box for each domain.
    [enChunk, enDomain] = find(enclosures);
    
    [~,smallestChunks] = max(sparse(enChunk, enDomain, 1./volumes(enChunk)), [], 1);
    
    enclosingChunks = smallestChunks;
end

function [disjointMeshes, nonPMLChunks] = processGeometry(meshes, srcMeasStructs, stepFile)
    % chunks = processGeometry(meshes, srcMeasStructs, stepFile)
    %
    % 

    %% Create mutually disjoint input meshes!
    % Each mesh subtracts off all previous meshes.
    
    disjointMeshes = makeDisjointInputs(meshes);
    %fprintf('Done with the difference operations.\n');
    
    % Make similar disjoint meshes but skip everything that does not reach PML.
    % This will really speed things up when intersecting every PML block with every
    % material block.
    
    pl = @(mesh) flatPatch('Vertices', mesh.vertices, 'Faces', mesh.faces,...
        'FaceColor', 'g', 'EdgeAlpha', 0.1, 'FaceAlpha', 0.2);
    
    srcMeasMeshes = gatherSrcMeasMeshes(srcMeasStructs);

    %% Structure not in PML!
    
    nonPMLChunks = uniteMaterials(disjointMeshes);
    
    %% Adjust for measurements!
    %
    
    if ~isempty(srcMeasStructs)
        nonPMLChunks = vennChunks(nonPMLChunks, srcMeasStructs);
    end
    
    %% Create the STEP file and set up STEP import.
    
    %fprintf('Got to the STEP file.\n');
    writeSTEP(nonPMLChunks, stepFile);
    
    % Just for testing!  If needed.
    function assertDisjoint(chunks)
        for cc = 1:numel(chunks)
            for dd = cc+1:numel(chunks)

                fprintf('testing %i and %i...\n', cc, dd);
                assert(~neflab.nefTestIntersection(...
                    chunks{cc}.vertices, chunks{cc}.faces, ...
                    chunks{dd}.vertices, chunks{dd}.faces));

            end
        end
        fprintf('ALL ARE DISJOINT\n');
    end
    %assertDisjoint(nonPMLChunks);
    
end

function srcMeasMeshes = gatherSrcMeasMeshes(srcMeasStructs)

    srcMeasMeshes = {};
    nOut = 0;
    for nn = 1:numel(srcMeasStructs)
        if srcMeasStructs{nn}.dimensions == 3
            srcMeasMeshes{nOut+1}.faces = srcMeasStructs{nn}.faces;
            srcMeasMeshes{nOut+1}.vertices = srcMeasStructs{nn}.vertices;
        end
        nOut = nOut + 1;
    end

end


function chunkFiles = writeSTEP(chunks, stepFile)
    
    if ~exist('importMeshes', 'dir')
        mkdir('importMeshes');
    end
    
    chunkFiles = arrayfun(@(ii) sprintf('importMeshes/chunk%i.stl',ii), ...
        1:numel(chunks), 'UniformOutput', false);
    
    for curChunk = 1:numel(chunks)
        fpath = [pwd filesep chunkFiles{curChunk}];
        writeSTL(chunks{curChunk}.vertices, chunks{curChunk}.faces, fpath);
    end
    
    spacedNames = cellfun(@(A) [A ' '], chunkFiles, 'UniformOutput', false);

    callMerge = ['mergeSTP -unit mm ',...
        spacedNames{:}, ...
        ' > mergeSTP.txt'];
    unix(callMerge);
    unix(sprintf('mv outStep.step %s', stepFile));
    
end

% Get a checksum for the whole geometry including bounds and materials.
% If I've already done all the hard work before, intersections and all, I don't
% need to do it again.
function checksum = geometryChecksum(meshes, bounds)
    
    ll.fletcher16(bounds);
    
    for mm = 1:numel(meshes)
        ll.fletcher16(meshes{mm}.vertices, true);
        ll.fletcher16(meshes{mm}.faces, true);
        ll.fletcher16(meshes{mm}.material, true);
    end
    
    checksum = ll.fletcher16(bounds);
end


function selections = selectBoxHelper(model, bbox)

    indexer = [1 4; 2 5; 3 6];
    selections = mphselectbox(model, 'geom1', bbox(indexer), 'domain');

end


% Create the Study node in the simulation.
function comsolStudy(model, doCalculateGradient)
    
%    global LL_MODEL;

    fprintf('Solution\n')

    model.study.create('std1');
    
    % Create forward study
    model.study('std1').create('stat', 'Stationary');
    model.study('std1').feature('stat').set('activate', {'es' 'on' 'es2' 'off'});
    model.study('std1').feature('stat').label('Forward');
    model.study('std1').feature('stat').set('notlistsolnum', 1);
    model.study('std1').feature('stat').set('notsolnum', '1');
    model.study('std1').feature('stat').set('listsolnum', 1);
    model.study('std1').feature('stat').set('solnum', '1');
    
    % Create the solver, which also creates the "dset1" dataset.
    model.sol.create('sol1');
    
    model.sol('sol1').study('std1');
    model.sol('sol1').attach('std1');
    
    model.sol('sol1').create('st1', 'StudyStep');
    model.sol('sol1').feature('st1').set('study', 'std1');
    model.sol('sol1').feature('st1').set('studystep', 'stat');

    model.sol('sol1').create('v1', 'Variables');
    model.sol('sol1').feature('v1').set('control', 'stat');
    
    model.sol('sol1').create('s1', 'Stationary');
    model.sol('sol1').feature('s1').create('fc1', 'FullyCoupled');
    model.sol('sol1').feature('s1').create('i1', 'Iterative');
    model.sol('sol1').feature('s1').feature('i1').set('linsolver', 'cg');
    model.sol('sol1').feature('s1').feature('i1').create('mg1', 'Multigrid');
    model.sol('sol1').feature('s1').feature('i1').feature('mg1').set('prefun', 'amg');
    model.sol('sol1').feature('s1').feature('fc1').set('linsolver', 'i1');
    model.sol('sol1').feature('s1').feature.remove('fcDef');
    
    % The Store Solution node makes the forward solution's fields
    % accessible to the dual solution for calculation of sources etc.
    % Creating this node also automatically creates the "dset2" dataset.
    model.sol('sol1').create('su1', 'StoreSolution');
    
    model.study('std1').create('stat1', 'Stationary');
    model.study('std1').feature('stat1').set('activate', {'es' 'off' 'es2' 'on'});
    model.study('std1').feature('stat1').label('Dual');
    model.study('std1').feature('stat1').set('notlistsolnum', 1);
    model.study('std1').feature('stat1').set('notsolnum', 'auto');
    model.study('std1').feature('stat1').set('listsolnum', 1);
    model.study('std1').feature('stat1').set('solnum', 'auto');
    
    if doCalculateGradient
        model.sol('sol1').create('st2', 'StudyStep');
        model.sol('sol1').feature('st2').set('study', 'std1');
        model.sol('sol1').feature('st2').set('studystep', 'stat1');

        model.sol('sol1').create('v2', 'Variables');
        model.sol('sol1').feature('v2').set('initmethod', 'sol');
        model.sol('sol1').feature('v2').set('initsol', 'sol1');
        model.sol('sol1').feature('v2').set('notsolmethod', 'sol');
        model.sol('sol1').feature('v2').set('notsol', 'sol1');
        model.sol('sol1').feature('v2').set('notsoluse', 'su1');
        model.sol('sol1').feature('v2').set('control', 'stat1');
        model.sol('sol1').feature('v2').set('solnum', 'auto');
        model.sol('sol1').feature('v2').set('solvertype', 'solnum');
        model.sol('sol1').feature('v2').set('listsolnum', {'1'});
        model.sol('sol1').feature('v2').set('solnum', 'auto');
        
        model.sol('sol1').create('s2', 'Stationary');
        model.sol('sol1').feature('s2').create('fc1', 'FullyCoupled');
        model.sol('sol1').feature('s2').create('i1', 'Iterative');
        model.sol('sol1').feature('s2').feature('i1').set('linsolver', 'cg');
        model.sol('sol1').feature('s2').feature('i1').create('mg1', 'Multigrid');
        model.sol('sol1').feature('s2').feature('i1').feature('mg1').set('prefun', 'amg');
        model.sol('sol1').feature('s2').feature('fc1').set('linsolver', 'i1');
        model.sol('sol1').feature('s2').feature.remove('fcDef');
    end
end

function comsolMesh(model, meshDomains, meshes, measurements, hmax, hmin, hgrad)
    theMesh = model.mesh.create('mesh1', 'geom1');

    nonElectrodeDomains = [];
    for mm = 1:length(meshes)
    if meshes{mm}.exclude == 0
        nonElectrodeDomains = [nonElectrodeDomains, meshDomains{mm}];
    end
    end
    comsolMeshParameters(theMesh, model, measurements, nonElectrodeDomains, ...
        meshes, hmin, hmax, hgrad);

    %globalSize = theMesh.feature('size');

    globalHmax = 20e-3;
    if ~isempty(hmax)
        globalHmax = str2double(hmax); % adjust for Mesh/Study
    end

    globalHmin = .1e-3;
    if ~isempty(hmin)
        globalHmin = str2double(hmin);
    end

    globalHgrad = 2;
    if ~isempty(hgrad)
        globalHgrad = str2double(hgrad);
    end

    maxMeshAttempts = 20;
    [meshingSucceeded, success.hmax, success.hgrad] = ...
        attemptMeshing(model, theMesh, globalHmax, globalHmin, globalHgrad, maxMeshAttempts);

    if ~meshingSucceeded
        error('Could not mesh!');
    end
end


function comsolMeshParameters(theMesh, model, measurements, domains, ...
    meshes, hmin, hmax, hgrad)
    %global LL_MODEL;
    
    globalSize = theMesh.feature('size');
    globalSize.set('custom', 'on');
    
    if ~isempty(hgrad)
        globalSize.set('hgrad', hgrad); % adjust for Mesh
    end
    
    if ~isempty(hmin)
        globalSize.set('hmin', hmin);
    end
    
    if ~isempty(hmax)
        globalSize.set('hmax', hmax);
    end
    
    for mm = 1:numel(meshes)
    if ~isempty(meshes{mm}.hmax) || ...
            ~isempty(meshes{mm}.hgrad) || ...
            ~isempty(meshes{mm}.hmin)

        bbox = calcBoundingBox(meshes{mm}.vertices) + ...
            [-1 -1 -1 1 1 1]*1e-6;

        comsolIndices = [1 4; 2 5; 3 6];

        boundaryDomains = mphselectbox(model, 'geom1', bbox(comsolIndices),...
            'boundary');

        szName = sprintf('size%i', mm);
        sz = theMesh.feature.create(szName, 'Size');
        sz.selection.geom('geom1', 2);
        sz.selection.set(boundaryDomains);
        sz.set('custom', 'on');
        
        if ~isempty(meshes{mm}.hmax)
            sz.set('hmaxactive', 'on');
            sz.set('hmax', meshes{mm}.hmax);
        end
        
        if ~isempty(meshes{mm}.hmin)
            sz.set('hminactive', 'on');
            sz.set('hmin', meshes{mm}.hmin);
        end
        
        if ~isempty(meshes{mm}.hgrad)
            sz.set('hgradactive', 'on');
            sz.set('hgrad', '1.5');
        end

    end 
    end
    
    % Change mesh size for measurement.
    sz = theMesh.feature.create('measSize', 'Size');
    sz.selection.geom('geom1', measurements{1}.dimensions);
    sz.selection.named('measSel');
    sz.set('custom', 'on');
    sz.set('hmaxactive', 'on');
    sz.set('hmax', 2e-3);
    %warning('Measurement hmax is 15');

    theMesh.feature.create('ftri1', 'FreeTri');
    theMesh.feature('ftri1').selection.named('movableMeshes');
    theMesh.feature.create('ftet1', 'FreeTet');
    theMesh.feature('ftet1').selection.geom('geom1', 3);
    theMesh.feature('ftet1').selection.set(domains);
end


function comsolPlots(X, model)
    
    model.result.dataset.create('movableDataset', 'Solution');
    model.result.dataset('movableDataset').name('Movable meshes');
    model.result.dataset('movableDataset').selection.named('movableMeshes');
    
    model.result.create('pg1', 'PlotGroup3D');
    model.result('pg1').label('Electric Potential (es)');
    model.result('pg1').set('oldanalysistype', 'noneavailable');
    model.result('pg1').set('frametype', 'spatial');
    model.result('pg1').set('data', 'dset1');
    model.result('pg1').feature.create('mslc1', 'Multislice');
    model.result('pg1').feature('mslc1').set('oldanalysistype', 'noneavailable');
    model.result('pg1').feature('mslc1').set('data', 'parent');
    
    model.result.create('pg2', 'PlotGroup3D');
    model.result('pg2').label('Electric Potential (es2)');
    model.result('pg2').set('oldanalysistype', 'noneavailable');
    model.result('pg2').set('frametype', 'spatial');
    model.result('pg2').set('data', 'dset1');
    model.result('pg2').feature.create('mslc1', 'Multislice');
    model.result('pg2').feature('mslc1').set('oldanalysistype', 'noneavailable');
    model.result('pg2').feature('mslc1').set('expr', 'V2');
    model.result('pg2').feature('mslc1').set('data', 'parent');

    if X.Gradient
        model.result.create('pg3', 'PlotGroup3D');
        model.result('pg3').feature.create('surf1', 'Surface');
        model.result('pg3').feature.create('arws1', 'ArrowSurface');
        model.result('pg3').feature.create('arws2', 'ArrowSurface');
        model.result('pg3').set('data', 'dset2');
        
        model.result('pg3').feature('surf1').set('expr', 'DF');
        model.result('pg3').feature('surf1').set('unit', 'W/m^3');
        model.result('pg3').feature('surf1').set('descr', '');
        model.result('pg3').feature('surf1').set('colortablerev', true);
        model.result('pg3').feature('surf1').set('colortablesym', true);
        model.result('pg3').feature('surf1').set('colortable', 'OrangeCrush');

        model.result('pg3').feature('arws1').set(...
            'expr', {'nx*DF*(DF<0)' 'ny*DF*(DF<0)' 'nz*DF*(DF<0)'});
        model.result('pg3').feature('arws1').set('arrowbase', 'head');
        %model.result('pg3').feature('arws1').set('scale', '3.2787942186813154E-10');
        %model.result('pg3').feature('arws1').set('scaleactive', false);

        model.result('pg3').feature('arws2').set(...
            'expr', {'nx*DF*(DF>0)' 'ny*DF*(DF>0)' 'nz*DF*(DF>0)'});
        model.result('pg3').feature('arws2').set('arrowbase', 'tail');
        model.result('pg3').feature('arws2').set('color', 'blue');
        %model.result('pg3').feature('arws2').set('scale', '4.240512841916214E-10');
        %model.result('pg3').feature('arws2').set('scaleactive', false);
    end
    

    %%
end


function comsolMeasurements(model, measurements, doCalculateGradient)
    
    bounds = measurements{1}.bounds;
    extent = bounds(4:6) - bounds(1:3);
    measDims = nnz(extent);
    
    %global LL_MODEL;
    
    % The surfaces data set should export the dual pressure DF at
    % all points on all surfaces.
    dsetSurfacesOnly = model.result.dataset.create('dsetSurfaces', 'Solution');
    dsetSurfacesOnly.name('Surfaces data set');
    dsetSurfacesOnly.selection.geom('geom1', 2);
    dsetSurfacesOnly.selection.named('movableMeshes');
%    dsetSurfacesOnly.selection.all;
    
    measData = model.result.dataset.create('dsetMeas', 'Solution');
    measData.name('Objective function data set');
    measData.selection.geom('geom1', measDims(1));
    measData.selection.named('measSel');

    assert(numel(measurements) == 1);
    tblF = model.result.table.create('tblF', 'Table');

    if measDims(1) == 0
        intPt = model.result.numerical.create('intF', 'EvalPoint');
        intPt.selection.named('measSel');
        intPt.set('probetag', 'none');
        intPt.set('table', 'tblF');
        intPt.set('expr', measurements{1}.F);
        
    elseif measDims(1) == 2
        intSurf = model.result.numerical.create('intF', 'IntSurface');
        intSurf.name('F');
        intSurf.selection.named('measSel');
        intSurf.set('probetag', 'none');
        intSurf.set('table', 'tblF');
        intSurf.set('expr', measurements{1}.F);
        
    elseif measDims(1) == 3
        intVol = model.result.numerical.create('intF', 'IntVolume');
        intVol.name('F');
        intVol.selection.named('measSel');
        intVol.set('probetag', 'none');
        intVol.set('table', 'tblF');
        intVol.set('expr', measurements{1}.F);
    end

    if doCalculateGradient
        model.result.export.create('exportSurfaceDF', 'Data');
        model.result.export('exportSurfaceDF').name('data set');
        model.result.export('exportSurfaceDF').set('data', 'dsetSurfaces');
        model.result.export('exportSurfaceDF').set('descr', {''});
        model.result.export('exportSurfaceDF').set('filename', 'DF_on_surfaces.txt');
        model.result.export('exportSurfaceDF').set('expr', {'mod1.DF', 'nx', 'ny', 'nz'});
        model.result.export('exportSurfaceDF').set('resolution', 'custom');
        model.result.export('exportSurfaceDF').set('lagorder', '5');
    end

    model.result.export.create('expTableF', 'tblF', 'Table');
    model.result.export('expTableF').set('filename', 'F.txt');
end


function comsolVariables(model)
    model.variable.create('var1');
    model.variable('var1').model('mod1');
    model.variable('var1').set('En', 'es.nx*es.Ex+es.ny*es.Ey+es.nz*es.Ez');
    model.variable('var1').set('E2n', 'es2.nx*es2.Ex+es2.ny*es2.Ey+es2.nz*es2.Ez');
    model.variable('var1').set('gg', '1', 'dual volume source');
    model.variable('var1').set('ee', '-D_eta*(-En)', 'primal boundary source');
    model.variable('var1').set('D_eta', 'motion_x*es.nx+motion_y*es.ny+motion_z*es.nz', 'boundary displ.');
    model.variable('var1').set('I_dual_bdy', '(-E2n)*ee', 'dual integrand, bdy');
    model.variable('var1').set('DF', '(-E2n)*(-En)', 'dual pressure');
    model.variable('var1').set('I_primal_vol', 'gg*V', 'primal integrand, vol');
    model.variable('var1').set('motion_x', '0', 'swept motion vec');
    model.variable('var1').set('motion_y', '0');
    model.variable('var1').set('motion_z', '1');
    model.variable('var1').set('a0', '1', 'const for F2');
    model.variable('var1').set('ax', '1', 'const for F3');
    model.variable('var1').set('ay', '0.5', 'const for F3');
    model.variable('var1').set('F_1', 'V', 'integrand for F1');
    model.variable('var1').set('F_2', 'a0*V', 'integrand for F2');
    model.variable('var1').set('F_3', '-ax*es.Ex-ay*es.Ey', 'integrand for F3');
    model.variable('var1').set('g_1', '1');
    model.variable('var1').set('g_2', 'a0');
    model.variable('var1').set('g_3', '0', '-div a = 0');
end


% Helper function.  When meshing or solving fails, step down through this
% master list of hmax and hgrad to try.
function [hmaxes, hmins, hgrads] = meshParameterAttempts(...
    hmax, hmin, hgrad, maxAttempts)
    
    n = ceil(sqrt(maxAttempts));
    hmaxes = hmax*linspace(1.0, 0.75, n);
    hgrads = hgrad*linspace(1.0, 0.5, n);
    hmins = hmin*linspace(1.0, 0.5, n);
    
    [opts.hmaxes, opts.hmins, opts.hgrads] = ndgrid(hmaxes, hmins, hgrads);
    [ii, jj, kk] = ndgrid(1:n, 1:n, 1:n);
    [~,ordering] = sort(ii(:)+jj(:)+kk(:));
    
    hmaxes = opts.hmaxes(ordering(1:maxAttempts));
    hmins = opts.hmins(ordering(1:maxAttempts));
    hgrads = opts.hgrads(ordering(1:maxAttempts));

end




function [meshingSucceeded, outHmax, outHgrad] = ...
    attemptMeshing(model, theMesh, hmax, hmin, hgrad, maxAttempts)
    
try
    [hmaxes, hmins, hgrads] = meshParameterAttempts(hmax, hmin, hgrad, maxAttempts);
catch exc
    keyboard
end
    globalSize = theMesh.feature('size');
    meshingSucceeded = false;
    attempts = 1;
    while ~meshingSucceeded && attempts < maxAttempts
        try
            attempts = attempts + 1;
            globalSize.set('hmax', num2str(hmaxes(attempts)));
            globalSize.set('hmin', num2str(hmins(attempts)));
            globalSize.set('hgrad', num2str(hgrads(attempts)));
            model.save([pwd filesep 'premesh.mph']);
            theMesh.run();
            meshingSucceeded = true;
        catch exc
            warning('Meshing attempt %i failed!\n', attempts);
            
            if attempts < maxAttempts
                fprintf('Retrying with hmax %f, hgrad %f\n', ...
                    hmaxes(attempts+1), hgrads(attempts+1));
            end
        end
    end

    if ~meshingSucceeded
        error('Meshing failed!');
    end
    
    outHmax = hmaxes(attempts);
    outHgrad = hgrads(attempts);
end





