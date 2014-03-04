function endSim(varargin)
% endSim

X.MPH = 'fromMatlab.mph';
X.StopEarly = false;
X.SaveFields = false;
X = parseargs(X, varargin{:});

global LL_MODEL;

import com.comsol.model.*
import com.comsol.model.util.*

model = ModelUtil.create('Model');

bounds = LL_MODEL.bounds;
pmlBounds = LL_MODEL.PMLBounds;
pmlThickness = LL_MODEL.PMLThickness;

if ~exist('importMeshes', 'dir')
    mkdir('importMeshes');
end

%%

fprintf('Making geometry node\n')

model.modelNode.create('mod1');
geom = model.geom.create('geom1', 3);
geom.lengthUnit('nm');

%%
% The new way I need to do this is to mirror all the Boolean operations
% that COMSOL carries out, and based on expected answers I can then
% determine whether or not to do anything explicitly in COMSOL.
%
% So I guess it's like this:
%
%   Create the structures (inputs).  User must provide background rect.
%   Make them mutually disjoint.
%       for each input structure, take all appropriate differences
%       if the end result is non-null, implement all non-null differences
%   Calculate non-PML components
%       take intersection with non-PML
%       if non-null, mirror in COMSOL
%   Calculate PML rects
%       for each PML rect, take all appropriate intersections
%           implement each non-null intersection
%
% "Mirroring" means do it in COMSOL and make a note of the material to
% apply.

%% Create suitable input meshes!
% Make sure they're all mutually disjoint, yo.

numMeshes = numel(LL_MODEL.meshes);
materialIndex = zeros(1, numMeshes);

% Each mesh subtracts off all previous meshes.

disjointVertices = cell(numMeshes,1);
disjointFaces = cell(numMeshes,1);

disjointStructureIndices = [];

fprintf('Got to the disjoint meshes.\n');
%%

% now figure out the disjointeries
for mm = 1:numMeshes
    
    v = LL_MODEL.meshes{mm}.vertices;
    f = LL_MODEL.meshes{mm}.faces;
    
    listOfDifferences = [];
    
    for nn = (mm+1):numMeshes
        
%         figure(1); clf
%         quickPatch(v, f);
%         axis image; view(3);
        
        if neflab.nefTestIntersection(v, f, ...
            LL_MODEL.meshes{nn}.vertices, ...
            LL_MODEL.meshes{nn}.faces)
            
            [v2 f2] = neflab.nefDifference(v, f, ...
                LL_MODEL.meshes{nn}.vertices, ...
                LL_MODEL.meshes{nn}.faces);
            
            %fprintf('MIRROR: append difference for COMSOL\n');
            
            listOfDifferences(end+1) = nn;
            
            v = v2;
            f = f2;
            
        end
        
    end
    disjointVertices{mm} = v;
    disjointFaces{mm} = f;
    
    if numel(f) > 0
        disjointStructureIndices(end+1) = mm;
        
        %addImport(geom, v, f, comsolStructureName(mm), comsolStructureName(mm));
        materialIndex(mm) = LL_MODEL.meshes{mm}.material;
    else
        %fprintf('MIRROR: do not add differences to COMSOL\n');
    end
end

fprintf('Done with the difference operations.\n');

%% Create suitable PML meshes!

pml = pmlBounds;
boundsX = [pml(1) bounds(1) bounds(4) pml(4)];
boundsY = [pml(2) bounds(2) bounds(5) pml(5)];
boundsZ = [pml(3) bounds(3) bounds(6) pml(6)];

numPMLs = 0;

pmlVertices = {};
pmlFaces = {};

for xx = 1:3
for yy = 1:3
for zz = 1:3
if xx ~= 2 || yy ~= 2 || zz ~= 2
    
    rectBounds = [boundsX(xx) boundsY(yy) boundsZ(zz) ...
            boundsX(xx+1) boundsY(yy+1) boundsZ(zz+1)];
    
    rectSize = rectBounds(4:6) - rectBounds(1:3);
    
    if all(rectSize > 0)
        numPMLs = numPMLs + 1;

        r = t6.model.Rect(@(p) rectBounds);

        m = r.meshes;
        
        pmlVertices{numPMLs} = m{1}.patchVertices;
        pmlFaces{numPMLs} = m{1}.faces;
    end
    
end
end
end
end

%% Structure in PML!

fprintf('Intersecting with PML.\n');

pmlIndex = 0;

pmlChunkIndex = 0;
pmlChunkMaterials = [];
pmlChunkVertices = {};
pmlChunkFaces = {};

maxTotalChunks = 26*numMeshes;

for pmlIndex = 1:numPMLs
    
    for ss = 1:numMeshes
        [chunkVertices chunkFaces] = neflab.nefIntersection(...
            disjointVertices{ss}, disjointFaces{ss},...
            pmlVertices{pmlIndex}, pmlFaces{pmlIndex});
        
%         figure(3); clf
%         quickPatch(pmlVertices{pmlIndex+1}, pmlFaces{pmlIndex+1});
%         quickPatch(disjointVertices{ss}, disjointFaces{ss}, ...
%             'g');
%         axis([-100 1100 -100 1100 -100 1100])
        
        if numel(chunkFaces) > 0
            
            pmlChunkIndex = pmlChunkIndex + 1;
            
            pmlChunkVertices{pmlChunkIndex} = chunkVertices;
            pmlChunkFaces{pmlChunkIndex} = chunkFaces;
            pmlChunkMaterials(pmlChunkIndex) = ...
                LL_MODEL.meshes{ss}.material;
            
%             figure(3);
%             quickPatch(chunkVertices, chunkFaces, 'b');
            
%            addPMLIntersection(geom, ss, pmlIndex, pmlChunkIndex);
            
            %fprintf('MIRROR: perform PML intersection\n');
            % make sure to append material to material list
            % also increment COMSOL PML block counter I guess
            
        end
        %pause
    end
    
end

numPMLChunks = numel(pmlChunkVertices);


%% Structure not in PML!

r = t6.model.Rect(@(p) bounds);
m = r.meshes;
nonPMLVertices = m{1}.patchVertices;
nonPMLFaces = m{1}.faces;

interiorStructureVertices = {};
interiorStructureFaces = {};

fprintf('Intersecting with non-PML.\n');

for mm = 1:numMeshes
    [interiorStructureVertices{mm} interiorStructureFaces{mm}] = ...
        neflab.nefIntersection(disjointVertices{mm}, disjointFaces{mm},...
        nonPMLVertices, nonPMLFaces);
    
    fracDone = mm/numMeshes;
end

%% Create the STEP file.

fprintf('Got to the STEP file.\n');

totalChunks = numel(disjointFaces) + numel(pmlChunkVertices);

chunkFiles = arrayfun(@(ii) sprintf('chunk%i.stl',ii), ...
    1:totalChunks, 'UniformOutput', false);

curChunk = 1;

for iIn = 1:numel(disjointFaces)
    fpath = [pwd filesep chunkFiles{curChunk}];
    ll.writeSTL(interiorStructureVertices{iIn}, ...
        interiorStructureFaces{iIn}, fpath);
    curChunk = curChunk + 1;
end

for iOut = 1:numel(pmlChunkVertices)
    fpath = [pwd filesep chunkFiles{curChunk}];
    ll.writeSTL(pmlChunkVertices{iOut}, pmlChunkFaces{iOut}, fpath);
    curChunk = curChunk + 1;
end

catCell = @(A) A{:};

spacedNames = cellfun(@(A) [A ' '], chunkFiles, 'UniformOutput', false);

callMerge = ['mergeSTP -unit mm ', spacedNames{:}];
unix(callMerge);
% this created a file called outStep.step.

stepImport = geom.feature.create('impSTEP', 'Import');
stepImport.set('createselection', true);
stepImport.set('type', 'cad');
stepImport.set('filename', [pwd filesep 'outStep.step']);
%stepImport.set('unit', 'source');

% The STEP file will be in millimeters.  STEP can't do nanometers.
% So, I'll scale everything on the way in!
% Pre-scaling results in internal geometry errors since COMSOL sucks.
geom.feature.create('sca1', 'Scale');
geom.feature('sca1').set('isotropic', '1e-6');
geom.feature('sca1').set('factor', '1e-6');
geom.feature('sca1').selection('input').named('impSTEP');

%% Output planes!!

fprintf('Source and output planes.\n');

numOutputs = numel(LL_MODEL.outputs);
numSources = numel(LL_MODEL.sources);
numMeasurements = numel(LL_MODEL.measurements);

outPlaneNames = arrayfun(@(aa) sprintf('wp_out%i', aa), 1:numOutputs, ...
    'UniformOutput', false);
srcPlaneNames = arrayfun(@(aa) sprintf('wp_src%i', aa), 1:numSources, ...
    'UniformOutput', false);
measPlaneNames = arrayfun(@(aa) sprintf('wp_meas%i', aa), 1:numMeasurements, ...
    'UniformOutput', false);

outputsSources = {LL_MODEL.outputs{:} LL_MODEL.sources{:} ...
    LL_MODEL.measurements{:}};
planeNames = {outPlaneNames{:} srcPlaneNames{:} measPlaneNames{:}};

for pp = 1:numel(planeNames)
    bounds = outputsSources{pp}.bounds;
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
    
    wp = geom.feature.create(planeNames{pp}, 'WorkPlane');
    wp.set('quickplane', plane);
    wp.set(quickPlane, quickPlaneCoord);
    wp.geom.feature.create('r1', 'Rectangle');
    wp.geom.feature('r1').set('base', 'center');
    wp.geom.feature('r1').set('size', {num2str(sz(1)), num2str(sz(2))});
    wp.geom.feature('r1').set('pos', {num2str(center(1)), num2str(center(2))});
    wp.set('createselection', true);
    
    %waitbar(pp/numel(planeNames), h, 'Creating source and output planes');
end

%% run that geom

model.save([pwd filesep 'intermediate.mph']);
fprintf('Running the geometry!\n');
geom.run;

%%
fprintf('Got to the intermediate spot.\n');
model.save([pwd filesep 'intermediate.mph']);

%% Materials!

tensorElems = @(T) arrayfun(@(a) sprintf('%2.8f+%2.8fi',real(a),imag(a)), T(:), ...
    'UniformOutput', false);

fprintf('Creating materials.\n');

numMaterials = numel(LL_MODEL.materials);

for mm = 1:numMaterials
    
    matName = sprintf('mat%i', mm);
    mat = model.material.create(matName);
    mat.name(matName);
    mat.propertyGroup('def').set('relpermittivity', ...
        tensorElems(LL_MODEL.materials{mm}.epsr * eye(3)));
    mat.propertyGroup('def').set('relpermeability', ...
        tensorElems(LL_MODEL.materials{mm}.mur * eye(3)));
    mat.propertyGroup('def').set('electricconductivity', ...
        tensorElems(LL_MODEL.materials{mm}.sigma * eye(3)));
    
end
%    mat.selection.named(sprintf('geom1_imp%i_dom', nn));

%% Figure out which materials go where

calcBoundingBox = @(A) [min(A) max(A)];

movableMeshDomains = [];

materialByDomain = zeros(geom.getNDomains, 1);

for ss = 1:numMeshes
        
    bbox = calcBoundingBox(interiorStructureVertices{ss});
    domainNum = identifyByBounds(model, bbox);
    materialByDomain(domainNum) = LL_MODEL.meshes{ss}.material;
    
    % while we're at it, find and mark boundaries...
    if any(LL_MODEL.meshes{ss}.jacobian(:))
        if numel(domainNum) == 1
            boundaryDomains = mphgetadj(model, 'geom1', 'boundary', ...
                'domain', domainNum);
            movableMeshDomains = [movableMeshDomains boundaryDomains];
        end
    end
end

for cc = 1:numPMLChunks
    bbox = calcBoundingBox(pmlChunkVertices{cc});
    domainNum = identifyByBounds(model, bbox);
    materialByDomain(domainNum) = pmlChunkMaterials(cc);
end


%% Assign domains to materials

for mm = 1:numMaterials
    matTag = sprintf('mat%i', mm);
    selTag = sprintf('selMat%i', mm);
    sel = model.selection.create(selTag, 'Explicit');
    sel.set(find(materialByDomain == mm));
    model.material(matTag).selection.named(selTag);
end

%% Make selection of movable meshes

sel = model.selection.create('movableMeshes', 'Explicit');
sel.name('Movable meshes');
sel.geom('geom1', 2);
sel.set(movableMeshDomains);

%% Mark PML!

everythingBox = model.selection.create('box1', 'Box');
everythingBox.name('Everything');

nonPMLRect = LL_MODEL.bounds + [-1 -1 -1 1 1 1];

nonPMLBox = model.selection.create('box2', 'Box');
nonPMLBox.name('Non PML');
nonPMLBox.set('condition', 'inside');
nonPMLBox.set('xmin', nonPMLRect(1));
nonPMLBox.set('ymin', nonPMLRect(2));
nonPMLBox.set('zmin', nonPMLRect(3));
nonPMLBox.set('xmax', nonPMLRect(4));
nonPMLBox.set('ymax', nonPMLRect(5));
nonPMLBox.set('zmax', nonPMLRect(6));

pmlSel = model.selection.create('pmlSel', 'Difference');
pmlSel.name('PML');
pmlSel.set('add', {'box1'});
pmlSel.set('subtract', {'box2'});

pml = model.coordSystem.create('pml1', 'geom1', 'PML');
pml.selection.named('pmlSel');

%% Assembly!  Yarr!

model.save([pwd filesep 'foobar.mph']);
%warning('Quitting early, hardcoded.\n');
%return

%% Forward physics!

fprintf('Forward physics\n')

model.physics.create('emw', 'ElectromagneticWaves', 'geom1');
model.physics('emw').prop('ShapeProperty').set('order_electricfield', '3');

if ~isempty(LL_MODEL.incidentField.Ex) ||...
    ~isempty(LL_MODEL.incidentField.Ey) || ...
    ~isempty(LL_MODEL.incidentField.Ez)
    
    model.physics('emw').prop('BackgroundField').set('SolveFor', 'scatteredField');
    model.physics('emw').prop('BackgroundField').set('Eb', ...
        {LL_MODEL.incidentField.Ex; ...
        LL_MODEL.incidentField.Ey; ...
        LL_MODEL.incidentField.Ez});
end

%% Forward Current sources!

numSources = numel(LL_MODEL.sources);

for ss = 1:numSources
    planeName = sprintf('geom1_wp_src%i_bnd', ss);
        
    weakName = sprintf('weakSource%i', ss);
    surfCurr = model.physics('emw').feature.create(weakName, ...
        'WeakContribution', 2);
    surfCurr.selection.named(planeName);
    surfCurr.set('weakExpression', ...
        sprintf('-(%s)*test(emw.curlEx)-(%s)*test(emw.curlEy)-(%s)*test(emw.curlEz)-emw.iomega*mu0_const*((%s)*test(emw.Ex)+(%s)*test(emw.Ey)+(%s)*test(emw.Ez))', ...
        LL_MODEL.sources{ss}.Mx, ...
        LL_MODEL.sources{ss}.My, ...
        LL_MODEL.sources{ss}.Mz, ...
        LL_MODEL.sources{ss}.Jx, ...
        LL_MODEL.sources{ss}.Jy, ...
        LL_MODEL.sources{ss}.Jz));
    surfCurr.name(sprintf('Surface current %i', ss));
end

%% Adjoint physics!

fprintf('Adjoint physics\n')

model.physics.create('emw2', 'ElectromagneticWaves', 'geom1');
model.physics('emw2').prop('ShapeProperty').set('order_electricfield', '3');
model.physics('emw2').prop('BackgroundField').set('SolveFor', 'fullField');
%% Adjoint current sources!

numMeasurements = numel(LL_MODEL.measurements);

measurementSel = model.selection.create('measSel', 'Union');
measurementSel.geom('geom1', 2);
measurementSel.name('Measurement selection');
measSel = {};

for ss = 1:numMeasurements
    planeName = sprintf('geom1_wp_meas%i_bnd', ss);
    probeName = sprintf('probe%i', ss);
    
    currName = sprintf('adjCurrentSource%i', ss);
    surfCurr = model.physics('emw2').feature.create(currName, ...
        'SurfaceCurrent', 2);
    surfCurr.selection.named(planeName);
    %surfCurr.set('weakExpression', '0');
    
    surfCurr.set('Js0', ...
        {LL_MODEL.measurements{ss}.Jx, ...
        LL_MODEL.measurements{ss}.Jy, ...
        LL_MODEL.measurements{ss}.Jz});
    surfCurr.name(sprintf('Objective current %i', ss));
    
    measSel = {measSel{:} planeName};
end

measurementSel.set('input', measSel);

%% View!

model.view('view1').set('renderwireframe', false);
hide1 = model.view('view1').hideEntities.create('hide1');
%gg = hide1.geom(3);
%hide1.add([5 10 15 20 25 30 36 41 46]);
%hide1.add([4 9 14 19 24 29 35 40 45]);
hide1.named('pmlSel');



%% Mesh!

fprintf('Meshing\n');

mesh1 = model.mesh.create('mesh1', 'geom1');
sz = mesh1.feature('size');
sz.set('custom', 'on');
%sz.set('hgradactive', 'on');
sz.set('hgrad', '2.5');
%sz.set('hgrad', '2.5');

if ~isempty(LL_MODEL.hmin)
    sz.set('hmin', LL_MODEL.hmin);
end

if ~isempty(LL_MODEL.hmax)
    sz.set('hmax', LL_MODEL.hmax);
end

warning('Ignoring per-object mesh size settings.');
%{
for mm = 1:numMeshes
if ~isempty(LL_MODEL.meshes{mm}.hmax)
    
    szName = sprintf('size%i', mm);
    sz = model.mesh('mesh1').feature.create(szName, 'Size');
    sz.selection.named(sprintf('geom1_%s_bnd', comsolStructureName(mm)));
    sz.set('custom', 'on');
    sz.set('hmaxactive', 'on');
    sz.set('hmax', LL_MODEL.meshes{mm}.hmax);
    
end 
end
%}

% Change mesh size for measurement.
sz = model.mesh('mesh1').feature.create('measSize', 'Size');
sz.selection.geom('geom1', 2);
sz.selection.named('measSel');
sz.set('custom', 'on');
sz.set('hmaxactive', 'on');
sz.set('hmax', 30);

model.mesh('mesh1').feature.create('ftet1', 'FreeTet');

model.mesh('mesh1').run;

model.save([pwd filesep X.MPH]);

%% Study!

fprintf('Study\n')

freqStr = sprintf('c_const/%s[nm]', num2str(3e8/LL_MODEL.frequency));

model.study.create('std1');
model.study('std1').feature.create('freq', 'Frequency');
model.study('std1').feature('freq').set('activate', {'emw' 'on' 'emw2' 'off'});
model.study('std1').feature('freq').set('plist', freqStr);

%model.study.create('std2');
model.study('std1').feature.create('freq1', 'Frequency');
model.study('std1').feature('freq1').set('activate', {'emw' 'off' 'emw2' 'on'});
model.study('std1').feature('freq1').set('plist', freqStr);

fprintf('Solution\n')

model.sol.create('sol1');
model.sol('sol1').study('std1');
model.sol('sol1').attach('std1');

% step 1: forward
model.sol('sol1').feature.create('st1', 'StudyStep');
model.sol('sol1').feature.create('v1', 'Variables');
model.sol('sol1').feature.create('s1', 'Stationary');
%model.sol('sol1').feature('s1').feature.create('p1', 'Parametric');
%model.sol('sol1').feature('s1').feature.create('fc1', 'FullyCoupled');
model.sol('sol1').feature('s1').feature.create('i1', 'Iterative');
model.sol('sol1').feature('s1').feature('i1').feature.create('mg1', 'Multigrid');
model.sol('sol1').feature('s1').feature('i1').feature('mg1').feature('pr').feature.create('sv1', 'SORVector');
model.sol('sol1').feature('s1').feature('i1').feature('mg1').feature('po').feature.create('sv1', 'SORVector');
%model.sol('sol1').feature('s1').feature.remove('fcDef');
% step 2: adjoint
model.sol('sol1').feature.create('st2', 'StudyStep');
model.sol('sol1').feature.create('v2', 'Variables');
model.sol('sol1').feature.create('s2', 'Stationary');
%model.sol('sol1').feature('s2').feature.create('p1', 'Parametric');
%model.sol('sol1').feature('s2').feature.create('fc1', 'FullyCoupled');
model.sol('sol1').feature('s2').feature.create('i1', 'Iterative');
model.sol('sol1').feature('s2').feature('i1').feature.create('mg1', 'Multigrid');
model.sol('sol1').feature('s2').feature('i1').feature('mg1').feature('pr').feature.create('sv1', 'SORVector');
model.sol('sol1').feature('s2').feature('i1').feature('mg1').feature('po').feature.create('sv1', 'SORVector');
%model.sol('sol1').feature('s2').feature.remove('fcDef');

% step 1: compile forward
model.sol('sol1').feature('st1').name('Compile Equations: Frequency Domain');
model.sol('sol1').feature('st1').set('studystep', 'freq');
model.sol('sol1').feature('v1').set('control', 'freq');
model.sol('sol1').feature('s1').set('control', 'freq');
model.sol('sol1').feature('s1').feature('dDef').active(true);
model.sol('sol1').feature('s1').feature('aDef').set('complexfun', true);
%model.sol('sol1').feature('s1').feature('p1').set('plistarr', {freqStr});
%model.sol('sol1').feature('s1').feature('p1').set('pname', {'freq'});
%model.sol('sol1').feature('s1').feature('p1').set('control', 'freq');
model.sol('sol1').feature('s1').feature('i1').set('linsolver', 'bicgstab');
% step 2: compile adjoint
model.sol('sol1').feature('st2').name('Compile Equations: Frequency Domain 2');
model.sol('sol1').feature('st2').set('studystep', 'freq1');
model.sol('sol1').feature('v2').set('control', 'freq1');
model.sol('sol1').feature('s2').set('control', 'freq1');
model.sol('sol1').feature('s2').feature('dDef').active(true);
model.sol('sol1').feature('s2').feature('aDef').set('complexfun', true);
%model.sol('sol1').feature('s2').feature('p1').set('plistarr', {freqStr});
%model.sol('sol1').feature('s2').feature('p1').set('pname', {'freq'});
%model.sol('sol1').feature('s2').feature('p1').set('control', 'freq1');
model.sol('sol1').feature('s2').feature('i1').set('linsolver', 'bicgstab');

% This step, preposterously, seems necessary to make the adjoint solver not
% blow away the forward fields.  This seems very counterintuitive.  :-/
model.study('std1').feature('freq1').set('usesol', 'on');

%% Surface outputs
dsetSurf = model.result.dataset.create('dsetSurfaces', 'Solution');
dsetSurf.name('Surfaces data set');
dsetSurf.selection.geom('geom1', 2);
dsetSurf.selection.all;

measData = model.result.dataset.create('dsetMeas', 'Solution');
measData.name('Objective function data set');
measData.selection.geom('geom1', 2);
measData.selection.named('measSel');

intSurf = model.result.numerical.create('intF', 'IntSurface');
intSurf.name('F');
intSurf.selection.named('measSel');
intSurf.set('probetag', 'none');

tblF = model.result.table.create('tblF', 'Table');
intSurf.set('table', 'tblF');
assert(numMeasurements == 1);
intSurf.set('expr', LL_MODEL.measurements{1}.F);

model.result.export.create('expSurfDF', 'Data');
model.result.export('expSurfDF').name('data set');
model.result.export('expSurfDF').set('data', 'dsetSurfaces');
model.result.export('expSurfDF').set('descr', {''});
model.result.export('expSurfDF').set('filename', 'DF_on_surfaces.txt');
model.result.export('expSurfDF').set('expr', {'mod1.DF'});
model.result.export('expSurfDF').set('resolution', 'custom');
model.result.export('expSurfDF').set('lagorder', '5');

model.result.export.create('expTableF', 'tblF', 'Table');
model.result.export('expTableF').set('filename', 'F.txt');

%%

model.variable.create('var1');
model.variable('var1').model('mod1');
model.variable('var1').set('Ex21', 'up(emw.Ex)-down(emw.Ex)');
model.variable('var1').set('Ey21', 'up(emw.Ey)-down(emw.Ey)');
model.variable('var1').set('Ez21', 'up(emw.Ez)-down(emw.Ez)');
model.variable('var1').set('Hx21', 'up(emw.Hx)-down(emw.Hx)');
model.variable('var1').set('Hy21', 'up(emw.Hy)-down(emw.Hy)');
model.variable('var1').set('Hz21', 'up(emw.Hz)-down(emw.Hz)');
model.variable('var1').set('Dx21', 'up(emw.Dx)-down(emw.Dx)');
model.variable('var1').set('Dy21', 'up(emw.Dy)-down(emw.Dy)');
model.variable('var1').set('Dz21', 'up(emw.Dz)-down(emw.Dz)');
model.variable('var1').set('Bx21', 'up(emw.Bx)-down(emw.Bx)');
model.variable('var1').set('By21', 'up(emw.By)-down(emw.By)');
model.variable('var1').set('Bz21', 'up(emw.Bz)-down(emw.Bz)');
model.variable('var1').set('En21', 'nx*mod1.Ex21+ny*mod1.Ey21+nz*mod1.Ez21');
model.variable('var1').set('Et121', 't1x*mod1.Ex21+t1y*mod1.Ey21+t1z*mod1.Ez21');
model.variable('var1').set('Et221', 't2x*mod1.Ex21+t2y*mod1.Ey21+t2z*mod1.Ez21');
model.variable('var1').set('Hn21', 'nx*mod1.Hx21+ny*mod1.Hy21+nz*mod1.Hz21');
model.variable('var1').set('Ht121', 't1x*mod1.Hx21+t1y*mod1.Hy21+t1z*mod1.Hz21');
model.variable('var1').set('Ht221', 't2x*mod1.Hx21+t2y*mod1.Hy21+t2z*mod1.Hz21');
model.variable('var1').set('Dn21', 'nx*mod1.Dx21+ny*mod1.Dy21+nz*mod1.Dz21');
model.variable('var1').set('Dt121', 't1x*mod1.Dx21+t1y*mod1.Dy21+t1z*mod1.Dz21');
model.variable('var1').set('Dt221', 't2x*mod1.Dx21+t2y*mod1.Dy21+t2z*mod1.Dz21');
model.variable('var1').set('Bn21', 'nx*mod1.Bx21+ny*mod1.By21+nz*mod1.Bz21');
model.variable('var1').set('Bt121', 't1x*mod1.Bx21+t1y*mod1.By21+t1z*mod1.Bz21');
model.variable('var1').set('Bt221', 't2x*mod1.Bx21+t2y*mod1.By21+t2z*mod1.Bz21');
model.variable('var1').set('en', 'nx*emw2.Ex+ny*emw2.Ey+nz*emw2.Ez');
model.variable('var1').set('et1', 't1x*emw2.Ex+t1y*emw2.Ey+t1z*emw2.Ez');
model.variable('var1').set('et2', 't2x*emw2.Ex+t2y*emw2.Ey+t2z*emw2.Ez');
model.variable('var1').set('hn', 'nx*emw2.Hx+ny*emw2.Hy+nz*emw2.Hz');
model.variable('var1').set('ht1', 't1x*emw2.Hx+t1y*emw2.Hy+t1z*emw2.Hz');
model.variable('var1').set('ht2', 't2x*emw2.Hx+t2y*emw2.Hy+t2z*emw2.Hz');
model.variable('var1').set('dn', 'nx*emw2.Dx+ny*emw2.Dy+nz*emw2.Dz');
model.variable('var1').set('dt1', 't1x*emw2.Dx+t1y*emw2.Dy+t1z*emw2.Dz');
model.variable('var1').set('dt2', 't2x*emw2.Dx+t2y*emw2.Dy+t2z*emw2.Dz');
model.variable('var1').set('bn', 'nx*emw2.Bx+ny*emw2.By+nz*emw2.Bz');
model.variable('var1').set('bt1', 't1x*emw2.Bx+t1y*emw2.By+t1z*emw2.Bz');
model.variable('var1').set('bt2', 't2x*emw2.Bx+t2y*emw2.By+t2z*emw2.Bz');
model.variable('var1').set('DF_E', '-emw.omega*imag(-conj(mod1.dn)*conj(mod1.En21)+conj(mod1.et1)*conj(mod1.Dt121)+conj(mod1.et2)*conj(mod1.Dt221))');
model.variable('var1').set('DF_H', 'emw.omega*imag(-conj(mod1.bn)*conj(mod1.Hn21)+conj(mod1.ht1)*conj(mod1.Bt121)+conj(mod1.ht2)*conj(mod1.Bt221))');
model.variable('var1').set('DF', 'mod1.DF_E+mod1.DF_H');

%% Plots

model.result.dataset.create('dset2', 'Solution');
model.result.dataset('dset2').name('Movable meshes');
%model.result.dataset('dset2').selection.geom('geom1', 2);
model.result.dataset('dset2').selection.named('movableMeshes');


model.result.create('pg1', 'PlotGroup3D');
model.result('pg1').name('Electric field');
model.result('pg1').set('frametype', 'spatial');
model.result('pg1').feature.create('mslc1', 'Multislice');
model.result('pg1').feature('mslc1').name('Multislice');

model.result.create('pg2', 'PlotGroup3D');
model.result('pg2').feature.create('mslc1', 'Multislice');
model.result('pg2').name('Electric field 1');
model.result('pg2').set('frametype', 'spatial');
model.result('pg2').feature('mslc1').name('Multislice');
model.result('pg2').feature('mslc1').set('expr', 'emw2.normE');

model.result.create('pg3', 'PlotGroup3D');
model.result('pg3').feature.create('surf1', 'Surface');
model.result('pg3').feature.create('arws1', 'ArrowSurface');
model.result('pg3').feature.create('arws2', 'ArrowSurface');
model.result('pg3').set('data', 'dset2');
model.result('pg3').feature('surf1').set('expr', 'DF');
model.result('pg3').feature('surf1').set('unit', 'W/m^3');
model.result('pg3').feature('surf1').set('descr', '');
model.result('pg3').feature('arws1').set('expr', {'nx*DF*(DF>0)' 'ny*DF*(DF>0)' 'nz*DF*(DF>0)'});
%model.result('pg3').feature('arws1').set('scale', '3.2787942186813154E-10');
model.result('pg3').feature('arws1').set('arrowbase', 'head');
%model.result('pg3').feature('arws1').set('scaleactive', false);
model.result('pg3').feature('arws2').set('expr', {'nx*DF*(DF<0)' 'ny*DF*(DF<0)' 'nz*DF*(DF<0)'});
model.result('pg3').feature('arws2').set('color', 'blue');
%model.result('pg3').feature('arws2').set('scale', '4.240512841916214E-10');
%model.result('pg3').feature('arws2').set('scaleactive', false);

model.result.create('pg4', 'PlotGroup3D');
model.result('pg4').run;
model.result('pg4').feature.create('slc1', 'Slice');
model.result('pg4').feature('slc1').set('quickplane', 'xy');
model.result('pg4').feature('slc1').set('quickzmethod', 'coord');
model.result('pg4').feature('slc1').set('expr', 'emw.Ex');
model.result('pg4').feature('slc1').set('quickz', '10');
%model.result('pg4').run;
%%

fprintf('Saving model.\n');

model.save([pwd filesep X.MPH]);
if X.StopEarly
    model.save([pwd filesep 'saveEarly.mph']);
    warning('quitting early')
    return
end
%model.save([pwd filesep X.MPH]);
%return

%% Run it all
succeeded = 0;
tries = 1;
while ~succeeded
    try
        model.sol('sol1').runAll;
        succeeded = 1;
        if X.SaveFields
            model.save([pwd filesep 'fields_' X.MPH]);
        end
    catch
        warning('Got some crappy problem.  Pausing and retrying.\n');
        pause(10);
        tries = tries + 1;
    end
end

if tries > 1
    fprintf('Succeded with %i tries\n', tries);
end

model.result.table('tblF').clearTableData;
model.result.numerical('intF').set('table', 'tblF');
model.result.numerical('intF').setResult;
model.result.export('expSurfDF').run;
model.result.export('expTableF').run;

%% Plots!!

figure(1); clf
mphplot(model, 'pg1', 'rangenum', 1);

figure(2); clf
mphplot(model, 'pg2', 'rangenum', 1);

figure(3); clf
mphplot(model, 'pg3', 'rangenum', 1);

figure(4); clf
mphplot(model, 'pg4', 'rangenum', 1);

end

function onOff = onOffString(boolval)

if boolval
    onOff = 'on';
else
    onOff = 'off';
end
end

function addRect(geom, rectBounds, impName)

fprintf('Making %s\n', impName);

rectCorner = rectBounds(1:3);
rectSize = rectBounds(4:6) - rectBounds(1:3);

cornerCell = rect2cell(rectCorner);
sizeCell = rect2cell(rectSize);

b = geom.feature.create(impName, 'Block');
b.name(impName);
b.set('pos', cornerCell);
b.set('size', sizeCell);
b.set('createselection', true);

end

function addImport(geom, v, f, impFileName, impName)

fprintf('Making %s\n', impName);
fname = [impFileName '.stl'];
fpath = [pwd filesep 'importMeshes' filesep fname];

ll.writeSTL(v, f, fpath);

%importSTL = geom.feature.create([impName '_imp'], 'Import');
importSTL = geom.feature.create(impName, 'Import');

importSTL.set('neighangle', '2.0');
importSTL.set('filename', fpath);
importSTL.set('facecurv', '10');
importSTL.set('type', 'stlvrml');
importSTL.set('curvedface', 'manual');
importSTL.set('facepartition', 'manual');
importSTL.set('createselection', true);
importSTL.set('faceangle', '1');
importSTL.set('facecleanup', '0');
%model.geom('geom1').feature('import_1').set('faceangle', '20');

%makeSolid = geom.feature.create(impName, 'ConvertToSolid');
%makeSolid.selection('input').named([impName '_imp']);
%makeSolid.set('createselection', true);

end

function addImportDifferences(geom, mm, listOfDifferences)

if ~isempty(listOfDifferences)

    subtractStructures = arrayfun(@comsolImportName,...
        listOfDifferences, 'UniformOutput', false);

    d = geom.feature.create(comsolStructureName(mm), 'Difference');
    d.name(sprintf('Disjoint chunk %i', mm));
    d.set('keep', true);
    d.selection('input').set(comsolImportName(mm));
    d.selection('input2').set(subtractStructures);
    d.set('createselection', true);
else
    % do nothing...
    d = geom.feature.create(comsolStructureName(mm), 'Union');
    d.name(sprintf('Original chunk %i', mm));
    d.set('keep', true);
    d.selection('input').set(comsolImportName(mm));
    d.set('createselection', true);
end

end

function addRawPML(geom, rectBounds, pmlIndex)

pmlCorner = rectBounds(1:3);
pmlSize = rectBounds(4:6) - rectBounds(1:3);

pmlCornerCell = rect2cell(pmlCorner);
pmlSizeCell = rect2cell(pmlSize);

pmlName = comsolPMLName(pmlIndex);
b = geom.feature.create(pmlName, 'Block');
b.name(pmlName);
b.set('pos', pmlCornerCell);
b.set('size', pmlSizeCell);
b.set('createselection', true);

end

function addNonPML(geom, bounds)

nonPMLPos = bounds(1:3);
nonPMLSize = bounds(4:6) - bounds(1:3);

b = geom.feature.create(comsolNonPMLName(), 'Block');
b.name('Non-PML');
b.set('pos', rect2cell(nonPMLPos));
b.set('size', rect2cell(nonPMLSize));


end

function itDoes = rectContains(verts, bounds)

p1 = repmat(bounds(1:3), [size(verts, 1) 1]);
p2 = repmat(bounds(4:6), [size(verts, 1) 1]);

itDoes = all(all(p1 <= verts & p2 >= verts));

end

function addNonPMLStructure(geom, mm)

intName = comsolNonPMLPartName(mm);

ii = geom.feature.create(intName, 'Intersection');
ii.name(intName);
ii.set('createselection', true);
ii.set('keep', true);
ii.selection('input').set({comsolNonPMLName(), comsolStructureName(mm)});

end

function copyNonPMLStructure(geom, mm)

intName = comsolNonPMLPartName(mm);

d = geom.feature.create(intName, 'Union');
d.name(sprintf('Original non-PML chunk %i', mm));
d.set('keep', true);
d.selection('input').set(comsolStructureName(mm));
d.set('createselection', true);

%{
ii = geom.feature.create(intName, 'Intersection');
ii.name(intName);
ii.set('createselection', true);
ii.set('keep', true);
ii.selection('input').set({comsolNonPMLName(), comsolStructureName(mm)});
%}

end

function addPMLIntersection(geom, ss, pmlIndex, pmlChunkIndex)

structureName = comsolStructureName(ss);
pmlName = comsolPMLName(pmlIndex);
pmlBitName = comsolPMLChunkName(pmlChunkIndex);

ii = geom.feature.create(pmlBitName, 'Intersection');
ii.name(sprintf('PML %s', pmlBitName));
ii.set('keep', true);
ii.selection('input').set({pmlName, structureName});
ii.set('createselection', true);

end

function deleteImports(geom, numImports)

d = geom.feature.create('deleteImports', 'Delete');
d.name('Delete imports');
d.selection('input').init;
importStructures = arrayfun(@comsolImportName, 1:numImports, ...
    'UniformOutput', false);
d.selection('input').set(importStructures);

end

function deleteDisjointStructures(geom, disjointStructureIndices)

d = geom.feature.create('deleteDisjointImports', 'Delete');
d.name('Delete disjoint structures');
d.selection('input').init;
disjointStructures = arrayfun(@comsolStructureName, disjointStructureIndices, ...
    'UniformOutput', false);
d.selection('input').set(disjointStructures);

end

function deleteNonPML(geom)

d = geom.feature.create('deleteNonPML', 'Delete');
d.name('Delete non-PML');
d.selection('input').init;
d.selection('input').set(comsolNonPMLName());

end

function deleteRawPML(geom, numPMLs)

d = geom.feature.create('deletePMLInputs', 'Delete');
d.name('Delete PML inputs');
d.selection('input').init;
d.selection('input').set(arrayfun(@comsolPMLName, 1:numPMLs, ...
    'UniformOutput', false));

end

function domainNumber = identifyByBounds(model, bbox)

domainNumber = selectBoxHelper(model, bbox + [-1 -1 -1 1 1 1]);
if numel(domainNumber) == 1
    return
end

% Now start trying to trim off everything that's not in those bounds

pushIn = diag([1 1 1 -1 -1 -1]);

for ii = 1:6
    removeThese = selectBoxHelper(model, bbox + pushIn(ii,:));
    domainNumber = setdiff(domainNumber, removeThese);
    if numel(domainNumber) < 2
        return;
    end
end

warning('There are still %i matching domains!', numel(domainNumber));

end

function selections = selectBoxHelper(model, bbox)

indexer = [1 4; 2 5; 3 6];

selections = mphselectbox(model, 'geom1', bbox(indexer), 'domain');

end

%% 

%% Utility function

function r = rect2cell(A)
r = arrayfun(@(aa) sprintf('%i', aa), A, 'UniformOutput', false);
end

%% Strings needed for COMSOL

function nom = comsolImportName(mm)
nom = sprintf('import_%i', mm);
end

function nom = comsolStructureName(mm)
nom = sprintf('structure_%i', mm);
end

function nom = comsolPMLName(pp)
nom = sprintf('blockPML_%i', pp);
end

function nom = comsolNonPMLName()
nom = 'blockNonPML';
end

function nom = comsolNonPMLPartName(mm)
nom = sprintf('intStructure_%i', mm);
end

function nom = comsolPMLChunkName(cc)
nom = sprintf('intPML_%i', cc);
end


