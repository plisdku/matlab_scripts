% writeSimulation2   Generate description files for Trogdor simulation
%   writeSimulation(simulationObject) takes a simulation from createSimulation
%   and generates the file 'params.xml' and other required files (HeightMap and
%   KeyImage images, Input files) for a simulation.  The simulation can be 
%   run with Trogdor.
%
%   writeSimulation(simulationObject, paramFileName) lets you specify a param
%   file other than 'params.xml'.
%
%   Example:
%
%   s = createSimulation([dx dy dz dt], 100, g);
%   writeSimulation(s, 'myParamFile.xml');
%   system('trogdor > trogdorOutput.txt');
%
%   See also: createSimulation, createGrid
%
%   version 4.5
%   July 29, 2008
function writeSimulationFile(simObject, varargin)

if nargin > 1
    paramFileName = varargin{1};
else
    paramFileName = 'params.xml';
end
saveTwister = rand('twister');
rand('twister', 5);
%assert(strcmp(simObject.type, 'Simulation'));
%{
numT = 1000;
dx = 10e-9;
dt = 0.99 * dx / (sqrt(3.0) * 2.99e8);

bb = createBlock([-60 -60 0 80 80 0], 'Air');
b = createBlock([-10 -10 0 15 2 0], 'Gold');
s = createEllipsoid([40 40 0 70 70 0], 'PEC');
a = createAssembly(bb, b, s);
mAu = createMaterial('Gold', 'DrudeMetalModel', { {'epsinf', '12.9898'}, ...
    {'omegap', '4.0217e15'}, {'tauc', '9.0441e-15'} });
mPEC = createMaterial('PEC', 'PECModel');
mAir = createMaterial('Air', 'StaticDielectricModel', { {'epsr', '1.0'}, ...
    {'mur', '1.0'} } );
o = createOneFieldOutput('ez', [-65 -65 0 85 85 0], 'outEz');

ll = createLink('SrcGrid', [-50 0 0 70 0 0], [-50 -50 0 70 70 0]);

g = createGrid('MainGrid', [10 10 0 10 10 0], a, mAu, mAir, mPEC, o, ll);

tfB = createBlock([-60 0 0 80 0 0], 'Air');
oo = createOneFieldOutput('ez', [-60 0 0 80 0 0], 'srcEz');
src = createSource('ez', [-55 0 0 -55 0 0], sin(2*pi*0.01*(1:numT)));
assemblery = createAssembly(tfB);
g2 = createGrid('SrcGrid', [10 0 0 10 0 0], assemblery, src, oo, mAir);



simObject = createSimulation([dx dx dx dt], numT, g, g2);

%}

%% Calculate all actual extents (grid cells)
% This entails walking through the dimensions and finding a center cell,
% then displacing everything accordingly on writing to file.
%
% Then save all the extents for all the children.


for gnum = 1:length(simObject.grids)
    extents = [inf, inf, inf, -inf, -inf, -inf];
    theGrid = simObject.grids{gnum};
    
    for cnum = 1:length(theGrid.pieces)
        theChild = theGrid.pieces{cnum};
        if strcmp(theChild.type, 'Assembly')
            for piece = theChild.children
                if length(piece{1}.dimensions) ~= 6
                    warning('Piece has wrong dimensions.');
                    piece{1}
                    error('Giving up.');
                end
                
                for (nn = 1:3)
                    extents(nn) = min(extents(nn), piece{1}.dimensions(nn));
                end
                for (nn = 4:6)
                    extents(nn) = max(extents(nn), piece{1}.dimensions(nn));
                end
                
            end
        end
    end
    
    origin = [extents(1:3) - theGrid.PML(1:3), extents(1:3) - theGrid.PML(1:3)];
    theGrid.origin = origin;
    theGrid.numCells = extents(4:6)-extents(1:3) + [1 1 1] + ...
        theGrid.PML(1:3) + theGrid.PML(4:6);
    theGrid.roi = [theGrid.PML(1:3), theGrid.numCells - [1 1 1] ...
        - theGrid.PML(4:6)];
    
    % write back to the original grid object since Matlab doesn't have
    % pointers
    simObject.grids{gnum} = theGrid;
end

%% Find the source origins of the link regions.

for gnum = 1:length(simObject.grids)
    theGrid = simObject.grids{gnum};
    
    for cnum = 1:length(theGrid.pieces)
        theChild = theGrid.pieces{cnum};
        if strcmp(theChild.type, 'Link')
            srcGridName = theChild.src;
            for sGrid = simObject.grids
                if strcmp(sGrid{1}.name, srcGridName)
                    theChild.srcOrigin = sGrid{1}.origin;
                end
            end
            theGrid.pieces{cnum} = theChild;
        end
    end
    % write back to the original grid object since Matlab doesn't have
    % pointers
    simObject.grids{gnum} = theGrid;
end

%% Traverse the simulation and generate the XML object

docNode = com.mathworks.xml.XMLUtils.createDocument('Simulation');
docRootNode = docNode.getDocumentElement;
docRootNode.setAttribute('dx', num2str(simObject.dx));
docRootNode.setAttribute('dy', num2str(simObject.dy));
docRootNode.setAttribute('dz', num2str(simObject.dz));
docRootNode.setAttribute('dt', num2str(simObject.dt));
docRootNode.setAttribute('numT', num2str(simObject.numT));


for theGrid = simObject.grids
    gridXML = docNode.createElement('Grid');
    gridXML.setAttribute('name', theGrid{1}.name);
    gridXML.setAttribute('nx', num2str(theGrid{1}.numCells(1)));
    gridXML.setAttribute('ny', num2str(theGrid{1}.numCells(2)));
    gridXML.setAttribute('nz', num2str(theGrid{1}.numCells(3)));
    gridXML.setAttribute('regionOfInterest', sprintf('%i ', theGrid{1}.roi));
    %gridXML.setAttribute('dumpGrid', '1');
    
    for thePiece = theGrid{1}.pieces
    switch thePiece{1}.type
        case 'Assembly'
            asXML = docNode.createElement('Assembly');
            for as = thePiece{1}.children
            switch as{1}.type
                case 'Block'
                    blockXML = docNode.createElement('Block');
                    blockXML.setAttribute('fillRect', sprintf('%i ', ...
                        as{1}.dimensions - theGrid{1}.origin));
                    blockXML.setAttribute('material', as{1}.material);
                    blockXML.setAttribute('fillStyle', as{1}.fillStyle);
                    asXML.appendChild(blockXML);
                case 'Ellipsoid'
                    sphXML = docNode.createElement('Ellipsoid');
                    
                    sphXML.setAttribute('fillRect', sprintf('%i ', ...
                        as{1}.dimensions - theGrid{1}.origin));
                    sphXML.setAttribute('material', as{1}.material);
                    sphXML.setAttribute('fillStyle', as{1}.fillStyle);
                    asXML.appendChild(sphXML);
                case 'HeightMap'
                    heightXML = docNode.createElement('HeightMap');
                    heightXML.setAttribute('fillRect', sprintf('%i ', ...
                        as{1}.dimensions - theGrid{1}.origin));
                    heightXML.setAttribute('material', as{1}.material);
                    heightXML.setAttribute('fillStyle', as{1}.fillStyle);
                    rowStr = makeAxisString(as{1}.row);
                    colStr = makeAxisString(as{1}.col);
                    upStr = makeAxisString(as{1}.up);
                    heightXML.setAttribute('row', rowStr);
                    heightXML.setAttribute('column', colStr);
                    heightXML.setAttribute('up', upStr);
                    
                    % now write the height map
                    imfilename = char('a' + randint(1, 4, [0, 25]));
                    imfilename = ['height_', imfilename, '.bmp'];
                    heightXML.setAttribute('file', imfilename);
                    
                    if isa(as{1}.image, 'double')
                        theImg = uint8(round(as{1}.image*255));
                    else
                        theImg = as{1}.image;
                    end
                    
                    if  max(theImg(:)) > 255
                        warning('Heightmap max intensity should be 255 or 1.0');
                    elseif min(theImg(:)) < 0
                        warning('Heightmap min intensity should be 0 or 0.0');
                    end
                    
                    imwrite(theImg, imfilename, 'bmp');
                    
                    asXML.appendChild(heightXML);
                case 'KeyImage'
                    keyXML = docNode.createElement('KeyImage');
                    keyXML.setAttribute('fillRect', sprintf('%i ', ...
                        as{1}.dimensions - theGrid{1}.origin));
                    rowStr = makeAxisString(as{1}.row);
                    colStr = makeAxisString(as{1}.col);
                    keyXML.setAttribute('row', rowStr);
                    keyXML.setAttribute('column', colStr);
                    
                    % now write the height map
                    imfilename = char('a' + randint(1, 4, [0, 25]));
                    imfilename = ['key_', imfilename, '.bmp'];
                    keyXML.setAttribute('file', imfilename);
                    
                    if isa(as{1}.image, 'double')
                        theImg = uint8(round(as{1}.image*255));
                    else
                        theImg = as{1}.image;
                    end
                    
                    if ndims(as{1}.image) == 3 % a color image
                        for nn = 1:length(as{1}.materials)
                            keyTag = docNode.createElement('Tag');
                            tag = as{1}.tags{nn};
                            if isa(tag, 'double')
                                tag = uint8(round(tag*255));
                            end
                            assert(length(tag) == 3);
                            hexred = dec2hex(tag(1), 2);
                            hexgreen = dec2hex(tag(2), 2);
                            hexblue = dec2hex(tag(3), 2);
                            
                            keyTag.setAttribute('color', ...
                                ['#', hexred, hexgreen, hexblue]);
                            keyTag.setAttribute('material', as{1}.materials{nn});
                            keyTag.setAttribute('fillStyle', as{1}.fillstyle{nn});
                            keyXML.appendChild(keyTag);
                        end
                    elseif ndims(as{1}.image) == 2 % a greyscale image
                        
                        as{1}.image = repmat(as{1}.image, [1 1 3]);
                        
                        for nn = 1:length(as{1}.materials)
                            keyTag = docNode.createElement('Tag');
                            tag = as{1}.tags{nn};
                            if isa(tag, 'double')
                                tag = uint8(round(tag*255));
                            end
                            assert(length(tag) == 1);
                            hexred = dec2hex(tag(1), 2);
                            hexgreen = dec2hex(tag(1), 2);
                            hexblue = dec2hex(tag(1), 2);
                            
                            keyTag.setAttribute('color', ...
                                ['#', hexred, hexgreen, hexblue]);
                            keyTag.setAttribute('material', as{1}.materials{nn});
                            keyTag.setAttribute('fillStyle', as{1}.fillstyle{nn});
                            keyXML.appendChild(keyTag);
                        end
                    else
                        error(sprintf('Key image has %i dimensions.', ...
                            ndims(as{1}.image)));
                    end
                    
                    asXML.appendChild(keyXML);
                    
                    % now write the key image
                    imwrite(theImg, imfilename, 'bmp');
            end
            end
            gridXML.appendChild(asXML);
        case 'Material'
            matXML = docNode.createElement('Material');
            
            matXML.setAttribute('name', thePiece{1}.name);
            matXML.setAttribute('class', thePiece{1}.class);
            
            if length(thePiece{1}.parameters) ~= 0
                parXML = docNode.createElement('Params');
                for nn = 1:length(thePiece{1}.parameters)
                    par = thePiece{1}.parameters{nn};
                    parXML.setAttribute(par{1}, par{2});
                end
                matXML.appendChild(parXML);
            end
            
            gridXML.appendChild(matXML);
        case 'Output'
            outXML = docNode.createElement('Output');
            outXML.setAttribute('class', thePiece{1}.class);
            outXML.setAttribute('period', num2str(thePiece{1}.period));
            outXML.setAttribute('filePrefix', thePiece{1}.filename);
            
            parXML = docNode.createElement('Params');
            parXML.setAttribute('field', thePiece{1}.field);
            parXML.setAttribute('region', sprintf('%i ', ...
                thePiece{1}.dimensions - theGrid{1}.origin));
            parXML.setAttribute('stride', sprintf('%i ', ...
                thePiece{1}.stride));
            outXML.appendChild(parXML);
            
            gridXML.appendChild(outXML);
        case 'Input'
            inXML = docNode.createElement('Input');
            inXML.setAttribute('class', thePiece{1}.class);
            inXML.setAttribute('filePrefix', thePiece{1}.filename);
            
            parXML = docNode.createElement('Params');
            parXML.setAttribute('field', thePiece{1}.field);
            parXML.setAttribute('region', sprintf('%i ', ...
                thePiece{1}.dimensions - theGrid{1}.origin));
            outXML.appendChild(parXML);
            gridXML.appendChild(outXML);
        case 'Source'
            inXML = docNode.createElement('Source');
            inXML.setAttribute('region', sprintf('%i ', ...
                thePiece{1}.dimensions - theGrid{1}.origin));
            inXML.setAttribute('field', thePiece{1}.fieldName);
            inXML.setAttribute('polarization', sprintf('%i ', ...
                thePiece{1}.polarization));
            
            parXML = docNode.createElement('Params');
            if (thePiece{1}.isSoft == 1)
                parXML.setAttribute('soft', '1');
            end
            
            if isfield(thePiece{1}, 'data')
                infile = char('a' + randint(1, 4, [0, 25]));
                infile = ['src_', infile];
                inXML.setAttribute('file', infile);
            
                % write the output file
                data = thePiece{1}.data;
                dsize = size(data);
                dsize = dsize(1:end-1);
                if length(dsize) == 1
                    dsize = [dsize(1), 1, 1];
                elseif length(dsize) == 2
                    dsize = [dsize(1:2), 1];
                end
                
                outFH = fopen(infile, 'w');
                try
                    count = fwrite(outFH, data(:), 'float32');
                catch
                    warning('Could not write source data file.');
                end
                fclose(outFH);
            
            elseif isfield(thePiece{1}, 'formula')
                inXML.setAttribute('formula', thePiece{1}.formula);
            end
            
            inXML.appendChild(parXML);
            
            gridXML.appendChild(inXML);
            
        case 'TFSFSource'
            tfXML = docNode.createElement('TFSFSource');
            tfXML.setAttribute('TFRect', sprintf('%i ', ...
                thePiece{1}.dimensions - theGrid{1}.origin));
            tfXML.setAttribute('class', thePiece{1}.class);
            tfXML.setAttribute('direction', sprintf('%i ', ...
                thePiece{1}.direction));
            
            parXML = docNode.createElement('Params');
            parXML.setAttribute('polarization', sprintf('%i ', ...
                thePiece{1}.polarization));
            
            if (isfield(thePiece{1}, 'formula'))
                parXML.setAttribute('formula', thePiece{1}.formula);
            elseif (isfield(thePiece{1}, 'data'))
                infile = char('a' + randint(1, 4, [0, 25]));
                infile = ['tfsfsrc_', infile];
                parXML.setAttribute('filename', infile);
            
                % write the output file
                data = thePiece{1}.data;
                dsize = size(data);
                dsize = dsize(1:end-1);
                if length(dsize) == 1
                    dsize = [dsize(1), 1, 1];
                elseif length(dsize) == 2
                    dsize = [dsize(1:2), 1];
                end
                
                outFH = fopen(infile, 'w');
                try
                    count = fwrite(outFH, data(:), 'float32');
                catch
                    warning('Could not write TFSF source data file.');
                end
                fclose(outFH);
            end
            
            tfXML.appendChild(parXML);
            
            if length(thePiece{1}.omitSides) ~= 0
                for nn = 1:length(thePiece{1}.omitSides)
                    omitXML = docNode.createElement('OmitSide');
                    omitXML.appendChild(docNode.createTextNode(sprintf('%i ', ...
                        thePiece{1}.omitSides{nn})));
                    tfXML.appendChild(omitXML);
                end
            end
            
            gridXML.appendChild(tfXML);
            
        case 'Link'
            linkXML = docNode.createElement('Link');
            linkXML.setAttribute('type', 'TF');  % why bother with anything else
            linkXML.setAttribute('sourceRect', sprintf('%i ', ...
                thePiece{1}.srcDimensions - thePiece{1}.srcOrigin));
            linkXML.setAttribute('destRect', sprintf('%i ', ...
                thePiece{1}.destDimensions - theGrid{1}.origin));
            linkXML.setAttribute('sourceGrid', thePiece{1}.src);
            
            if length(thePiece{1}.omitSides) ~= 0
                for nn = 1:length(thePiece{1}.omitSides)
                    omitXML = docNode.createElement('OmitSide');
                    omitXML.appendChild(docNode.createTextNode(sprintf('%i ', ...
                        thePiece{1}.omitSides{nn})));
                    linkXML.appendChild(omitXML);
                end
            end
            gridXML.appendChild(linkXML);
        otherwise
    end
    end
    
    
    docRootNode.appendChild(gridXML);
end
rand('twister', saveTwister);
%xmlFileName = 'params.xml';
%xmlFileName = [tempname, '.xml'];
xmlwrite(paramFileName, docNode);




