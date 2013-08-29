function writeMultiOFF(fid_or_fname, vertices, faces)

if ischar(fid_or_fname)
    fh = fopen(fid_or_fname, 'w');
else
    fh = fid_or_fname;
end

shellFaces = divideIntoShells(faces);

numShells = numel(shellFaces);
shellVolumes = zeros(numShells,1);

for ss = 1:numShells
    shellVolumes(ss) = polyhedronVolume(vertices, shellFaces{ss});
end

numPositiveShells = sum(shellVolumes > 0);
numNegativeShells = sum(shellVolumes < 0);
assert(numPositiveShells + numNegativeShells == numShells);

% Write how many positive and negative shells I've got.
fprintf(fh, '%i %i\n', numPositiveShells, numNegativeShells);

% Write the positive shells
for ss = 1:numShells
if shellVolumes(ss) > 0
    writeOneOFF(fh, vertices, shellFaces{ss});
end
end

for ss = 1:numShells
if shellVolumes(ss) < 0
    writeOneOFF(fh, vertices, shellFaces{ss});
end
end

if ischar(fid_or_fname)
    fclose(fh);
end

function writeOneOFF(fh, vertices, faces)

fprintf(fh, 'OFF\n');

numVertices = size(vertices, 1);
numFaces = size(faces, 1);
numEdges = countEdges(vertices, faces);

fprintf(fh, '%i %i %i\n', numVertices, numFaces, numEdges);

fprintf(fh, '%2.8f %2.8f %2.8f\n', vertices');
fprintf(fh, '3 %i %i %i\n', faces' - 1);



% Sometimes they say the number of edges can be ignored, so I'll try
% that...
function n = countEdges(vertices, faces)

n = 0;