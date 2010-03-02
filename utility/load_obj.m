function [vertices, objects] = load_obj(filename);
% load_obj  Limited loading of Alias Wavefront OBJ files
%   [vertices, objects] = load_obj(filename) will read an OBJ file and
%   return the vertices and a cell array of groups from that file, suitable
%   for plotting as a patch.
%
%   To inspect the file contents,
%      [v, obj] = load_obj('myfile.obj');
%       obj
%   will show how many objects are in the file; they may be accessed using
%   cell notation, e.g
%      obj{5}.name
%   will tell the name of the current object.  Trogdor will group materials
%   into OBJ groups, so the name will likely be 'Gold' or similar.
%
%   The data may be plotted as
%
%   patch('Vertices', verts, 'Faces', objects{1}.faces)
%   
%   and additional plotting arguments may follow, such as the pairs
%       'FaceColor', 'g'     or 'FaceColor', [0 1 0]
%       'EdgeColor', 'none'
%       'FaceAlpha', 0.5
%   et cetera.
%
%   See also: trogRect, trogEllipsoid, patch

fid = fopen(filename, 'r');

done = 0;

groupName = 'Default';
materialName = 'DefaultMaterial';
vertices = [];
normals = [];
faceverts = [];
facenorms = [];

objects = [];
numObjects = 0;

while (done == 0)
    ll = fgets(fid);
    if ~ischar(ll)
        done = 1;
    else
       [token, remainder] = strtok(ll);
       
       switch token
           case -1
           case '#'
               
           case 'v'
               [dat, count] = sscanf(remainder, '%f');
               if count == 3
                   vertices = [vertices; dat'];
               else
                   error('Cannot parse line %s', ll);
               end
               
           case 'vt'
               [dat, count] = sscanf(remainder, '%f')
               if count == 3
                   normals = [normals; dat'];
               else
                   error('Cannot parse line %s', ll);
               end
               
           case 'f'
               [dat2, count2] = sscanf(remainder, '%f//%f');
               
               if count2 == 8
                   faceverts = [faceverts; dat2([1 3 5 7])'];
                   facenorms = [facenorms; dat2([2 4 6 8])'];
               else
                   error('Cannot parse line %s', ll);
               end
               
           case 'g'
               [name, count] = sscanf(remainder, '%s');
               
               if count
                   if numObjects ~= 0
                       objects{numObjects}.faces = faceverts;
                       %objects{numObjects}.normals = facenorms;
                       faceverts = [];
                       facenorms = [];
                   end
                   numObjects = numObjects+1;
                   objects{numObjects}.name = name;
               else
                   error('Cannot parse line %s', ll);
               end
               
           case 'usemtl'
               disp(sprintf('Ignoring material directive %s', remainder));
       end
    end
end

if numObjects ~= 0
   objects{numObjects}.faces = faceverts;
   %objects{numObjects}.normals = facenorms;
   faceverts = [];
   facenorms = [];
end
    

% patch('Vertices', verts, 'Faces', obj{2}.vertices, 'FaceColor', 'g',
% 'EdgeColor', 'None')




