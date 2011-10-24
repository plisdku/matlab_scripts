function [cc, movableVertices] = readDeps()

%fclose all
fid = fopen('Depsilon', 'r');
if fid == -1
    error('Cannot open Depsilon');
end

cc = cell(0,3);

for tensor_ij = 0:2
    
    % Read TENS marker and (i,j)
    str = fread(fid, 4, 'char=>char')';
    
    if ~strcmp(str, 'TENS')
        error('Expected TENS marker');
    end
    
    ij = fread(fid, 2, 'uint32')';
    if any(ij ~= [tensor_ij tensor_ij])
        error('Expected tensor direction %i but got %i', tensor_ij, ij);
    end
    
    % Read vertices
    numVertices = fread(fid, 1, 'uint32');
    %fprintf('%i vertices coming up.\n', numVertices);
    
    done = 0;
    while ~done
        
        str = fread(fid, 4, 'char=>char')';
        if strcmp(str, 'TEND')
            done = 1;
            continue
        elseif ~strcmp(str, 'VERT')
            error('Expected VERT marker');
        end
        %if ~strcmp(str, 'VERT')
        %    error('Expected VER marker');
        %end
        
        vNum = fread(fid, 1, 'uint32');
        %if vNum ~= vert
        %    error('Expected vertex %i but got %i', vert, vNum);
        %end
        
        %fprintf('\tVertex %i\n', vNum);
        
        for freeDir = 0:2
            
            % For each free direction:
            str = fread(fid, 4, 'char=>char')';
            if ~strcmp(str, 'DIRE')
                error('Expected DIRE marker');
            end

            xyz = fread(fid, 1, 'uint32');
            if xyz ~= freeDir
                error('Expected vertex direction %i but got %i', ...
                    freeDir, xyz);
            end
            
            %fprintf('\t\tDirection %i\n', xyz);
            
            % Read numerator
            
            numerDone = 0;
            while ~numerDone
                str = fread(fid, 4, 'char=>char')';
                if strcmp(str, 'NEND')
                    numerDone = 1;
                    continue;
                elseif ~strcmp(str, 'NUME')
                    error('Expected NUME marker');
                end
                
                lag = fread(fid, 1, 'uint32');
                %if lag ~= nn
                %    error('Expected lag %i, got %i', nn, lag);
                %end
                %fprintf('\t\t\tLag %i\n', lag);
                
                numValues = fread(fid, 1, 'uint32');
                % Get the values
                vals = fread(fid, numValues, 'float64');
                
                % Get the corresponding indices
                inds = fread(fid, numValues, 'uint32');
                
                 %fprintf('\t\t\tIndices: ');
                 %fprintf('%i ', inds);
                 %fprintf('\n\t\t\tValues: ');
                 %fprintf('%2.4f ', vals);
                 %fprintf('\n');
                 
	cc{vNum+1, freeDir+1}.tensor{tensor_ij+1, tensor_ij+1}.DB{lag+1}.coefficients ...
        = vals;
    cc{vNum+1, freeDir+1}.tensor{tensor_ij+1, tensor_ij+1}.DB{lag+1}.indices...
        = inds;
            end
            
            % Read denominator
            
            denomDone = 0;
            while ~denomDone
                str = fread(fid, 4, 'char=>char')';
                if strcmp(str, 'DEND')
                    denomDone = 1;
                    continue;
                elseif ~strcmp(str, 'DENO');
                    error('Expected DENO marker');
                end
                
                lag = fread(fid, 1, 'uint32');
                %if lag ~= nn
                %    error('Expected lag %i, got %i', nn, lag);
                %end
                %fprintf('\t\t\tLag %i\n', nn);
                
                numValues = fread(fid, 1, 'uint32');
                % Get the values
                vals = fread(fid, numValues, 'float64');
                
                % Get the corresponding indices
                inds = fread(fid, numValues, 'uint32');
                
                 %fprintf('\t\t\tIndices: ');
                 %fprintf('%i ', inds);
                 %fprintf('\n\t\t\tValues: ');
                 %fprintf('%2.4f ', vals);
                 %fprintf('\n');
                 
	cc{vNum+1, freeDir+1}.tensor{tensor_ij+1, tensor_ij+1}.EH{lag+1}.coefficients ...
        = vals;
    cc{vNum+1, freeDir+1}.tensor{tensor_ij+1, tensor_ij+1}.EH{lag+1}.indices...
        = inds;
            end
            
        end
    end
        
    %fprintf('Done with %i\n', tensor_ij);
    
end


% Figure out which vertices have data in them.
movableVertices = [];
if exist('cc')
for nn = 1:size(cc,1)
   movable = 0;
   for xyz = 1:size(cc,2)
       if ~isempty(cc{nn,xyz})
           movable = 1;
       end
   end
   if movable
       movableVertices = [movableVertices, nn];
   end
end
end
