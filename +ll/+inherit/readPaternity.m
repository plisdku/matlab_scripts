outFile = '/private/tmp/nefOut736429.6919944.txt';

indicesOld = [];
indicesNew = [];

[fh,message] = fopen(outFile);
assert(fh ~= -1);

iOld = 0;
while 1
    tline = fgetl(fh);
    if ~ischar(tline)
        break
    end
    
    [indices, count, errmsg, nextindex] = sscanf(tline, '%*d: %d', [1,inf]); 
    
    indicesOld = [indicesOld; iOld*ones(length(indices),1)];
    indicesNew = [indicesNew; reshape(indices, [], 1)];
    
    iOld = iOld + 1;
    
end

fclose(fh)