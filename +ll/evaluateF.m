function F = evaluateF

fid = fopen('F.txt');
AA = cell2mat(textscan(fid, '%n', 'CommentStyle', '%'));
fclose(fid);

F = AA(2);

