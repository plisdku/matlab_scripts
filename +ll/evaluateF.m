function F = evaluateF(varargin)

X.Column = 2;
X = parseargs(X, varargin{:});

fid = fopen('F.txt');
AA = cell2mat(textscan(fid, '%n', 'CommentStyle', '%'));
fclose(fid);

F = AA(X.Column);

