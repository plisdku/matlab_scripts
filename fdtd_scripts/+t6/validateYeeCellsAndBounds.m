function valid = validateYeeCellsAndBounds(X)

valid = 1;

if ~isempty(X.YeeCells) && ~isempty(X.Bounds)
    error('YeeCells and Bounds are mutually exclusive options');
end

if isempty(X.YeeCells) && isempty(X.Bounds)
    error('Either YeeCells or Bounds must be given');
end

if isempty(X.Bounds) && size(X.YeeCells, 2) ~= 6
    error('YeeCells must have six columns.');
end

if isempty(X.YeeCells) && size(X.Bounds, 2) ~= 6
    error('Bounds must have six columns.');
end

