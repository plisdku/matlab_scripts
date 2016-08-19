function boundaries = outerDomainBoundaryEntities(model, domains)
% boundaries = outerDomainBoundaryEntities(model, domains)
%
% Return entity indices of boundaries neighboring only one of the input
% domains, i.e. outer boundaries of the union of domains.

% Collect a list of all boundary entity indices.
% Boundaries will appear in the list once for each domain they bound,
% so zero times, once, or twice depending on topology.
boundaries_redundant = [];

for dd = domains
    domainBoundaries = mphgetadj(model, 'geom1', 'boundary', 'domain', dd);
    
    n0 = length(boundaries_redundant)+1;
    n1 = n0 + length(domainBoundaries)-1;
    boundaries_redundant(n0:n1) = domainBoundaries;
end

% Now return just the elements that occur once.

A = sparse(boundaries_redundant, 1:length(boundaries_redundant), ...
    ones(size(boundaries_redundant)));
occurrences = sum(A, 2);  % occurrences(boundary_id) = number of occurrences

boundaries = find(occurrences == 1);


