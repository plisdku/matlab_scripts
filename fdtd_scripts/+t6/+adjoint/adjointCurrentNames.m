function [str, permutation] = adjointCurrentNames(forwardFileName)
% For each measured forward field, return the corresponding adjoint current
% name.  Return as well the proper reordering of the fields in the adjoint
% current (because in the adjoint sim, H is before E and not vice-versa).%
%
% Example: file "outFields" contains ex, dy, and bz.
% The corresponding adjoint currents are jex, jy, and mz.
% They need to be written to the adjoint source file in order as mz, jy,
% jex.
%
% [str, permutation] = adjointCurrentNames('outFields')
%
% returns
%
% str = 'mz jy jex '
% and
% permutation = [3 2 1]

of = t6.OutputFile(forwardFileName);

adjointFieldOrder = [];
adjointCurrents = {};

for ff = 1:length(of.Fields)
    fi = of.Fields{ff};
    
    number = fi.Name(2) - 'w';
    
    switch fi.Name(1)
        case 'd'
            adjointCurrents = [adjointCurrents, ['j', fi.Name(2)]];
            adjointFieldOrder = [adjointFieldOrder, 6+number];
        case 'e'
            adjointCurrents = [adjointCurrents, ['je', fi.Name(2)]];
            adjointFieldOrder = [adjointFieldOrder, 9+number];
        case 'b'
            adjointCurrents = [adjointCurrents, ['m', fi.Name(2)]];
            adjointFieldOrder = [adjointFieldOrder, 0+number];
        case 'h'
            adjointCurrents = [adjointCurrents, ['mh', fi.Name(2)]];
            adjointFieldOrder = [adjointFieldOrder, 3+number];
    end
end

[unused, permutation] = sort(adjointFieldOrder);
str = cell2mat(cellfun(@(s) [s, ' '], adjointCurrents(permutation), ...
    'UniformOutput', false));

