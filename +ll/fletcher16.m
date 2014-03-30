function checksum = fletcher16(A, keepFlag)
% fletcher16   Calculate Fletcher-16 checksum of data
% checksum = fletcher16(A)
% checksum = fletcher16(A, keepFlag)
%
% If keepFlag is true, the checksum will continue to accumulate from
% previous values.

persistent sum1 sum2

if nargin < 2
    keepFlag = false;
end

if ~keepFlag
    sum1 = uint16(0);
    sum2 = uint16(0);
end

A = typecast(A(:), 'uint8');

% slow... barf.

for ii = 1:numel(A)
    %fprintf('Byte(%i) = %x\n', ii, A(ii));
    sum1 = mod(sum1 + uint16(A(ii)), 255);
    %fprintf('C0(%i) = %x\n', ii, sum1);
    sum2 = mod(sum2 + sum1, 255);
    %fprintf('C1(%i) = %x\n', ii, sum2);
end

assert(isa(sum1, 'uint16'));
assert(isa(sum2, 'uint16'));

checksum = bitor(bitshift(sum2, 8), sum1);