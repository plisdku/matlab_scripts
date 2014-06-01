function c = allWords(numSymbols, wordLength)
% c = allWords(numSymbols, wordLength)
% Return an exhaustive list of all "words" of a given length using symbols
% 0:numSymbols-1.
%
% EXAMPLE: The digits of the binary numbers from 0 to 7.
%
%   >> allWords(2,3)
% 
%     ans =
% 
%          0     0     0
%          0     0     1
%          0     1     0
%          0     1     1
%          1     0     0
%          1     0     1
%          1     1     0
%          1     1     1
%

% slow, limited, built-in way...
%c = dec2base(0:numSymbols^wordLength-1, numSymbols) - '0';

numWords = numSymbols^wordLength;

c = zeros(numSymbols^wordLength, wordLength);

sequence = (0:numWords-1)';

divisor = 1;
for dim = wordLength:-1:1
    c(:,dim) = floor(rem(sequence/divisor, numSymbols));
    divisor = divisor * numSymbols;
end