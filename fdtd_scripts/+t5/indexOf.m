function index = indexOf(nameOfStruct, inCellArray)
%indexOf Search inCellArray for a struct with the given .name or .Name field
%   This function is internal to Trogdor.
index = -1;
for nn = 1:length(inCellArray)
    if isfield(inCellArray{nn}, 'name')
        if strcmp(inCellArray{nn}.name, nameOfStruct)
            index = nn;
        end
    elseif strcmp(inCellArray{nn}.Name, nameOfStruct)
            index = nn;
        end
    end
end
