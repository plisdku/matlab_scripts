function n = getIndexOfStringInCellArray(cellArray, stringToFind)

n = -1;
for nn = 1:length(cellArray)
    if strcmp(cellArray{nn}.Name, stringToFind)
        n = nn;
    end
end

