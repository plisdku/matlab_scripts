function rr = uniqueRows(m)

if isempty(m)
    rr = [];
    return;
end

mSorted = sortrows(m);

lastUniqueRow = 1;
rr = mSorted(lastUniqueRow,:);

for row = 2:size(mSorted,1)
    if any(mSorted(row,:) ~= mSorted(lastUniqueRow,:))
        rr(end+1,:) = mSorted(row,:);
        lastUniqueRow = row;
    end
end
