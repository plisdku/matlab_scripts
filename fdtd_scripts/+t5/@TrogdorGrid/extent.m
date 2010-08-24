function yeeCells = extent(obj)

yeeCells = [Inf Inf Inf -Inf -Inf -Inf];

for aa = 1:length(obj.Assembly)
    if isfield(obj.Assembly{aa}, 'yeeCells')
        objCells = obj.Assembly{aa}.yeeCells;
        yeeCells(1:3) = min(yeeCells(1:3), obj.Assembly{aa}.yeeCells(1:3));
        yeeCells(4:6) = max(yeeCells(4:6), obj.Assembly{aa}.yeeCells(4:6));
    end
end

if ~t5.validateRect(yeeCells)
    error('Cannot calculate proper extent of grid.');
end
