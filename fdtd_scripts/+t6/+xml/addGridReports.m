function addGridReports(grid, gridXML, doc, originTrogdor)

originTwice = [originTrogdor originTrogdor];

for gg = 1:length(grid.GridReports)
    rep = grid.GridReports{gg};
    
    elemXML = doc.createElement('GridReport');
    elemXML.setAttribute('file', rep.filename);
    % no further options are permitted at present
    
    if length(rep.yeeCells) ~= 0
        regionXML = doc.createElement('Region');
        regionXML.setAttribute('yeeCells', ...
            sprintf('%i ', rep.yeeCells + originTwice));
        elemXML.appendChild(regionXML);
    end
    gridXML.appendChild(elemXML);
end
