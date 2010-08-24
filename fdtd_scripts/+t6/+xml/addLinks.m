function addLinks(grid, gridXML, doc, originTrogdor)

originTwice = [originTrogdor originTrogdor];

for ll = 1:length(grid.Links)
    link = grid.Links{ll};
    
    srcOriginTwice = repmat(link.sourceTrogdorOrigin, [1 2]);
    
    elemXML = doc.createElement('Link');
    elemXML.setAttribute('type', 'TF');  % should I expose SF links too?
    elemXML.setAttribute('sourceGrid', link.sourceGrid);
    elemXML.setAttribute('fromYeeCells', ...
        sprintf('%i ', link.sourceYeeCells + srcOriginTwice));
    elemXML.setAttribute('toYeeCells', ...
        sprintf('%i ', link.destYeeCells + originTwice));
    
    for oo = 1:length(link.omitSides)
        omitXML = doc.createElement('OmitSide');
        omitXML.appendChild(doc.createTextNode(...
            sprintf('%i ', link.omitSides{oo})));
        elemXML.appendChild(omitXML);
    end
    
    gridXML.appendChild(elemXML);
end
