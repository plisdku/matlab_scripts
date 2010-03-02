function addMaterials(documentNode, simHandle)
doc = documentNode;
root = doc.getDocumentElement;
sim = simHandle;

for mm = 1:length(sim.Materials)
    matXML = doc.createElement('Material');
    matXML.setAttribute('name', sim.Materials{mm}.name);
    matXML.setAttribute('model', sim.Materials{mm}.model);
    paramXML = doc.createElement('Params');
    fields = fieldnames(sim.Materials{mm});
    for ff = 1:length(fields)
    if ~strcmp(fields(ff), 'name') && ~strcmp(fields(ff), 'model') &&...
        ~strcmp(fields(ff), 'PMLParams')
        paramXML.setAttribute(fields{ff}, ...
            num2str(getfield(sim.Materials{mm}, fields{ff}), '%2.8g '));
    end
    end
    matXML.appendChild(paramXML);
    root.appendChild(matXML);
end
