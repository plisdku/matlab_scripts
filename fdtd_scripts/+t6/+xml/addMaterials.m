function addMaterials(doc, sim)

root = doc.getDocumentElement;

for mm = 1:length(sim.Materials)
    matXML = doc.createElement('Material');
    matXML.setAttribute('name', sim.Materials{mm}.name);
    matXML.setAttribute('zNumerator', num2str(sim.Materials{mm}.numerator, 15));
    matXML.setAttribute('zDenominator', num2str(sim.Materials{mm}.denominator, 15));
    
    root.appendChild(matXML);
end
