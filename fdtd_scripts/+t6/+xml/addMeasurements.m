function addMeasurements(sim, gridXML, doc, mode)

import t6.*

%directory = sim.outputDirectoryString;

%for oo = 1:length(grid.Measurements)
    %meas = grid.Measurements{oo};
    
meas = sim.Grid.Measurement;

if isempty(meas)
    return
end

if strcmpi(mode, 'forward')
    
    output = meas;
    output.type = 'Output';
    xml.addOutput(output, sim, gridXML, doc);
    
elseif strcmpi(mode, 'adjoint')
    
    fname = ['output' filesep meas.filename];
    %warning('Need to get file name in a better way');
    %[f Df] = adjoint.evalQuadFormFile(fname, meas.function);
    [f Df] = evalQuadraticFormFile(fname, meas.function);
    [adjCurrents permutedFields] = adjoint.adjointCurrentNames(fname);
    of = OutputFile(fname);
    fwdTimesteps = of.timesteps();
    adjCells = of.Regions.YeeCells(1,:);
    adjSource = squeeze(Df(:,:,:,permutedFields,end:-1:1));
    
    fwdDuration = fwdTimesteps([1 end]);
    adjDuration = sim.NumT - 1 - fwdDuration([2 1]);
    
    fieldTokens = tokenizeFields(adjCurrents, 'j m je mh');
    
    src = struct('type', 'CurrentSource', ...
        'yeeCells', adjCells, ...
        'bounds', [], ...
        'timesteps', adjDuration, ...
        'fieldFunction', [], ...
        'fieldFunctor', [], ...
        'timeData', [], ...
        'spaceTimeData', adjSource, ...
        'mode', 'adjoint', ...
        'overwrite', true);
    src.field = fieldTokens; % we have to do it separately because MATLAB SUCKS
    
    xml.addCurrentSource(src, sim, gridXML, doc);
    
else
    error('Mode must be forward or adjoint');
end
    
%    t6.xml.addOutput(output, sim, grid, gridXML, doc);
%end
