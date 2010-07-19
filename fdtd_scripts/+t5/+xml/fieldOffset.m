function xyzt = fieldOffset(fieldString)

xzyt = [0 0 0 0];

switch fieldString
    case 'ex'
        xyzt = [0.5 0 0 0];
    case 'ey'
        xyzt = [0 0.5 0 0];
    case 'ez'
        xyzt = [0 0 0.5 0];
    case 'hx'
        xyzt = [0 0.5 0.5 0.5];
    case 'hy'
        xyzt = [0.5 0 0.5 0.5];
    case 'hz'
        xyzt = [0.5 0.5 0 0.5];
    case 'dx'
        xyzt = [0.5 0 0 0];
    case 'dy'
        xyzt = [0 0.5 0 0];
    case 'dz'
        xyzt = [0 0 0.5 0];
    case 'bx'
        xyzt = [0 0.5 0.5 0.5];
    case 'by'
        xyzt = [0.5 0 0.5 0.5];
    case 'bz'
        xyzt = [0.5 0.5 0 0.5];
    case {'jx', 'jex'}
        xyzt = [0.5 0 0 0.5];
    case {'jy', 'jey'}
        xyzt = [0 0.5 0 0.5];
    case {'jz', 'jez'}
        xyzt = [0 0 0.5 0.5];
    case {'mx', 'mhx'}
        xyzt = [0 0.5 0.5 0];
    case {'my', 'mhy'}
        xyzt = [0.5 0 0.5 0];
    case {'mz', 'mhz'}
        xyzt = [0.5 0.5 0 0];
end
