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
    case 'jx'
        xyzt = [0.5 0 0 0.5];
    case 'jy'
        xyzt = [0 0.5 0 0.5];
    case 'jz'
        xyzt = [0 0 0.5 0.5];
    case 'kx'
        xyzt = [0 0.5 0.5 0];
    case 'ky'
        xyzt = [0.5 0 0.5 0];
    case 'kz'
        xyzt = [0.5 0.5 0 0.5];
end
