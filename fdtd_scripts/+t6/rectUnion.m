function r = rectUnion(r1, r2)
%rectUnion Calculate the rect encompassing both rects

r = [min(r1(1:3), r2(1:3)), max(r1(4:6), r2(4:6))];