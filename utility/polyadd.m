function psum = polyadd(p1, p2)
% psum = polyadd(p1,p2)

if length(p1) >= length(p2)
    psum = p1;
    
    psum((1+length(p1)-length(p2)):end) = ...
        psum((1+length(p1)-length(p2)):end) + p2;
else
    psum = polyadd(p2,p1);
end