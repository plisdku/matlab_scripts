function cm = blueDark()

v = linspace(0,1,1000)';
cm = ([v.^2, v.^2, sqrt(v)]).^0.4;