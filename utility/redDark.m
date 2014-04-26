function cm = redDark()

v = linspace(0,1,1000);
cm = ([sqrt(v); v.^2; v.^2]').^0.8;