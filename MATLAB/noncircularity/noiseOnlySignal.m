y1 = wgn(1000,1,-30,'complex');
var(y1) % \approx 0.001

y2 = wgn(1000,1,10,'complex');
var(y2) % \approx 10.2073

y3 = wgn(1000,1,0,'real');
var(y3)  % \approx 1

y4 = wgn(1000,1,20,'real');
var(y4)

y5 = wgn(1000,1,30,'real');
var(y5)

y6 = wgn(1000,1,25,'real');
var(y6)

y7 = wgn(1000,1,60,'real');
var(y7)