fs = 100;
t = (0:1/fs:100)';
fc = 10;
x = sin(2*pi*t);
ydouble = ammod(x,fc,fs);
ysingle = ssbmod(x,fc,fs);
sa = dsp.SpectrumAnalyzer('SampleRate',fs, ...
    'PlotAsTwoSidedSpectrum',false, ...
    'YLimits',[-60 40]);
step(sa,ydouble);
set(gca,'FontWeight','bold','FontSize',12);
title('Spectrum of Amplitude-Modulated Double-Sideband Signal, carrier frequency fc=10 Hz','FontWeight','bold','FontSize',12,'Color','k')

step(sa,ysingle)
title('Spectrum of Amplitude-Modulated Single-Sideband Signal, carrier frequency fc=10 Hz','FontWeight','bold','FontSize',12,'Color','k')
set(gca,'FontWeight','bold','FontSize',12);