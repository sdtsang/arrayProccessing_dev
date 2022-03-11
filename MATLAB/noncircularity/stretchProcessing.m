[rx_pulses, waveform] = helperStretchSimulate;
fs = waveform.SampleRate;

helperStretchSignalSpectrogram(pulsint(rx_pulses,'coherent'),fs,...
    8,4,'Received Signal');

z_stretch = pulsint(rx_pulses,'coherent');

refrng = 6700;
rngspan = 500;
prop_speed = physconst('lightspeed');
stretchproc = getStretchProcessor(waveform,refrng,rngspan,prop_speed);

y_stretch = stretchproc(rx_pulses);
y = pulsint(y_stretch,'coherent');

helperStretchSignalSpectrogram(y,fs,16,12,'Deramped Signal');
