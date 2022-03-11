% waveform1 = {'Rectangular','PRF',10e3 'PulseWidth',50e-6};
% waveform2 = {'LinearFM','PRF',10e3,'PulseWidth',50e-6,'SweepBandwidth',1e5, ...
%     'SweepDirection','Up','SweepInterval', 'Positive'};
waveform1 = {'PhaseCoded','PRF',1e4,'Code','Zadoff-Chu', ...
    'SequenceIndex',3,'ChipWidth',5e-6,'NumChips',8};
waveform2 = {'PhaseCoded','PRF',1e4,'Code','Zadoff-Chu', ...
    'SequenceIndex',3,'ChipWidth',5e-6,'NumChips',8};
pulsesib = phased.PulseWaveformLibrary('SampleRate',1e6,...
    'WaveformSpecification',{waveform1,waveform2});

coeff1 = getMatchedFilter(pulsesib,1,1);
subplot(2,1,1)
stem(real(coeff1))
title('Matched filter coefficients, real part')
coeff2 = getMatchedFilter(pulsesib,2,1);
subplot(2,1,2)
stem(real(coeff2))
title('Matched filter coefficients, real part')

filter = phased.MatchedFilter('Coefficients',wav);
taylorfilter = phased.MatchedFilter('Coefficients',wav,...
    'SpectrumWindow','Taylor');

test_nc = [real(pulsesib(1)); imag(pulsesib(1))];
wavelib = pulsesib(1);

filter_nc = phased.MatchedFilter('Coefficients',wavelib);

y_nc=filter_nc(test_nc);

coeff_circ = filter_nc(pulsesib(1));

figure(1)
stem(real(coeff_circ),'k')
hold on
stem(imag(coeff_circ),'g')
figure(2)
stem(real(y_nc),'m')
hold on
stem(imag(y_nc),'r')

figure(3)
plot(pulsesib,1,'PlotType','complex')

figure
antenna = phased.IsotropicAntennaElement('FrequencyRange',[5e9 15e9]);
transmitter = phased.Transmitter('Gain',20,'InUseOutputPort',true);
fc = 10e9;
target = phased.RadarTarget('Model','Nonfluctuating',...
   'MeanRCS',1,'OperatingFrequency',fc);
txloc = [0;0;0];
tgtloc = [5000;5000;10];
transmitterplatform = phased.Platform('InitialPosition',txloc);
targetplatform = phased.Platform('InitialPosition',tgtloc);
[tgtrng,tgtang] = rangeangle(targetplatform.InitialPosition,...
   transmitterplatform.InitialPosition);
% 
% waveform = phased.RectangularWaveform('PulseWidth',25e-6,...
%    'OutputFormat','Pulses','PRF',10e3,'NumPulses',1);
% waveform = {'PhaseCoded','PRF',1e4,'Code','Zadoff-Chu', ...
%     'SequenceIndex',3,'ChipWidth',5e-6,'NumChips',8};
waveform = phased.PhaseCodedWaveform('Code','Zadoff-Chu',...
    'ChipWidth',5e-6,'NumChips',16,...
    'OutputFormat','Pulses','NumPulses',2);

c = physconst('LightSpeed');
maxrange = c/(2*waveform.PRF);
SNR = npwgnthresh(1e-6,1,'noncoherent');
Pt = radareqpow(c/fc,maxrange,SNR,...
   waveform.ChipWidth,'RCS',target.MeanRCS,'Gain',transmitter.Gain);

transmitter.PeakPower = Pt;

radiator = phased.Radiator('PropagationSpeed',c,...
   'OperatingFrequency',fc,'Sensor',antenna);
channel = phased.FreeSpace('PropagationSpeed',c,...
   'OperatingFrequency',fc,'TwoWayPropagation',false);
collector = phased.Collector('PropagationSpeed',c,...
   'OperatingFrequency',fc,'Sensor',antenna);
receiver = phased.ReceiverPreamp('NoiseFigure',0,...
   'EnableInputPort',true,'SeedSource','Property','Seed',2e3);
filter = phased.MatchedFilter(...
   'Coefficients',getMatchedFilter(waveform),...
   'GainOutputPort',true);

wf = waveform();
[wf,txstatus] = transmitter(wf);
wf = radiator(wf,tgtang);
wf = channel(wf,txloc,tgtloc,[0;0;0],[0;0;0]);
wf = target(wf);
wf = channel(wf,tgtloc,txloc,[0;0;0],[0;0;0]);
wf = collector(wf,tgtang);
rx_puls = receiver(wf,~txstatus);
[mf_puls,mfgain] = filter(rx_puls);

Gd = length(filter.Coefficients)-1;

mf_puls=[mf_puls(Gd+1:end); mf_puls(1:Gd)];
subplot(2,1,1)
t = unigrid(0,1e-6,2e-4,'[)');
rangegates = c.*t;
rangegates = rangegates/2;
plot(rangegates,abs(rx_puls))
title('Received Pulse')
ylabel('Amplitude')
hold on
plot([tgtrng, tgtrng], [0 max(abs(rx_puls))],'r')
subplot(2,1,2)
plot(rangegates,abs(mf_puls))
title('With Matched Filtering')
xlabel('Meters')
ylabel('Amplitude')
hold on
plot([tgtrng, tgtrng], [0 max(abs(mf_puls))],'r')
hold off