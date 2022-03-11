% ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
% FSK modulation and demodulation in AWGN
M=8;
freqSep=100;

% FSK modulator and demodulator system objects
fskMod = comm.FSKModulator(M,freqSep);
fskDemod=comm.FSKDemodulator(M,freqSep);

% additive white Gaussian noise channel
ch = comm.AWGNChannel('NoiseMethod', ...
    'Signal to noise ratio (SNR)','SNR',-2);
err = comm.ErrorRate; % error rate calculator object

% transmit one hundred 50-symbol frames using 8-FSK in AWGN channel
for counter = 1:100
    data = randi([0 M-1],50,1);
    modSignal = step(fskMod,data);
    noisySignal = step(ch,modSignal);
    receivedData = step(fskDemod,noisySignal);
    errorStats = step(err,data,receivedData);
end

% error statistics
es = 'Error rate = %4.2e\nNumber of errors = %d\nNumber of symbols = %d\n';
fprintf(es,errorStats);
% +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
% PSK Modulation; AWGN to 8-PSK signal example

% Create psk modulator system object; default order is 8
pskModulator = comm.PSKModulator;   
% pskModulator.ModulationOrder=2;  % change order

% Modulate the signal
tempArr = randi([0 7],2000,1);
modData_PSK = step(pskModulator,tempArr);

% AWGN to the signal by passing it through a noise channel
channelPsk = comm.AWGNChannel('EbNo',20,'BitsPerSymbol',3);

% Transmit the signal through the AWGN channel
channelOutputPsk = step(channelPsk,modData_PSK); %channelPsk(modData);

% plot noiseless and noisy data using scatterplot
figure;
scatterplot(modData_PSK); grid on;
set(gca,'FontWeight','bold','FontSize',12);

% figure;
scatterplot(channelOutputPsk); grid on;
set(gca,'FontWeight','bold','FontSize',12);

% figure;
constellation(pskModulator);
set(gca,'FontWeight','bold','FontSize',12);

% change noise level and compare
channelPsk.EbNo = 10;
channelOutputPsk10 = step(channelPsk,modData_PSK); %channelPsk(modData);

% figure;
scatterplot(channelOutputPsk10); legend('Noise: 10 dB');
set(gca,'FontWeight','bold','FontSize',12);