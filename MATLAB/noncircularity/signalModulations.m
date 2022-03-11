% modData_BPSK and y_64
%qam_64 = y_64;
% sigType = input('User, select signal type?\n 1.) BPSK, 2.) QPSK 3.) PSK, M=8 4.) PSK,M=16 5.) QAM, M=16 6.) QAM, M=64 7.) QAM, M=256, 8.) pammod 9.) aamod 10.) ask\n');
% 
stringCase = 'BPSK';
% +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
% bpsk modulation scheme
% Create a BPSK modulator System Object
bpskModulator = comm.BPSKModulator;
% change phase offset to pi/16
bpskModulator.PhaseOffset = 0;
% create binary data symbols
data = randi([0 1],100,1);
% modulate the data
% modData = bpskModulator(data);
modData_BPSK = step(bpskModulator,data);
channelBPsk = comm.AWGNChannel('EbNo',20,'BitsPerSymbol',64);
channelOutputBpsk = step(channelBPsk,modData_BPSK); %channelPsk(modData);

%%%%
%sig = modData_BPSK;

stringCase = 'QPSK';
% +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
% QPSK Modulation Scheme
mod = comm.QPSKModulator;
refC = constellation(mod);
constellation(mod);

% Create a QPSK demodulator having a 0 phase offset
demod = comm.QPSKDemodulator('PhaseOffset',0);
constellation(demod);

% create a QPSK modulator object and a phase noise object
qpskModulator = comm.QPSKModulator;
phNoise = comm.PhaseNoise('Level',-55,'FrequencyOffset',20,'SampleRate',1000);
channelNoise = comm.AWGNChannel('EbNo',10,'BitsPerSymbol',16);
d = randi([0 3],1000,1);
modData_Qpsk = step(qpskModulator,d);
%         channelOutputQpsk = step(phNoise,modData_Qpsk);
channelOutputQpsk = step(channelNoise,modData_Qpsk);
%sig= modData_Qpsk;

stringCase = 'PSK, M=8';
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

%%%%%
%sig = modData_PSK;

stringCase = 'PSK, M=16';
% 16 PSK example
M=16;
data = randi([0 M-1],1000,1);
modData_PSK_16 = pskmod(data,M,pi/M);
channelOutputPSK_16 = awgn(modData_PSK_16,25);
%scatterplot(rxSig)
%sig= modData_PSK_16;

stringCase = 'QAM, M=16';
% Modulation order of 16 =4^2
M_16 = 16;
x_16 = (0:M_16-1)';
y_16 = qammod(x_16,M_16,'bin','UnitAveragePower',true,'PlotConstellation',true);
%scatterplot(y_16); grid on;
%sig = y_16;

stringCase = 'QAM, M=64';
% M = 64 4^3
M_64 = 64;
x_64 = (0:M_64-1)';
y_64 = qammod(x_64,M_64,'bin','UnitAveragePower',true,'PlotConstellation',true);
%scatterplot(y_64); grid on;
%sig = y_64;

stringCase = 'QAM, M=256';
% M = 256 = 4^4
M_256 = 256;
x_256 = (0:M_256-1)';
y_256 = qammod(x_256,M_256,'bin','UnitAveragePower',true,'PlotConstellation',true);
%scatterplot(y_256); grid on;
%sig = y_256;

stringCase = 'PAM, M=8';
M = 8;
data = randi([0 M-1],100,1);
modData_PAM = pammod(data,M,pi/8);
%sig = modData_PAM;

stringCase = 'AAMOD';
fs = 100;
t = (0:1/fs:100)';
fc = 10;
x = sin(2*pi*t);
ydouble = ammod(x,fc,fs);
ysingle = ssbmod(x,fc,fs);
sa = dsp.SpectrumAnalyzer('SampleRate',fs, ...
    'PlotAsTwoSidedSpectrum',false, ...
    'YLimits',[-60 40]);
step(sa,ysingle);
%         
%         
% %         b = [0 1 0 1 1 1 0];
% %         n = length(b);
% %         t = 0:.01:n;
% %         x = 1:1:(n+1)*100;
% %         for i = 1:n
% %             for j = i:.1:i+1
% %                 bw(x(i*100:(i+1)*100)) = b(i); %#ok<SAGROW>
% %             end
% %         end
% %         bw = bw(100:end);
% %         sint = sin(2*pi*t);
% %         st = bw.*sint;
% %         sig = st;
% end