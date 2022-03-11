limVals = 1.5;

% +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
% bpsk modulation scheme
% Create a BPSK modulator System Object
bpskModulator = comm.BPSKModulator;
% change phase offset to pi/16
bpskModulator.PhaseOffset = pi/16;
%  bpskModulator.PhaseOffset = 0;

% create binary data symbols
data = randi([0 1],1000,1);
% modulate the data
% modData = bpskModulator(data);
modData_BPSK = step(bpskModulator,data);

channelBPsk_0 = comm.AWGNChannel('EbNo',0,'BitsPerSymbol',16);
channelOutputBpsk_0 = step(channelBPsk_0,modData_BPSK); %channelPsk(modData);
channelBPsk_3 = comm.AWGNChannel('EbNo',3,'BitsPerSymbol',16);
channelOutputBpsk_3 = step(channelBPsk_3,modData_BPSK); %channelPsk(modData);
channelBPsk_10 = comm.AWGNChannel('EbNo',10,'BitsPerSymbol',16);
channelOutputBpsk_10 = step(channelBPsk_10,modData_BPSK); %channelPsk(modData);
channelBPsk_15 = comm.AWGNChannel('EbNo',15,'BitsPerSymbol',16);
channelOutputBpsk_15 = step(channelBPsk_15,modData_BPSK); %channelPsk(modData);

% BPSK
figure(1)
scatter(real(channelOutputBpsk_3),imag(channelOutputBpsk_3),'k.','LineWidth',2 ); hold on;
scatter(real(channelOutputBpsk_10),imag(channelOutputBpsk_10),'g.','LineWidth',2 ); hold on;
scatter(real(channelOutputBpsk_15),imag(channelOutputBpsk_15),'m.','LineWidth',2 ); hold on;
scatter(real(modData_BPSK),imag(modData_BPSK),'b.','LineWidth',4); hold on;

set(gca,'FontWeight','bold','FontSize',12);
ylim([-limVals limVals]); xlim([-limVals limVals]);
axis equal; grid on; box on; hold on;
xlabel('In-Phase Amplitude'); ylabel('Quadrature Amplitude');
leg = legend('3','10','15','0'); title(leg,'E_b/N_0, dB')
title('BPSK with Various Levels of White Gaussian Noise')

% BPSK
figure(2)
scatter(real(modData_BPSK),imag(modData_BPSK),'b.','LineWidth',8); hold on;
set(gca,'FontWeight','bold','FontSize',12);
ylim([-limVals limVals]); xlim([-limVals limVals]);
axis equal; grid on; box on; hold on;
xlabel('In-Phase Amplitude'); ylabel('Quadrature Amplitude');
title('BPSK Modulation Scheme, with 0 phase offset')
%set(gca,'XTick',-limVals:0.5:limVals)


% PSK, M = 8
% Create psk modulator system object; default order is 8
pskModulator = comm.PSKModulator;
% pskModulator.ModulationOrder=2;  % change order

% Modulate the signal
tempArr = randi([0 7],2000,1);
modData_PSK = step(pskModulator,tempArr);

% AWGN to the signal by passing it through a noise channel
channelPsk = comm.AWGNChannel('EbNo',8,'BitsPerSymbol',16);

% Transmit the signal through the AWGN channel
channelOutputPsk = step(channelPsk,modData_PSK); %channelPsk(modData);
% M = 8
figure(3)
scatter(real(channelOutputPsk),imag(channelOutputPsk),'k.','LineWidth',1 ); hold on;
scatter(real(modData_PSK),imag(modData_PSK),'g.','LineWidth',8); hold on;
set(gca,'FontWeight','bold','FontSize',12);
ylim([-limVals limVals]); xlim([-limVals limVals]);
axis equal; grid on; box on; hold on;
xlabel('In-Phase Amplitude'); ylabel('Quadrature Amplitude');
title('PSK Modulation Scheme, M=16')
leg = legend('E_b/N_0=20 dB' , 'reference'); %title(leg,'E_b/N_0, dB')

% QPSK
figure(4)
scatter(real(channelOutputQpsk),imag(channelOutputQpsk),'k.','LineWidth',2 ); hold on;
scatter(real(modData_Qpsk),imag(modData_Qpsk),'g.','LineWidth',8); hold on;
set(gca,'FontWeight','bold','FontSize',12);
ylim([-limVals limVals]); xlim([-limVals limVals]);
axis equal; grid on; box on; hold on;
xlabel('In-Phase Amplitude'); ylabel('Quadrature Amplitude');
title('QPSK Modulation Scheme')
leg = legend('E_b/N_0=10 dB','reference'); %title(leg,'E_b/N_0, dB')

% PSK, M = 16
figure(5)
scatter(real(channelOutputPSK_16),imag(channelOutputPSK_16),'k.','LineWidth',2 ); hold on;
scatter(real(modData_PSK_16),imag(modData_PSK_16),'g.','LineWidth',8); hold on;
set(gca,'FontWeight','bold','FontSize',12);
ylim([-limVals limVals]); xlim([-limVals limVals]);
axis equal; grid on; box on; hold on;
xlabel('In-Phase Amplitude'); ylabel('Quadrature Amplitude');
title('PSK Modulation Scheme, M = 16')
leg = legend('25','0'); title(leg,'E_b/N_0, dB')


% QAM_M_16
figure(6)
%scatter(real(x_16),imag(x_16),'k.','LineWidth',2 ); hold on;
scatter(real(y_16),imag(y_16),'k.','LineWidth',16); hold on;
set(gca,'FontWeight','bold','FontSize',12);
ylim([-limVals limVals]); xlim([-limVals limVals]);
axis equal; grid on; box on; hold on;
xlabel('In-Phase Amplitude'); ylabel('Quadrature Amplitude');
title('QAM Modulation Scheme, M = 16')
%leg = legend('25','0'); title(leg,'E_b/N_0, dB')


% QAM_M_64
figure(7)
%scatter(real(x_16),imag(x_16),'k.','LineWidth',2 ); hold on;
scatter(real(y_64),imag(y_64),'k.','LineWidth',16); hold on;
set(gca,'FontWeight','bold','FontSize',12);
ylim([-limVals limVals]); xlim([-limVals limVals]);
axis equal; grid on; box on; hold on;
xlabel('In-Phase Amplitude'); ylabel('Quadrature Amplitude');
title('QAM Modulation Scheme, M = 64')
%leg = legend('25','0'); title(leg,'E_b/N_0, dB')


% QAM_M_64
figure(8)
%scatter(real(x_16),imag(x_16),'k.','LineWidth',2 ); hold on;
scatter(real(y_256),imag(y_256),'ko','LineWidth',2); hold on;
scatter(real(y_64),imag(y_64),'b+','LineWidth',2); hold on;
scatter(real(y_16),imag(y_16),'rx','LineWidth',2); hold on;

%scatter(real(dataModSO),imag(dataModSO),'k*','LineWidth',2); hold on;

scatter(real(channelOutputQAM),imag(channelOutputQAM),'ms','LineWidth',2); hold on;


set(gca,'FontWeight','bold','FontSize',12);
ylim([-limVals limVals]); xlim([-limVals limVals]);
axis equal; grid on; box on; hold on;
xlabel('In-Phase Amplitude'); ylabel('Quadrature Amplitude');
title('QAM Modulation Scheme at different Modulation Orders, M')
leg = legend('256','64','16'); title(leg,'M')


% PAMMOD, M = 8
figure(9)
%scatter(real(channelOutputPSK_16),imag(channelOutputPSK_16),'k.','LineWidth',2 ); hold on;
scatter(real(modData_PAM),imag(modData_PAM),'k.','LineWidth',8); hold on;
%scatterplot(modData_PAM);
set(gca,'FontWeight','bold','FontSize',12);
%ylim([-limVals limVals]); xlim([-limVals limVals]);
axis equal; grid on; box on; hold on;
xlabel('In-Phase Amplitude'); ylabel('Quadrature Amplitude');
title({'Pulse Amplitude Modulation Scheme', 'M = 8, \pi/4 Offset'})
%leg = legend('25','0'); title(leg,'E_b/N_0, dB')