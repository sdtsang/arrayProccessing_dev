% +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
% bpsk modulation scheme
% Create a BPSK modulator System Object
bpskModulator = comm.BPSKModulator;
% change phase offset to pi/16
bpskModulator.PhaseOffset = pi/8;
% create binary data symbols
data = randi([0 1],100,1);
% modulate the data
% modData = bpskModulator(data);
modData_BPSK = step(bpskModulator,data);
channelBPsk = comm.AWGNChannel('EbNo',20,'BitsPerSymbol',64);
channelOutputBpsk = step(channelBPsk,modData_BPSK); %channelPsk(modData);

% figure;
scatterplot(modData_BPSK);
set(gca,'FontWeight','bold','FontSize',12);

% figure;
scatterplot(channelOutputBpsk);
set(gca,'FontWeight','bold','FontSize',12);

% figure;
constellation(bpskModulator);

