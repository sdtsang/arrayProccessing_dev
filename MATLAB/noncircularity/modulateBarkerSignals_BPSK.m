%  Create general BPSK modulation signals for all signals
bpskModulator = comm.BPSKModulator;
% change phase offset to pi/16
bpskModulator.PhaseOffset = pi/8;


% Barker 2
barker_2 = comm.BarkerCode('SamplesPerFrame',2,'Length',2);
seq_2 = barker_2();  % equivalent to your randi
ind_seq_2 = find(seq_2<0); seq_2(ind_seq_2) = seq_2(ind_seq_2)+1;

modData_BPSK_barker_2 = step(bpskModulator,seq_2);
channelNoise = comm.AWGNChannel('EbNo',20,'BitsPerSymbol',64);
channelOutputBpsk_barker_2 = step(channelNoise,modData_BPSK_barker_2); %channelPsk(modData);


% Barker 3
barker_3 = comm.BarkerCode('SamplesPerFrame',3,'Length',3);
seq_3 = barker_3();  % equivalent to your randi
ind_seq_3 = find(seq_3<0); seq_3(ind_seq_3) = seq_3(ind_seq_3)+1;

modData_BPSK_barker_3 = step(bpskModulator,seq_3);
channelNoise = comm.AWGNChannel('EbNo',20,'BitsPerSymbol',64);
channelOutputBpsk_barker_3 = step(channelNoise,modData_BPSK_barker_3); %channelPsk(modData);

% Barker 4
barker_4 = comm.BarkerCode('SamplesPerFrame',4,'Length',4);
seq_4 = barker_4();  % equivalent to your randi
ind_seq_4 = find(seq_4<0); seq_4(ind_seq_4) = seq_4(ind_seq_4)+1;

modData_BPSK_barker_4 = step(bpskModulator,seq_4);
channelNoise = comm.AWGNChannel('EbNo',20,'BitsPerSymbol',64);
channelOutputBpsk_barker_4 = step(channelNoise,modData_BPSK_barker_4); %channelPsk(modData);

% Barker 5
barker_5 = comm.BarkerCode('SamplesPerFrame',5,'Length',5);
seq_5 = barker_5();  % equivalent to your randi
ind_seq_5 = find(seq_5<0); seq_5(ind_seq_5) = seq_5(ind_seq_5)+1;

modData_BPSK_barker_5 = step(bpskModulator,seq_5);
channelNoise = comm.AWGNChannel('EbNo',20,'BitsPerSymbol',64);
channelOutputBpsk_barker_5 = step(channelNoise,modData_BPSK_barker_5); %channelPsk(modData);

% Barker 7
barker_7 = comm.BarkerCode('SamplesPerFrame',7,'Length',7);
seq_7 = barker_7();  % equivalent to your randi
ind_seq_7 = find(seq_7<0); seq_7(ind_seq_7) = seq_7(ind_seq_7)+1;

modData_BPSK_barker_7 = step(bpskModulator,seq_7);
channelNoise = comm.AWGNChannel('EbNo',20,'BitsPerSymbol',64);
channelOutputBpsk_barker_7 = step(channelNoise,modData_BPSK_barker_7); %channelPsk(modData);

% Barker 11
barker_11 = comm.BarkerCode('SamplesPerFrame',11,'Length',11);
seq_11 = barker_11();  % equivalent to your randi
ind_seq_11 = find(seq_11<0); seq_11(ind_seq_11) = seq_11(ind_seq_11)+1;

modData_BPSK_barker_11 = step(bpskModulator,seq_11);
channelNoise = comm.AWGNChannel('EbNo',20,'BitsPerSymbol',64);
channelOutputBpsk_barker_11 = step(channelNoise,modData_BPSK_barker_11);

% Barker 13
barker_13 = comm.BarkerCode('SamplesPerFrame',13,'Length',13);
seq_13 = barker_13();  % equivalent to your randi
ind_seq_13 = find(seq_13<0); seq_13(ind_seq_13) = seq_13(ind_seq_13)+1;

modData_BPSK_barker_13 = step(bpskModulator,seq_13);
channelNoise = comm.AWGNChannel('EbNo',20,'BitsPerSymbol',64);
channelOutputBpsk_barker_13 = step(channelNoise,modData_BPSK_barker_13);

% Barker 7+11=18
barker_18 = [seq_7 ; seq_11];
%modData_QPSK_barker_18 = step(qpskModulator,barker_18);
modData_BPSK_barker_18 = step(bpskModulator,barker_18);
channelNoise = comm.AWGNChannel('EbNo',20,'BitsPerSymbol',64);
channelOutputBpsk_barker_18 = step(channelNoise,modData_BPSK_barker_18);

% Barker 11+13=24
barker_24 = [seq_11 ; seq_13];
modData_BPSK_barker_24 = step(bpskModulator,barker_24);
channelNoise = comm.AWGNChannel('EbNo',20,'BitsPerSymbol',64);
channelOutputBpsk_barker_24 = step(channelNoise,modData_BPSK_barker_24);

% Barker 18+13=31
barker_31 = [barker_18 ; seq_13];
%modData_QPSK_barker_18 = step(qpskModulator,barker_18);
modData_BPSK_barker_31 = step(bpskModulator,barker_31);
channelNoise = comm.AWGNChannel('EbNo',20,'BitsPerSymbol',64);
channelOutputBpsk_barker_31 = step(channelNoise,modData_BPSK_barker_31);
