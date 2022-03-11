% generate maximal length pn sequences
data = randi([0 1],1000,1);
%modData = qammod(data,64);
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



% define orders of the generator polynomial
orderPoly = [1 2 3 4 5];

calculateMaximalSequenceLengths = 2.^orderPoly-1;

% determine samples per frame for each
samplesPerFrame = calculateMaximalSequenceLengths*2;

% define set of initial conditions for each
initCond_1 = [0];
initCond_2 = [0 1];
initCond_3 = [0 0 1];
initCond_4 = [0 0 0 1];
initCond_5 = [0 0 0 0 1];
initCond_6 = [0 0 0 0 0 1];

% poly3Type = {[3 2 0]
%              [2 3 0]
%              [3 1 0]
%              [1 3 0]
%              [0 1 3]
%              [0 3 1]
%              [0 3 2]
%              [0 2 3]
%              [1 0 3]
%              [3 0 1]
%              [2 0 3]
%              [3 0 2]
%             };

poly3Type = {
             [3 2 0]
             [3 1 0]
             [1 3 0]
             [0 1 3]
             [0 3 1]
             [0 3 2]
             [0 2 3]
             [1 0 3]
             [3 0 1]
             [2 0 3]
             [3 0 2]
            };

pnSequence = comm.PNSequence('Polynomial',[1 0], ...
    'SamplesPerFrame',samplesPerFrame(1),'InitialConditions',initCond_1);
x1 = pnSequence();
        
pnSequence = comm.PNSequence('Polynomial',[ 2 1 0], ...
    'SamplesPerFrame',samplesPerFrame(2),'InitialConditions',initCond_2);
x2 = pnSequence();
        
pnSequence = comm.PNSequence('Polynomial',[3 2 1 0], ...
    'SamplesPerFrame',samplesPerFrame(3),'InitialConditions',initCond_3);
x3 = pnSequence();
 
pnSequence = comm.PNSequence('Polynomial',[4 3 2 1 0], ...
    'SamplesPerFrame',samplesPerFrame(4),'InitialConditions',initCond_4);
x4 = pnSequence();

pnSequence = comm.PNSequence('Polynomial',[5 4 3 2 1 0], ...
    'SamplesPerFrame',samplesPerFrame(5),'InitialConditions',initCond_5);
x5 = pnSequence();


for i = 1:1
pnSequence = comm.PNSequence('Polynomial',[5 4 3 2 1 0], ...
    'SamplesPerFrame',samplesPerFrame(5),'InitialConditions',initCond_5);
x5(:,i) = pnSequence();
dlmwrite('x5.txt',x5(:,i)','delimiter','','-append');
end


pnSequence = comm.PNSequence('Polynomial',[3 2 0], ...
    'SamplesPerFrame',14,'InitialConditions',[0 0 1]);
x1 = pnSequence();

% generate Gold sequence
goldseq = comm.GoldSequence('FirstPolynomial','x^5+x^2+1',...
    'SecondPolynomial','x^5+x^4+x^3+x^2+1',...
    'FirstInitialConditions',[0 0 0 0 1],...
    'SecondInitialConditions',[0 0 0 0 1],...
    'Index',4,'SamplesPerFrame',10);
x = goldseq();

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

modData=channelOutputQpsk;

gold_mod = modData.*x';
x2_mod = modData.*x2';
x3_mod = modData.*x3';
x4_mod = modData.*x4';
x5_mod = modData.*x5';


% pnSequence2 = comm.PNSequence('Polynomial','x^4+x+1', ...
%     'InitialConditions',[0 0 0 1],'SamplesPerFrame',30);
% x2 = pnSequence2();

% PN Sequences are "deterministic" in that they always return the same results if given the same input values.
