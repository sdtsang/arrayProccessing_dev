% Modulate Data Using QAM:

% Signals with M-QAM and M = 4^k are circular
% Modulation order of 16 =4^2
M_16 = 16;
x_16 = randi([0 M_16-1], 1000,1)';
y_16 = qammod(x_16,M_16,'bin','UnitAveragePower',true,'PlotConstellation',true);
scatterplot(y_16); grid on;

% M = 64 4^3
M_64 = 64;
x_64 = randi([0 M_64-1], 1000,1)';
y_64 = qammod(x_64,M_64,'bin','UnitAveragePower',true,'PlotConstellation',true);
scatterplot(y_64); grid on;

% M = 256 = 4^4
M_256 = 256;
x_256 = randi([0 M_256-1], 1000,1)';
y_256 = qammod(x_256,M_256,'bin','UnitAveragePower',true,'PlotConstellation',true);
scatterplot(y_256); grid on;

% Add Phase Offset and Noise
bps = 6;
M = 2^bps; % 64-QAM
data = randi([0 1],1008,1);
H = comm.RectangularQAMModulator('ModulationOrder',64,'BitInput',true);
H.PhaseOffset = 0;
H.ModulationOrder = 64;
dataMod = qammod(data,M,'InputType','bit');
 dataModSO = H(data); %         sig= dataModSO
%  scatterplot(dataModSO);
%  scatterplot(dataMod);
channelNoise = comm.AWGNChannel('EbNo',3,'BitsPerSymbol',256);
channelOutputQAM = step(channelNoise,dataModSO); %channelPsk(modData);
sig = channelOutputQAM
 