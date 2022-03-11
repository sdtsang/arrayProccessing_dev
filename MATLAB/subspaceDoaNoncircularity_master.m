close all; clearvars; clc;
% Master script for noncircular sources analysis

%% ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
% Initialization
theMathsPath = '/home/sdtsang/Documents/Research/2021_2022/maths/';
% Doesn't work -- FIX
% s       = pathsep;
% pathStr = [s, path, s];
% onPath  = contains(pathStr, [s, theMathsPath, s], 'IgnoreCase', ispc);
% if ~onPath, addpath(genpath(theMathsPath)); end

%% ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
% 0.) Generate noncircular signals
% 1) BPSK 2) QPSK
NC_sigType = 1; % 2
M = 3; % number of sources of incoming signals

switch NC_sigType
    case 1
        % a) +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        % bpsk modulation scheme
        NcStringCase = 'BPSK';
        % Create a BPSK modulator System Object
        bpskModulator = comm.BPSKModulator;
        % change phase offset to pi/16
        bpskModulator.PhaseOffset = 0;
        % create binary data symbols
        data = randi([0 1],1000,M);
        channelBPsk = comm.AWGNChannel('EbNo',20,'BitsPerSymbol',64);
        % modulate the data - for M = 3 sources
        % modData = bpskModulator(data);

        modData_BPSK = zeros(length(data),M);
        for ii = 1:M
            modData_BPSK(:,ii)= step(bpskModulator,data(:,ii));
        end

        channelOutputBpsk = step(channelBPsk,modData_BPSK); %channelPsk(modData);
        x_sig = channelOutputBpsk;

    case 2
        % b) +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        stringCase = 'QPSK';
        % QPSK Modulation Scheme
        % mod = comm.QPSKModulator;
        % refC = constellation(mod);
        % constellation(mod);
        %
        % % Create a QPSK demodulator having a 0 phase offset
        % demod = comm.QPSKDemodulator('PhaseOffset',0);
        % constellation(demod);

        % create a QPSK modulator object and a phase noise object
        qpskModulator = comm.QPSKModulator;
        phNoise = comm.PhaseNoise('Level',-55,'FrequencyOffset',20,'SampleRate',1000);
        channelNoise = comm.AWGNChannel('EbNo',10,'BitsPerSymbol',16);
        d = randi([0 3],1000,M);

        modData_QPSK = zeros(length(d),M);
        for ii = 1:M
            modData_QPSK(:,ii)= step(qpskModulator,d(:,ii));
        end
        %         channelOutputQpsk = step(phNoise,modData_Qpsk);
        channelOutputQpsk = step(channelNoise,modData_QPSK);
        x_sig = channelOutputQpsk;
        %sig= modData_Qpsk;
end

%% ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
% 1.) Specify signal parameters
fc = 10e9;  % Operating frequency, K-band (12-18 GHz), X (8-12 GHz)
c = physconst('LightSpeed');
lambda = c / fc ;

%% ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
% 2.) Define ULA parameters
N_elements = 10;  % number of elements
wavelengthFactor = 1/2;
d_elementSpacing = lambda*wavelengthFactor; % 0.4;  % wavelength, lambda
elementPositions = (0:N_elements-1)*d_elementSpacing;
D_incomingAng = [0 -25 30 12.5];   % truth source locations
M_sources = numel(D_incomingAng);

%% ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
% 3.) Calculate covariance matrix. Calcuate steering vectors
% This will be the input argument to the subspace DOA functions

% Define noise matrix.  The noise snapshot N(w) is a zero-mean random
% vector with spectral matrix, S_n(w) = S_c(w) + sigma_w^2*I
% n = 1:N_elements-1; %(0:299)';
% N_omega = ncov_in*randn(size(n));
ncov_in = db2pow(-25);
S_n_omega = ncov_in*eye(N_elements); % N x N

% Calculate incoming signal covariance matrix. Default value can be 1.
% 1.) incoming signal covariance is a scalar
% 2.) matrix signal covariance matrix; need for noncircular sources
SigCovTypes = {' Scalar incoming SNR (1)', 'Calculate incoming signal covariance matrix (2)'};
str_sigCovTypes = ['Signal  Types:  ']; disp(str_sigCovTypes); disp(SigCovTypes);
promptSigCovTypes = 'User, please indicate which incoming signal covariance type from list above (1 or 2):    ' ;
sigCovUserInput = input(promptSigCovTypes);
tmpSCUI = sigCovUserInput;

if( (isempty(tmpSCUI)) || (tmpSCUI > 2) || (tmpSCUI < 1) || (ischar(tmpSCUI)) )
    disp('User, your input for the "incoming signal covariance type" was invalid.  Please select again (1 or 2).');
    return
end

switch sigCovUserInput
    case 1
        % SNR of incoming plane wave sources
        D_incomingSNR = -5; %0.8; % SNR of incoming sources
        R_sigCov = 10^(D_incomingSNR./10);  %   db2pow(D_incomingSNR);
    case 2
        % Define K snapshots
        K_snapshots = length(x_sig);

        %         R_sigCov = zeros(M_sources,M_sources);
        %         for ii = 1:M
        %             for jj = 1:K_snapshots
        %             R_sigCov(ii)=  (x_sig(jj,ii)*conj(x_sig(jj,ii)'));  % sum over K_snapshots
        %             R_sigCov = R_sigCov + R_sigCov_temp;
        %             clear R_sigCov_temp
        %         end

        R_sigCov = (1/K_snapshots).*x_sig'*conj(x_sig);
end

% Default for R_sigCov
if ~exist("R_sigCov",'var') || isnan(any(R_sigCov,'all')), R_sigCov = 1; end

% function returns "sensor spatial covariance matrix"
%[R_cov,N,D] = sensorSpatialCovarianceMatrix(elementPositions,D_incomingAng,S_n_omega,R_sigCov,fc);
[R_cov,N,D] = sensorSpatialCovarianceMatrix(elementPositions,D_incomingAng,0,1,10e9);

% if isscaler(R_sigCov) && (SigCovTypes == 1)
%     [R_cov,N,D] = sensorSpatialCovarianceMatrix(elementPositions,D_incomingAng,S_n_omega,R_sigCov,fc);
% elseif ismatrix(R_sigCov) && (SigCovTypes == 2)
%     [R_cov,N,D] = sensorSpatialCovarianceMatrix(elementPositions,D_incomingAng,S_n_omega,R_sigCov,fc);
% else
%     R_sigCov = 1;
% end

%% ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
% 4.) Query User input for subspace method type:
% 1.) Unitary TLS-ESPRIT
% 2.) Root MuSiC
subspaceMethodTypes = {' Unitary TLS-ESPRIT (1)', 'Unitary Root MuSiC (2)'};
str_subMethodType = ['Subspace Method Types:  ']; disp(str_subMethodType); disp(subspaceMethodTypes);
promptSubspaceMethodType = 'User, please indicate which subspace method to execute from list above (1 or 2):    ' ;
subspaceMethodsUserInput = input(promptSubspaceMethodType);
tmpSMUI = subspaceMethodsUserInput;

if( (isempty(tmpSMUI)) || (tmpSMUI > 2) || (tmpSMUI < 1) || (ischar(tmpSMUI)) )
    disp('User, your input for the Subspace Method Type was invalid.  Please select again (1 or 2).');
    return
end
% function [DOA_ESPRIT] = tlsUnitaryEsprit(R_cov,D,numElements,d_Spacing,subarrayOffset,rowWeighting)% elemPos,D_incomingAng,noiseAddl,D_incomingSNR,d_elemSpacing,subarrayOffset,rowWeighting)

switch subspaceMethodsUserInput
    case 1
        DOA = tlsUnitaryEsprit(R_cov,M_sources,N_elements,wavelengthFactor);
    case 2
        DOA = rootMusic(R_cov,M_sources,wavelengthFactor); 
end
%% ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
