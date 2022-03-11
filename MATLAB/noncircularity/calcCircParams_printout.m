% modData_BPSK and y_64
%qam_64 = y_64;
% sigType = input('User, select signal type?\n 1.) BPSK, 2.) QPSK 3.) PSK, M=8 4.) PSK,M=16 5.) QAM, M=16 6.) QAM, M=64 7.) QAM, M=256, 8.) pammod 9.) aamod 10.) ask\n');
% 
% switch sigType
%     case 1
%         stringCase = 'BPSK';
%         % +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
%         % bpsk modulation scheme
%         % Create a BPSK modulator System Object
%         bpskModulator = comm.BPSKModulator;
%         % change phase offset to pi/16
%         bpskModulator.PhaseOffset = 0;
%         % create binary data symbols
%         data = randi([0 1],100,1);
%         % modulate the data
%         % modData = bpskModulator(data);
%         modData_BPSK = step(bpskModulator,data);
%         channelBPsk = comm.AWGNChannel('EbNo',20,'BitsPerSymbol',64);
%         channelOutputBpsk = step(channelBPsk,modData_BPSK); %channelPsk(modData);
%         
%         %%%%
%         sig = modData_BPSK;
%         
%     case 2
%         stringCase = 'QPSK';
%         % +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
%         % QPSK Modulation Scheme
%         mod = comm.QPSKModulator;
%         refC = constellation(mod);
%         constellation(mod);
%         
%         % Create a QPSK demodulator having a 0 phase offset
%         demod = comm.QPSKDemodulator('PhaseOffset',0);
%         constellation(demod);
%         
%         % create a QPSK modulator object and a phase noise object
%         qpskModulator = comm.QPSKModulator;
%         phNoise = comm.PhaseNoise('Level',-55,'FrequencyOffset',20,'SampleRate',1000);
%         channelNoise = comm.AWGNChannel('EbNo',10,'BitsPerSymbol',16);
%         d = randi([0 3],1000,1);
%         modData_Qpsk = step(qpskModulator,d);
% %         channelOutputQpsk = step(phNoise,modData_Qpsk);
%         channelOutputQpsk = step(channelNoise,modData_Qpsk);
%         sig= modData_Qpsk;
%         
%     case 3
%         stringCase = 'PSK, M=8';
%         % PSK Modulation; AWGN to 8-PSK signal example
%         
%         % Create psk modulator system object; default order is 8
%         pskModulator = comm.PSKModulator;
%         % pskModulator.ModulationOrder=2;  % change order
%         
%         % Modulate the signal
%         tempArr = randi([0 7],2000,1);
%         modData_PSK = step(pskModulator,tempArr);
%         
%         % AWGN to the signal by passing it through a noise channel
%         channelPsk = comm.AWGNChannel('EbNo',20,'BitsPerSymbol',3);
%         
%         % Transmit the signal through the AWGN channel
%         channelOutputPsk = step(channelPsk,modData_PSK); %channelPsk(modData);
%         
%         %%%%%
%         sig = modData_PSK;
%         
%     case 4
%         stringCase = 'PSK, M=16';
%         % 16 PSK example
%         M=16;
%         data = randi([0 M-1],1000,1);
%         modData_PSK_16 = pskmod(data,M,pi/M);
%         channelOutputPSK_16 = awgn(modData_PSK_16,25);
%         %scatterplot(rxSig)
%         sig= modData_PSK_16;
%         
%     case 5
%         stringCase = 'QAM, M=16';
%         % Modulation order of 16 =4^2
%         M_16 = 16;
%         x_16 = (0:M_16-1)';
%         y_16 = qammod(x_16,M_16,'bin','UnitAveragePower',true,'PlotConstellation',true);
%         %scatterplot(y_16); grid on;
%         sig = y_16;
%         
%     case 6
%         stringCase = 'QAM, M=64';
%         % M = 64 4^3
%         M_64 = 64;
%         x_64 = (0:M_64-1)';
%         y_64 = qammod(x_64,M_64,'bin','UnitAveragePower',true,'PlotConstellation',true);
%         %scatterplot(y_64); grid on;
%         sig = y_64;
%         
%     case 7
%         stringCase = 'QAM, M=256';
%         % M = 256 = 4^4
%         M_256 = 256;
%         x_256 = (0:M_256-1)';
%         y_256 = qammod(x_256,M_256,'bin','UnitAveragePower',true,'PlotConstellation',true);
%         %scatterplot(y_256); grid on;
%         sig = y_256;
%         
%     case 8
%         stringCase = 'PAM, M=8';
%         M = 8;
%         data = randi([0 M-1],100,1);
%         modData_PAM = pammod(data,M,pi/8);
%         sig = modData_PAM;
%         
%     case 9
%         stringCase = 'AAMOD';
%         fs = 100;
%         t = (0:1/fs:100)';
%         fc = 10;
%         x = sin(2*pi*t);
%         ydouble = ammod(x,fc,fs);
%         ysingle = ssbmod(x,fc,fs);
%         sa = dsp.SpectrumAnalyzer('SampleRate',fs, ...
%     'PlotAsTwoSidedSpectrum',false, ...
%     'YLimits',[-60 40]);
%         step(sa,ysingle);
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

% 1.) For a complex random variable, z = x + iy, calculate the variance,
% E[|z|^2] = sigma_x^2 + sigma_y^2
xi = real(sig);
yi = imag(sig);

fprintf('%++++++++++++++++++++++++++++++%')
str = [sprintf('\n Printing out signal properties for: %s ', stringCase)]; disp(str);

avg_xi = sum(xi)./numel(xi);
sigma_xi = sum((xi-avg_xi).^2)./(numel(xi)-1);

avg_yi = sum(yi)./numel(yi);
sigma_yi = sum((yi-avg_yi).^2)./(numel(yi)-1);

fprintf('\n sigma_xi^2 = %0.8f ',sigma_xi);
fprintf('\n sigma_yi^2 = %0.8f ',sigma_yi);

% 2.) calculate sigma_xy, which is the joint probability of two RVs
% cov(X,Y) = cov(Y,X) = sum(X-xbar)(Y-ybar)/n-1
% evident that sigma_xy = sigma_yx
sigma_xyi = sum((xi-avg_xi).*(yi-avg_yi))./(numel(xi)-1);
sigma_yxi = sum((yi-avg_yi).*(xi-avg_xi))./(numel(yi)-1);

fprintf('\n sigma_yxi^2 = %0.8f ',sigma_yxi);
fprintf('\n sigma_xyi^2 = %0.8f ',sigma_xyi);

% 3.)calculate pseudo-variance: \tau_z = sigma_x^2 - sigma_y^2 + j2sigma_xy
tau_zi = sigma_xi.^2 - sigma_yi.^2 + 1j.*2.*sigma_xyi;

fprintf('\n tau_zi^2 =  %0.8f ',tau_zi);

% Note: A circular r.va. z, has the property that its pseudo-variance
% vanishes, tau_z = 0  (i.e. sigma_x^2 = sigma_y^2 and sigma_xy = 0) PROPER

% 4.) form the covariance matrix:
cov_rvzi = [sigma_xi.^2 sigma_xyi; sigma_yxi sigma_yi.^2];

fprintf('\n cov_rvzi^2 = %0.8f ',cov_rvzi);

% 5.) eigenvector decomposition of covariance matrix
eigVals = eig(cov_rvzi);
lambda_1 = max(eigVals);
lambda_2 = min(eigVals);

fprintf('\n lambda_1 = %0.8f ',lambda_1);
fprintf('\n lambda_2 = %0.8f ',lambda_2);

% 6.) define the first eigenvector to have positive first coordinate
% If lambda_1 > lambda_2, then the EVD is unique.  If lambda_1 = lambda_2,
% then sigma_x = sigma_y and sigma_xy = 0, then lambda cannot be determined
% and is arbitrary.

% 7.) calculate eccentricity. Note that epsilon = 0, denotes a circle
epsilon = sqrt( (lambda_1-lambda_2) ./ (lambda_1+lambda_2) );
fprintf('\n epsilon = %0.8f ',epsilon);

% 8.) geometric mean of the eigenvalues and mean of eigenvalues, q <= 1
q = (sqrt(lambda_1.*lambda_2)) ./ (0.5*(lambda_1 + lambda_2));
fprintf('\n q = %0.8f ',q);

% 9.) variance of sigma_z = lambda_1 + lambda_2 measures the scale of ellipse
% complex covariance sigma_z^2 = cov(z,z), tau_z = cov(z,z*)
sigma_zi = sqrt(lambda_1 + lambda_2);
fprintf('\n sigma_zi = %0.8f  ',sigma_zi);

% 10.) circularity quotient is defined as the quotient between the
% pseudo-variance and the variance
% Its unique polar representation rho_z = r_zexp(jtheta) induces quantities
% r_z |rho_z| as the circularity coefficient of z and theta = arg[rho_z]
% which is the circularity angle of z. Note
% arg(rho_z)=atan(imag(rho_z)./real(rho_z))
circularityQuotient = tau_zi ./ sigma_zi.^2;
fprintf('\n (circularity quotient) rho_z = %0.8f ',circularityQuotient);

% circularity coefficient
rho_z = circularityQuotient;
r_z = sqrt(real(rho_z).^2 + imag(rho_z).^2);

epsilon_rz = sqrt(abs(r_z));  %sqrt( (lambda_1-lambda_2) ./ (lambda_1+lambda_2) );
fprintf('\n epsilon_rz = %0.8f ',epsilon_rz);

alphaOrientation = atan2(imag(circularityQuotient),real(circularityQuotient));
fprintf('\n alphaOrientation = %0.8f ',alphaOrientation);

% link circularity quotient with the correlation coefficient
rho_corr = sigma_xyi ./ (sigma_xi.*sigma_yi);
fprintf('\n rho_corr = %0.8f ',rho_corr);

% 11.) hypothesis testing H0 : C = 0 is proper; H1 : C neq 0 is improper
% statistical hypothesis test of circularity of the sample {z=x_i+jyi}
% is equivalent with the test of sphericity of the composite sample
% {vi=(xi,yi)^T}; tau_z = 0 or tau_z neq 0.

% If the covariance matrix vanishes, it is called proper, C_zz = 0
% Complex variable z is represented as a vector [z_R,z_I]^T composed of its
% real and imaginary parts and described by a bivariate CDF: F(z_R,z_I);
% F(z) = F(z_R,z_I), with the real bivariate PDF defined as f(z_R,z_I) =
% d^2F(z_r,z_I)./dzzrdzi; PDF f(z) = f(z_R,z_I)

% 256 x 2Sigmass
vi = [xi yi]';

% calculate Sigma_hat
numVi = size(vi,2);
Sigma_hat = zeros(2,2,numVi);
for i = 1:numVi
    Sigma_hat(:,:,i) = (vi(:,i)*vi(:,i)'); %disp(Sigma_hat);
end

% sum over columns
sumTerm = sum(Sigma_hat,3)./numVi;

% check the same values as vi*vi^T
testSumTerm = (vi*vi')./numVi;

ln = (sqrt(det(sumTerm))./(0.5*trace(sumTerm))).^numVi;

% if H0 is true, then -2nlogq --> xhi^2
% ln = (1-r_z)^(n/2)
% H0 : tau_z  = 0 --> Proper
% H1 : tau_z != 0 --> Improper, and testStat_H0 exceeds tolerance_chi2
testStat_H0 = -numVi.*log10(ln.^(2/numVi));
chiSquared_quantile = 5.991;  % tolerance

testStatistic = testStat_H0;
fprintf('\n test statistic, testStat = %0.8f \n',testStatistic);
%atan2(imag(rho_z),real(rho_z))

fprintf('\n %0.8f ',sigma_xi);
fprintf('\n %0.8f ',sigma_yi);
fprintf('\n %0.8f ',sigma_xyi);
fprintf('\n %0.8f ',sigma_yxi);
fprintf('\n %0.8f ',real(tau_zi));
fprintf('\n %0.8f ',imag(tau_zi));
fprintf('\n %0.8f ',cov_rvzi(1,1));
fprintf('\n %0.8f ',cov_rvzi(1,2));
fprintf('\n %0.8f ',cov_rvzi(2,1));
fprintf('\n %0.8f ',cov_rvzi(2,2));
fprintf('\n %0.8f ',lambda_1);
fprintf('\n %0.8f ',lambda_2);
fprintf('\n %0.8f ',(epsilon));
fprintf('\n %0.8f ',q);
fprintf('\n %0.8f ',sigma_zi);
fprintf('\n %0.8f ',rho_z);
fprintf('\n %0.8f ',epsilon_rz);
fprintf('\n %0.8f ',alphaOrientation);
fprintf('\n %0.8f ',rho_corr);
fprintf('\n %0.8f ',ln);
fprintf('\n %0.8f \n' ,testStatistic);
