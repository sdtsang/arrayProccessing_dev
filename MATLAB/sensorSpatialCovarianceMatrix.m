function [R_cov,N,D] = sensorSpatialCovarianceMatrix(elemPos,angIncoming,noiseAddl,sigCov,fc)

if nargin < 3, noiseAddl = 0; end
if nargin < 4, sigCov = 1; end

pos_inp = elemPos; D_srcAng_inp = angIncoming;

[sV, N, D] = steeringVector(pos_inp,D_srcAng_inp,fc);

%N_omega = noiseAddl*randn(length(elemPos));
S_n_omega = noiseAddl*eye(length(elemPos));

%   sigCov is either a positive scalar, a 1xM vector, or
%   an MxM positive semidefinite matrix. If sigCov is a scalar, it represents
%   the power (in watts) of the incoming signals. All incoming signals are
%   assumed to be uncorrelated and share the same power level. If sigCov is a
%   1xM matrix, then it represents the power of individual incoming
%   signals. However, all incoming signals are still uncorrelated to each
%   other. If sigCov is an MxM matrix, then it represents the covariance
%   matrix for all incoming signals. The default value of sigCov is 1.

% sigCov needs to be an Mx M matrix to represent covariance matrix for
% incoming signals and to extend to noncircular
% sV and sV' are the N x M and M x N steering vectors, respectively.
Rx = sV*sigCov*sV' + S_n_omega;
Rx = (Rx+Rx')/2;  % ensure Hermitian

R_cov = Rx;

end