function [w_MVDR, arrayPattern_MVDR, N, D] = mvdrBeamform(elemPos,angIncoming,noiseAddl,sigCov)

if nargin < 3, noiseAddl = 0; end
if nargin < 4, sigCov = 1; end

% ncov_in = db2pow(-10);
ncov_in = db2pow(noiseAddl);
[R_cov,N,D] = sensorSpatialCovarianceMatrix(elemPos,angIncoming,ncov_in,sigCov);

% Constraint matrix specified as a complex-valued, N-by-K, complex-valued matrix.
% In this matrix N represents the number of elements in the sensor array while K
% represents the number of constraints. Each column of the matrix specifies a
% constraint on the beamformer weights. The number of K must be less than or equal to N.
constraintMatrix = steeringVector(elemPos,angIncoming,1);
convinvAS = R_cov\constraintMatrix;
w_MVDR = convinvAS ./ sum(conj(constraintMatrix).*convinvAS);

% Extent of angles over which to plot over for array pattern
plotAngles = -90:0.1:90;
steerVecPlt_MVDR = steeringVector(elemPos,plotAngles,1);

arrayPattern_MVDR = w_MVDR'*steerVecPlt_MVDR;

end