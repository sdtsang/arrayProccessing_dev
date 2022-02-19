function [w_LCMV, arrayPattern_LCMV, N, D] = lcmvBeamform(elemPos,angIncoming,resp,noiseAddl,sigCov)
%   Calculates the linea constraint minimum variance (LCMV) weights using sample matrix
%   inversion (SMI) method. X is an MxN data matrix whose rows are
%   snapshots in time. W is a length N column vector representing the LCMV
%   weights. C is an NxL matrix whose columns are constraints and L stands
%   for the number of constraints. Note that L cannot exceed N. G is a
%   length L column vector whose elements are the desired responses
%   corresponding to the constraints specified in C. DL is the diagonal
%   loading factor for the data covariance matrix. DL must be a
%   non-negative scalar.

if nargin < 3, noiseAddl = 0; end
if nargin < 4, sigCov = 1; end

%   Desired response specified as complex-valued, K-by-1 column vector
%   where K is the number of constraints. The value of each element in
%   the vector is the desired response to the constraint specified in the
%   corresponding column of constr.
C = resp;

% ncov_in = db2pow(-10);
ncov_in = db2pow(noiseAddl);
[R_cov,N,D] = sensorSpatialCovarianceMatrix(elemPos,angIncoming,ncov_in,sigCov);
G = R_cov;

sV_lcmv = steeringVector(elemePos,angIncoming,1);

[m,n] = size(sV_lcmv);
sV_lcmv = sV_lcmv/sqrt(m);

% Add diagonal loading+
delta = sqrt(dl);
if delta ~= 0
    x = [sV_lcmv;delta*eye(n)];
else
    x = sV_lcmv;
end

% The covariance matrix is defined as E{x.'*conj(x)}
if size(C,2) == 1
    % MVDR
    temp = qrlinsolve(x.',C);
    w_LCMV = G*temp/(C'*temp);
else
    % LCMV
    if m >= n
        [temp,F] = qrlinsolve(x.',C);
        w_LCMV = temp*qrlinsolve(F',G);
    else
        % when matrix is fat, F is no longer square and we cannot play the
        % trick of thin matrix. Therefore, we have to form R2 and use LU.
        temp = qrlinsolve(x.',C);
        R2 = C'*temp;   % R2 = C'*R^(-1)*C
        [L2, U2] = lu(R2);
        temp2 = U2\(L2\G);
        w_LCMV = temp*temp2;
    end
end

% Extent of angles over which to plot over for array pattern
plotAngles = -90:0.1:90;  % assume ULA for now
steerVecPlt_LCMV = steeringVector(elemPos,plotAngles,1);

arrayPattern_LCMV = w_MVDR'*steerVecPlt_LCMV;

end