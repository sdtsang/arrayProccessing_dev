function [DOA_ESPRIT] = tlsUnitaryEsprit(elemPos,D_incomingAng,noiseAddl,D_incomingSNR,d_elemSpacing,subarrayOffset,rowWeighting)
% 1.) Three Signals Arriving at HALF-wavelength-Spaced ULA
% Assume a half-wavelength spaced uniform line array with 10 elements.
% Three plane waves arrive from the 0°, –25°, and 30° azimuth directions.
% Elevation angles are 0°. The noise is spatially and temporally white.
% The SNR for each signal is 5 dB. Find the arrival angles.

if nargin < 3, noiseAddl = 0; end
if nargin < 4, D_incomingSNR = 1; end
if nargin < 5, d_elemSpacing = 0.5; end
if nargin < 6, subarrayOffset = 1; end
if nargin < 7, rowWeighting = 1; end

% SNR of incoming plane wave sources
SNR_pw = log10(D_incomingSNR./10);  %   db2pow(D_incomingSNR);
% ncov_in = db2pow(-10);
ncov_in = db2pow(noiseAddl);
[R_cov,N,D] = sensorSpatialCovarianceMatrix(elemPos,D_incomingAng,ncov_in,SNR_pw);

% [V,D] = eig(A,B) produces a diagonal matrix D of generalized
% eigenvalues and a full matrix V whose columns are the corresponding
% eigenvectors so that A*V = B*V*D.
[eigenvects,sED] = eig((R_cov+R_cov')/2);  % ensure Hermitian

% Check for positive semi definite
sED = diag(sED);
tol = eps(max(abs(sED))); % based on stats cholcov
sED(abs(sED)<=tol)=0;
cond = any(sED<0);
if cond
    disp('User, not positive semi-definite, exiting now.');
    disp('Check input parameters and covariance matrix');
    return;
end

diagEigenVals = sED;

% Sort eigenvectors
[~,indx] = sort(diagEigenVals,'descend');
eigenvects = eigenvects(:,indx);

% Check eigenvalues against source dimension. ESPRIT does not work
% when there are less then D non-zero eigenvalues.
D_sumDiagEig = sum(diagEigenVals>eps(max(abs(diagEigenVals))));

if D_sumDiagEig < D
    str_checkEig = ['User, check rank.  It seems there are less than D (src) non-zero eigenvalues.']; disp(str_checkEig);
    sprintf('(D_sumDiagEig = %s) < (D_srcNum = %s)',num2str(D_sumDiag),num2str(D));
    return;
end

elSpacing = d_elemSpacing;
saSpacing = subarrayOffset;
rweight   = rowWeighting;

% Row weighting
Ns = N-saSpacing; %number of elements in a subarray
ms = rweight;
w = min(ms,Ns-ms+1);                             % Eq 9.133 in [1]
weights = diag(sqrt([1:w-1 w*ones(1,Ns-2*(w-1)) w-1:-1:1])); % Eq 9.132 in [1]
O = zeros(Ns,saSpacing);

% Selection Matrices
Js1 = [weights O]; % Eq 9.134 in [1]
Js2 = [O weights]; % Eq 9.135 in [1]

% Selecting subarray signal subspaces
Us1 = Js1*eigenvects(:,1:D);
Us2 = Js2*eigenvects(:,1:D);

% TLS-ESPRIT
C = [Us1';Us2']*[Us1 Us2];    % Eq. (9.123) in Optimum Array Processing [1]
[U,~,~] = svd(C);             % C is 2*D x 2*D
V12 = U(1:D,D+1:2*D);         % D x D
V22 = U(D+1:2*D,D+1:2*D);     % D x D
psi = -V12/V22;               % Eq. (9.122) in Optimum Array Processing [1]
psieig = eig(psi);

%   Extract angle information estimated from two subarrays based on the
%   distance of the phase center between subarrays.
doas = 1/saSpacing*angle(psieig);

%Convert estimated angle in sin-space to degrees. This method is valid for
%ULA only.

u = doas/(2*pi*elSpacing);

% check whether all elements of u are within [-1,1]
idx = find(abs(u)<=1);
if  length(idx) < D
    warning(message('User, invalid psi; check input parameters ...',D));
end

if isempty(idx)
    ang = zeros(1,0);
else
    ang = asind(u(idx));
    % convert to row vector
    ang = ang(:).';
end

    DOA_ESPRIT = ang; str_esprit = ['DOA (ESPRIT) =  ' num2str(ang)];  disp(str_esprit);

end