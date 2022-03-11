function DOA_rootMuSiC =  rootMuSiC(R,D,lambdaFactor) % elemPos,D_incomingAng,noiseAddl,D_incomingSNR)

elSpacing = lambdaFactor;
Msig = D;
%%++++++++++++++++++++++++++++++++++
% Check for positive semi definite
[eigenvects,sEDArg] = eig((R+R')/2);  % ensure Hermitian
sED = diag(real(sEDArg));
diagEigenVals =sED;
tol = eps(max(abs(sED))); % based on stats cholcov
sED(abs(sED)<=tol)=0;
cond = any(sED<0);
if cond
    coder.internal.errorIf(cond,...
        'phased:phased:notPositiveSemiDefinite','R');
end

%%++++++++++++++++++++++++++++++++++
N = size(R,1);
%Sort eigenvectors
[~,indx1] = sort(diagEigenVals,'descend');
eigenvects = eigenvects(:,indx1);
% Separate the signal and noise eigenvectors
noise_eigenvects = eigenvects(:,Msig+1:end);

% Form a polynomial D
% D consists of a sum of polynomials given by the product of the noise
% subspace eigenvectors and the reversed and conjugated version.
D = complex(zeros(2*N-1,1));
for i = 1:N-Msig
    D = D + conv(noise_eigenvects(:,i),conj(flipud(noise_eigenvects(:,i))));
end

% Take the angle of the MSIG roots of D closest to the unit circle.
roots_D = roots(D);
roots_D1 = roots_D(abs(roots_D) < 1);
[~,indx] = sort(abs(abs(roots_D1)-1));
sorted_roots = roots_D1(indx);
doas = angle(sorted_roots(1:Msig));

u = doas/(2*pi*elSpacing);
% check whether all elements of u are within [-1,1]
idx = find(abs(u)<=1);
if  length(idx) <Msig && isempty(coder.target)
    warning(message('phased:phased:internal:AbstractULASubspaceDOA:InvalidPsi',Msig));
end
if isempty(idx)
    ang = zeros(1,0);
else
    ang = asind(u(idx));
    % convert to row vector
    ang = ang(:).';
end

DOA_rootMuSiC = ang; str_rootmusic = ['DOA (root-MuSiC) =  ' num2str(DOA_rootMuSiC,'%0.6f')];  disp(str_rootmusic);


end
