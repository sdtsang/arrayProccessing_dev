% modData_BPSK and y_64
% qam_64 = y_64;
function [testStatistic] =  calcCircParams(sig) 
% 1.) For a complex random variable, z = x + iy, calculate the variance,
% E[|z|^2] = sigma_x^2 + sigma_y^2
xi = real(sig);
yi = imag(sig);

avg_xi = sum(xi)./numel(xi);
sigma_xi = sum((xi-avg_xi).^2)./(numel(xi)-1);

avg_yi = sum(yi)./numel(yi);
sigma_yi = sum((yi-avg_yi).^2)./(numel(yi)-1);

% 2.) calculate sigma_xy, which is the joint probability of two RVs
% cov(X,Y) = cov(Y,X) = sum(X-xbar)(Y-ybar)/n-1
% evident that sigma_xy = sigma_yx
sigma_xyi = sum((xi-avg_xi).*(yi-avg_yi))./(numel(xi)-1);
sigma_yxi = sum((yi-avg_yi).*(xi-avg_xi))./(numel(yi)-1);

% 3.)calculate pseudo-variance: \tau_z = sigma_x^2 - sigma_y^2 + j2sigma_xy
tau_zi = sigma_xi.^2 - sigma_yi.^2 + 1j.*2.*sigma_xyi;

% Note: A circular r.va. z, has the property that its pseudo-variance 
% vanishes, tau_z = 0  (i.e. sigma_x^2 = sigma_y^2 and sigma_xy = 0) PROPER

% 4.) form the covariance matrix:
cov_rvzi = [sigma_xi.^2 sigma_xyi; sigma_yxi sigma_yi.^2];

% 5.) eigenvector decomposition of covariance matrix
eigVals = eig(cov_rvzi);
lambda_1 = eigVals(1);
lambda_2 = eigVals(2);

% 6.) define the first eigenvector to have positive first coordinate
% If lambda_1 > lambda_2, then the EVD is unique.  If lambda_1 = lambda_2,
% then sigma_x = sigma_y and sigma_xy = 0, then lambda cannot be determined
% and is arbitrary.  

% 7.) calculate eccentricity. Note that epsilon = 0, denotes a circle
epsilon = sqrt( (lambda_1-lambda_2) ./ (lambda_1+lambda_2) );

% 8.) geometric mean of the eigenvalues and mean of eigenvalues, q <= 1
q = (sqrt(lambda_1.*lambda_2)) ./ (0.5*(lambda_1 + lambda_2));

% 9.) variance of sigma_z = lambda_1 + lambda_2 measures the scale of ellipse
% complex covariance sigma_z^2 = cov(z,z), tau_z = cov(z,z*)
sigma_zi = sqrt(lambda_1 + lambda_2);

% 10.) circularity quotient is defined as the quotient between the
% pseudo-variance and the variance
% Its unique polar representation rho_z = r_zexp(jtheta) induces quantities
% r_z |rho_z| as the circularity coefficient of z and theta = arg[rho_z]
% which is the circularity angle of z. Note
% arg(rho_z)=atan(imag(rho_z)./real(rho_z))
circularityQuotient = tau_zi ./ sigma_zi.^2;
alphaOrientation = atan2(imag(circularityQuotient),real(circularityQuotient));

% link circularity quotient with the correlation coefficient
rho_corr = sigma_xyi ./ (sigma_xi.*sigma_yi);

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
end