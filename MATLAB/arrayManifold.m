function [v_kz, v_psi] = arrayManifold(fc,N_elem)
%% +++++++++++++++++++++++++++++++++++++++++++++
% Calculates array manifold vector following Van Trees formulation
% Returns vector (ULA) or matrix (planar) v_k and v_psi

% For vector case, the 0th sensor is assumed to be the reference
% For matrix case, the center sensor in the planar array is
% assumed to be the reference

% a.) fc: frequency of propagating wave

%% ++++++++++++++++++++++++++++++++++++++++++++++

% Define constants and additional parameters
c = physconst('Lightspeed');
lambda = c/fc ;
d = lambda/2;
theta_vec = -2*pi:0.1:2*pi;
% omega = 2*pi*fc; %-2*pi:0.1:2*pi;
kz = (2*pi./lambda)*cos(theta_vec);
d_interelement = d;  % interelement spacing of sensors, 10 total
psi_array = -kz*d; % -kz*d=(2*pi/lamba)*cos(theta)*d=(2*pi/lambda)*uz*d, uz = cos(theta)

% Place 0th sensor as the reference sensor
N = N_elem; % number of elements in array
sensor_dist = zeros(1,N);
sensor_dist(1) = 0;
for i = 2:N
    sensor_dist(i) = sensor_dist(i-1)+d;
end

% Calculate the locations of the elements
p0 = (0-((N-1)/2))*d_interelement;
p_1toN = zeros(N-1,1);
for ii = 1:N-1
    p_1toN(ii) = (ii-((N-1)/2))*d_interelement;
end

% Concatenate first element of position vector
p_zn = [p0; p_1toN];

% Array manifold vector for a ULA exhibits conjugate symmetry (2.73).
% This form emphasizes the Vandermonde structure (A.163) of v_psi(psi)
% calculate array manifold vector

% Pre-allocate memory
v_kz  = zeros(length(kz),N);
v_psi = zeros(length(psi_array),N);

% Calculate array manifold vector; columns in array correspond to each
% element
for jj = 1:N
    v_kz(:,jj)  = exp( -1j*kz'.*p_zn(jj) );  % 126 x 10 - example
    v_psi(:,jj) = exp( -1j*psi_array.*(p_zn(jj)/d_interelement) );  %  126 x 10 - example
end


end