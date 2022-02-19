function [sV, N, D] = steeringVector(pos_inp,D_srcAng_inp,fc)
%% +++++++++++++++++++++++++++++++++++++++++++++++++++++
% a.) pos_inp - element positions in array; diff(pos_inp) 
% yields the interelement spacing
% b.) D_srcAng_inp - the angle(s) of the plane wave propagating
% into the arrays, specified in spherical AZ and EL coordinates
% c.) fc:  "Operating frequency".  In other words, this is the
% frequency of the wave(s) propagating into the array


%% +++++++++++++++++++++++++++++++++++++++++++++++++++++

% Define constants and additional parameters
% fc = 12e9;                               % Operating frequency

if fc ~= 1
    c = physconst('Lightspeed');
else
    c = 1;
end

% Parse the inputs of the arguments, ensuring proper format 
[pos, D_srcAng,N_elem, D_src] = parseInputsSteerVec(pos_inp,D_srcAng_inp);

% Declare positions of incoming sources in Az and El, [deg]
Az = D_srcAng(:,1);
El = D_srcAng(:,2);

% Assign position angles
azang = Az;
elang = El;

% Calculate incident angles; this is direction cosines
% angles defined in local coordinate system
incidentdir = [-cosd(elang).*cosd(azang);...
    -cosd(elang).*sind(azang);...
    -sind(elang)];

% calculate delay, tau [sec], between the elements in the array
tau = pos.'*incidentdir/c;  % N x D, for D plane waves

% Need to add in quantized phase portion for frequency array in 
% steering vector calculation.  For now, just leave freq as a scalar.
freq = fc;
sV = exp(-1i*2*pi*freq*tau);  % 10 x 1, N x D, for D plane-waves
N = N_elem; D = D_src;

end

% [EOF]
