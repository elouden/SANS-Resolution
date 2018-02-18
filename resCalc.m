function [dqR dpA dwL] = resCalc(B, lam, dlam_lam, w1, w2, L) 
% 10/31/2017 - E R Louden 
%
% function [dqR dpA dwL] = resCalc(B, lam, dlam_lam, w1, w2, L) 
% RESCAL:  Calculates the SANS experimental resolution in DEGREES (for dpA and dwL)
%          Note, this script requires the FINAL pinhole separation
%          make sure you adjust for any offsets based on the particular instrument
%   
%   dqR         -    radial resolution
%   dpA         -    azimuthal resolution
%   dwL         -    longitudinal resolution  
%
%   B           -    magnetic field (T)
%   lam         -    wavlength (AA)
%   dlam_lam    -    wavlength resolution (percent)
%   w1          -    source aperture
%   w2          -    sample aperture
%   L           -    pinhole separation

% These expressions are based on those found in M R Eskildsen's Thesis

% q-value, determined from B
    fluxQ = 2070; %nm^2 * T,  2.067833831*10^(-15) Wb = m^2 * T
    q_nm = 2*pi*sqrt((2*(B)) / (sqrt(3)*fluxQ)); % 1/nm
    q = q_nm/10; % conver to inv ang

% two main components of resolution: 
    % wavelength spread
    d2Theta_WS = ((q*lam)/(4*pi)) * sqrt((4*log(2))/3) * (dlam_lam);
    % beam divergence
    d2Theta_BD = 0.5 * sqrt(2*log(2)) * sqrt((w1/L)^2 + (w2/L)^2);

% resolution in each direction
    dqR = sqrt( (((2*pi)/lam)*d2Theta_WS)^2 + (((2*pi)/lam)*d2Theta_BD)^2 );
    dpA = ((2*pi)/lam)/q * d2Theta_BD;
    dwL = sqrt(d2Theta_WS^2 + d2Theta_BD^2);

% convert dpA and dwL to degrees
    dpA = dpA * (180/pi);
    dwL = dwL * (180/pi);
end