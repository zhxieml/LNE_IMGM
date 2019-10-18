function [ featBin rho_bins theta_bins ] = makeFeatBin( attrE, nRho_win, nTheta_win )
% make a polar feature representaion using bins
% output type: single

nRho_win = 3;
nTheta_win = 3;

rho_bins = sqrt(2).^([-3:8 Inf]);
theta_bins = [ -pi + pi/18 : pi/9 : pi - pi/18 ]; % notice the last overlap!

rho_win = gausswin(nRho_win);   % odd number of elements
theta_win = gausswin(nTheta_win); % odd number of elements

featBin = computePolarBins( attrE, rho_bins', theta_bins', rho_win, theta_win ); 





