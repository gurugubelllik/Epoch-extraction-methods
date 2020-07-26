function [n,sigma] = GaussianWhiteNoise(signal,SNR)
% Purpose: 
%          Generates white Gaussian noise (WGN) with a specific SNR.
%
% Arguments:
%           1) signal: the clean signal (should be zero mean!!!!).
%           2) SNR: signal to noise ratio (in dB).
% Return:  
%           1) n : the noise vector which has the same length with signal.
%           2) sigma: standard deviation of noise.
%
% Reference: 
%       1. A.I. Koutrouvelis, G.P. Kafentzis, N.D. Gaubitch, R. Heusdens,
%          "A Fast Method for High-Resolution Voiced/Unvoiced Detection
%          and Glottal Closure/Opening Instant Estimation of Speech"
%
% Last modified: 11/07/2015
%
% Author: Andreas Koutrouvelis
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Copyright 2015: Delft University of Technology. The software is free for
% non-commercial use. This program comes WITHOUT ANY WARRANTY.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
noise_energy = ((norm(signal)^2)/(10^(SNR/10)));
noise_power = noise_energy/length(signal);
sigma = sqrt(noise_power);
n = sigma*randn(length(signal),1);