function [residual,GFD,A] = all_pole_speech_analysis(sig,p,N)
% Purpose: 
%       This function computes the residual and GFD of the entire speech
%       signal via the autocorrelation method combined with a second-order
%       pre-emphasis filter. It also gives the corresponding lpc vector for 
%       each frame.
%
% Arguments:
%       1) sig [Lsig x 1]: The speech signal to be analyzed.
%       2) p [1 x 1]: Linear prediction order.
%       3) N [1 x 1]: Number of samples per frame.
% Return:
%       1) residual [Shift*Nfr x 1]: It is the total residual,
%       where Shift and Nfr are defined in the code. Approximetely
%       Shift*Nfr = Lsig.
%       2) GFD [Shift*Nfr x 1] : It is the GFD estimate.
%       3) A [p+1 x Nfr]: It is the matrix containing the lpc vectors of
%                         all frames.
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

signal = sig(:);

Win = hanning(N);              % analysis window.
Shift = N/2;                   % OLA step.
Lsig = length(signal);         % signal length.
Nfr = floor(Lsig/Shift);       % number of frames.
slice = 1:N;
Buffer_re = zeros(Shift,1);    % it is required for the OLA procedure
Buffer_GFD = zeros(Shift,1);
tosave = 1:Shift;
residual = zeros(Shift*Nfr,1);
GFD = zeros(Shift*Nfr,1);


residual_frame = 0;
GFD_frame = 0;
pre_emph = 0.99;
A = zeros(p+1,Nfr);            % lpc matrix
signal = [signal; zeros(Shift,1)];


% speech analysis frame by frame.
for l=1:Nfr
    
  sig_slice = signal(slice);
  sig_slice = sig_slice - mean(sig_slice);
  
  x_frame_no_preEmphasis = Win.*sig_slice;
  
%  second-order pre-emphasis
  x_frame_withoutWindow = filter([1 -pre_emph],1,filter([1 -pre_emph],1,sig_slice));
  x_frame = Win.*x_frame_withoutWindow;
  

  x_frame = x_frame(:); 
  
  a = lpc(x_frame,p);  
  
  A(1:p+1,l) = a(1:p+1);
  residual_frame = filter(a,1,x_frame);  
  
  % inverse filtering with the non-pre-emphasized signal
  % (Estimation of Glottal flow derivative)
  GFD_frame = filter(a(1:p+1),1,x_frame_no_preEmphasis);
  
  
  % Overlap and add
  residual_frame(1:Shift) = residual_frame(1:Shift) + Buffer_re;
  GFD_frame(1:Shift) = GFD_frame(1:Shift) + Buffer_GFD;
  residual(tosave) = residual_frame(1:Shift);
  GFD(tosave) = GFD_frame(1:Shift);
  Buffer_re = residual_frame(Shift+1:N);
  Buffer_GFD = GFD_frame(Shift+1:N);

  
  slice = slice+Shift;      % move to next frame
  tosave = tosave+Shift;
  
end

end