% demo1: A simple demonstration of GEFBBA.
%
% Reference: 
%       1. A.I. Koutrouvelis, G.P. Kafentzis, N.D. Gaubitch, R. Heusdens,
%          "A Fast Method for High-Resolution Voiced/Unvoiced Detection
%          and Glottal Closure/Opening Instant Estimation of Speech"
%
% Last modified: 25/01/2018
%
% Author: Andreas Koutrouvelis
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Copyright 2015: Delft University of Technology. The software is free for
% non-commercial use. This program comes WITHOUT ANY WARRANTY.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Speech from APWALD
[sig,fs] = audioread('H.22.16k.wav');
sig = sig - mean(sig);
start_sig = sig;



WGN_select = 0; % 1: with WGN, 0 without WGN.
if(WGN_select)
    SNR = 20;   % Specify the noise SNR.
    [nnnn,sigmaaaa] = GaussianWhiteNoise(sig,SNR);
    sig = sig + nnnn;
end




%%%%%%%%%%%%%%%%%%%%%%%%%%%% Run VGEFBA %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
MinP = 80;
tic;
[EstimatedGOIs, EstimatedGCIs, EstimatedEes,...
 EstimatedFZCIs, EstGFD, EstRes, A] = GEFBA(sig,fs);

% find the boundaries of the voiced segments of GEFBA
[startssGEFBA,finishessGEFBA] =...
    voiced_segments(EstimatedGCIs, EstimatedGOIs,fs,MinP);
GEFBA_speed = toc





%%%%%%%%%%%%%%%%%%%%%%%%% PLOTTING %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% subplot(3,1,1) contains the following information:
%        a) the speech signal.
%        b) the estimated voiced segments boundaries with black lines.
%        c) the estimated GCIs and GOIs with black '+' and red 'o'.
%
% subplot(3,1,2) contains the following information:
%        a) the residual.
%        b) the estimated voiced segments boundaries with black lines.
%        c) the estimated GCIs and GOIs with black '+' and red 'o'.
%
% subplot(3,1,3) contains the following information:
%        a) the GFD signal.
%        b) the estimated voiced segments boundaries with black lines.
%        c) the estimated GCIs and GOIs with black '+' and red 'o'.

figure(1);
skip_Samples = 500;
xaxis_start = 1;
xaxis_end = length(sig)-skip_Samples;
xaxis_ms = xaxis_start:xaxis_end;


subplot(3,1,1);
plot(sig(xaxis_start:xaxis_end)./max(abs(sig(xaxis_start:xaxis_end))),'k');
hold on;
plot(EstimatedGOIs(EstimatedGOIs>=xaxis_start & EstimatedGOIs<=xaxis_end)-(xaxis_start-1),0,'c--o','MarkerSize',7);
hold on;
plot(EstimatedGCIs(EstimatedGCIs>=xaxis_start & EstimatedGCIs<=xaxis_end)-(xaxis_start-1),0,'c--+','MarkerSize',7);
hold on;
hold on;
for i=1:length(startssGEFBA)
    plot(ones(2,1).*startssGEFBA(i),[0,1],'k')
    hold on;
end
for i=1:length(finishessGEFBA)
    plot(ones(2,1).*finishessGEFBA(i),[0,1],'k')
    hold on;
end
ylabel('speech');


subplot(3,1,2);
plot(EstRes(xaxis_start:xaxis_end)./max(abs(EstRes(xaxis_start:xaxis_end))),'k');
hold on;
plot(EstimatedGOIs(EstimatedGOIs>=xaxis_start & EstimatedGOIs<=xaxis_end)-(xaxis_start-1),0,'c--o','MarkerSize',7);
hold on;
plot(EstimatedGCIs(EstimatedGCIs>=xaxis_start & EstimatedGCIs<=xaxis_end)-(xaxis_start-1),0,'c--+','MarkerSize',7);
hold on;
for i=1:length(startssGEFBA)
    plot(ones(2,1).*startssGEFBA(i),[0,1],'k')
    hold on;
end
for i=1:length(finishessGEFBA)
    plot(ones(2,1).*finishessGEFBA(i),[0,1],'k')
    hold on;
end
ylabel('speech');
ylabel('residual');

subplot(3,1,3);
plot(EstGFD(xaxis_start:xaxis_end)./max(abs(EstGFD(xaxis_start:xaxis_end))),'k');
hold on;
plot(EstimatedGOIs(EstimatedGOIs>=xaxis_start & EstimatedGOIs<=xaxis_end)-(xaxis_start-1),0,'c--o','MarkerSize',7);
hold on;
plot(EstimatedGCIs(EstimatedGCIs>=xaxis_start & EstimatedGCIs<=xaxis_end)-(xaxis_start-1),0,'c--+','MarkerSize',7);
hold on;
for i=1:length(startssGEFBA)
    plot(ones(2,1).*startssGEFBA(i),[0,1],'k')
    hold on;
end
for i=1:length(finishessGEFBA)
    plot(ones(2,1).*finishessGEFBA(i),[0,1],'k')
    hold on;
end
ylabel('GFD');