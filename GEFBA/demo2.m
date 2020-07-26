% demo2: A complete demonstration and evaluation of GEFBBA.
%        Needs two VOICEBOX functions (i.e., readsfs, sigma algorithm) + 
%        APLAWD database
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

% Speech from APWALD
[sig,fs] = readsfs('as01b0.sfs',1);
[lx,fs] = readsfs('as01b0.sfs',2); % lx is EGG
sig = sig - mean(sig);
start_sig = sig;



WGN_select = 0; % 1: with WGN, 0 without WGN.
if(WGN_select)
    SNR = 20;   % Specify the noise SNR.
    [nnnn,sigmaaaa] = GaussianWhiteNoise(sig,SNR);
    sig = sig + nnnn;
end

    
% propagation time from glottis to lips
propagation_time = 19;   % in APLAWD is 19, in SAM is 14.

lx = [zeros(propagation_time-1,1); lx];
DEGG = filter([1 -1],1,lx);
 

%%%%%%%%%%%%%%%%%%%%%%% Reference Glottal parameters %%%%%%%%%%%%%%%%%%%%
[GCIs,GOIs] = v_sigma(lx,fs); % Extract the reference GCIs/GOIs

% remove short voiced spurts
[GCIs,GOIs] = RemoveShortVoicedGCIsGOIs(GCIs,GOIs,fs);

% find the reference boundaries of VUD
MinP = 80;
[startss,finishess] = voiced_segments(GCIs, GOIs, fs, MinP); 




%%%%%%%%%%%%%%%%%%%%%%%%%%%% Run VGEFBA %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tic;
[EstimatedGOIs, EstimatedGCIs, EstimatedEes,...
 EstimatedFZCIs, EstGFD, EstRes, A] = GEFBA(sig,fs);

% find the boundaries of the voiced segments of GEFBA
[startssGEFBA,finishessGEFBA] =...
    voiced_segments(EstimatedGCIs, EstimatedGOIs,fs,MinP);
GEFBA_speed = toc


%%%%%%%%%%%%%%%%%%%%%%%%%%% EVALUATION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[identsGCI, identsGOI, missesGCI, missesGOI, fasGCI, fasGOI, errorsGCI, ...
 errorsGOI, VUDE_samples] = Eval_GCI_GOI(GOIs,GCIs,EstimatedGOIs,...
                                         EstimatedGCIs,MinP,fs);

Total_NO_LarCyclesGCI = length(GCIs);
Total_NO_LarCyclesGOI = length(GOIs);

GEFBA_IDR_GCI = 100*(identsGCI./Total_NO_LarCyclesGCI)
GEFBA_IDR_GOI = 100*(identsGOI./Total_NO_LarCyclesGOI)
GEFBA_MR_GCI = 100*(missesGCI./Total_NO_LarCyclesGCI)
GEFBA_MR_GOI = 100*(missesGOI./Total_NO_LarCyclesGOI)
GEFBA_FAR_GCI = 100*(fasGCI./Total_NO_LarCyclesGCI)
GEFBA_FAR_GOI = 100*(fasGOI./Total_NO_LarCyclesGOI)
GEFBA_Bias_GCI = mean(errorsGCI)
GEFBA_Bias_GOI = mean(errorsGOI)
GEFBA_std_GCI  = std(errorsGCI)
GEFBA_std_GOI  = std(errorsGOI)
GEFBA_MSE_GCI = (1/length(errorsGCI))*sum(errorsGCI.^2)
GEFBA_MSE_GOI = (1/length(errorsGOI))*sum(errorsGOI.^2)
GEFBA_VUDE = 100*(VUDE_samples./length(sig))






%%%%%%%%%%%%%%%%%%%%%%%%% PLOTTING %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% subplot(4,1,1) contains the following information:
%        a) the speech signal.
%        b) the reference voiced segments boundaries with red lines.
%        c) the reference GCIs and GOIs with red 'x' and red star,
%           respectively.
%
% subplot(4,1,2) contains the following information:
%        a) the dEGG signal.
%        b) the reference voiced segments boundaries with red lines.
%        c) the reference GCIs and GOIs with red 'x' and red star.
%
% subplot(4,1,3) contains the following information:
%        a) the residual.
%        b) the reference voiced segments boundaries with red lines.
%        c) the reference GCIs and GOIs with red 'x' and red star.
%
% subplot(4,1,4) contains the following information:
%        a) the GFD signal.
%        b) the reference and estimated voiced segments boundaries with red
%           lines and black lines, respectively.
%        c) the reference GCIs and GOIs with red 'x' and red star.
%        d) the estimated GCIs and GOIs with black '+' and red 'o'.

figure(1);
skip_Samples = 500;
xaxis_start = 1;
xaxis_end = length(sig)-skip_Samples;
xaxis_ms = xaxis_start:xaxis_end;


subplot(4,1,1);
plot(sig(xaxis_start:xaxis_end)./max(abs(sig(xaxis_start:xaxis_end))),'k');
hold on;
plot(GOIs(GOIs>=xaxis_start & GOIs<=xaxis_end)-(xaxis_start-1),0,'r--p','MarkerSize',7);
hold on;
plot(GCIs(GCIs>=xaxis_start & GCIs<=xaxis_end)-(xaxis_start-1),0,'r--x','MarkerSize',7);
hold on;
for i=1:length(startss)
    plot(ones(2,1).*startss(i),[0,1],'r')
    hold on;
end
for i=1:length(finishess)
    plot(ones(2,1).*finishess(i),[0,1],'r')
    hold on;
end
ylabel('speech');

subplot(4,1,2);
plot(DEGG(xaxis_start:xaxis_end)./max(abs(DEGG(xaxis_start:xaxis_end))),'k');
hold on;
for i=1:length(startss)
    plot(ones(2,1).*startss(i),[0,1],'r')
    hold on;
end
for i=1:length(finishess)
    plot(ones(2,1).*finishess(i),[0,1],'r')
    hold on;
end
hold on;
plot(GOIs(GOIs>=xaxis_start & GOIs<=xaxis_end)-(xaxis_start-1),0,'r--p','MarkerSize',7);
hold on;
plot(GCIs(GCIs>=xaxis_start & GCIs<=xaxis_end)-(xaxis_start-1),0,'r--x','MarkerSize',7);
ylabel('dEGG');

subplot(4,1,3);
plot(EstRes(xaxis_start:xaxis_end)./max(abs(EstRes(xaxis_start:xaxis_end))),'k');
hold on;
plot(GOIs(GOIs>=xaxis_start & GOIs<=xaxis_end)-(xaxis_start-1),0,'r--p','MarkerSize',7);
hold on;
plot(GCIs(GCIs>=xaxis_start & GCIs<=xaxis_end)-(xaxis_start-1),0,'r--x','MarkerSize',7);
hold on;
for i=1:length(startss)
    plot(ones(2,1).*startss(i),[0,1],'r')
    hold on;
end
for i=1:length(finishess)
    plot(ones(2,1).*finishess(i),[0,1],'r')
    hold on;
end
ylabel('speech');
ylabel('residual');

subplot(4,1,4);
plot(EstGFD(xaxis_start:xaxis_end)./max(abs(EstGFD(xaxis_start:xaxis_end))),'k');
hold on;
plot(GOIs(GOIs>=xaxis_start & GOIs<=xaxis_end)-(xaxis_start-1),0,'r--p','MarkerSize',7);
hold on;
plot(GCIs(GCIs>=xaxis_start & GCIs<=xaxis_end)-(xaxis_start-1),0,'r--x','MarkerSize',7);
hold on;
plot(EstimatedGOIs(EstimatedGOIs>=xaxis_start & EstimatedGOIs<=xaxis_end)-(xaxis_start-1),0,'c--o','MarkerSize',7);
hold on;
plot(EstimatedGCIs(EstimatedGCIs>=xaxis_start & EstimatedGCIs<=xaxis_end)-(xaxis_start-1),0,'c--+','MarkerSize',7);
hold on;
for i=1:length(startss)
    plot(ones(2,1).*startss(i),[0,1],'r')
    hold on;
end
for i=1:length(finishess)
    plot(ones(2,1).*finishess(i),[0,1],'r')
    hold on;
end
for i=1:length(startssGEFBA)
    plot(ones(2,1).*startssGEFBA(i),[0,1],'k')
    hold on;
end
for i=1:length(finishessGEFBA)
    plot(ones(2,1).*finishessGEFBA(i),[0,1],'k')
    hold on;
end
ylabel('GFD');