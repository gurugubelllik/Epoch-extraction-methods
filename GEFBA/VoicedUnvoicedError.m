function VUDE_samples = VoicedUnvoicedError(Sg,Fg,Sr,Fr)
% purpose: 
%         Computes the number of samples that are wrong (in VUD sense)
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

Sg = Sg(:);
Fg = Fg(:);
Sr = Sr(:);
Fr = Fr(:);

N = length(Sg); % Number of estimated voiced segments
M = length(Sr); % Number of reference voiced segments

voicedSpeechGEFBA = [];
for i = 1:N
    voicedSpeechGEFBA = [voicedSpeechGEFBA; [Sg(i):Fg(i)]'];
end


voicedSpeechReference = [];
for i = 1:M
    voicedSpeechReference = [voicedSpeechReference; [Sr(i):Fr(i)]'];
end

VUDE_samples = length(setxor(voicedSpeechGEFBA,voicedSpeechReference));