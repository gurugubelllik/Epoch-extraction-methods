function [newGCIs,newGOIs] = RemoveShortVoicedGCIsGOIs(GCIs,GOIs,fs)
% This function removes short voiced spurts. Short voiced spurts are
% unified voiced segments with very short duration (i.e., <= 4 pitch 
% periods).
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

min_pitch  = 80;
min_PP = fs/min_pitch;
% Estimate the initiations of voiced segments
voiced_sepGCI = find(diff(GCIs) > 2*min_PP);
voiced_sepGCI = [1; voiced_sepGCI(:)];

voiced_sepGOI = find(diff(GOIs) > 2*min_PP);
voiced_sepGOI = [1; voiced_sepGOI(:)];

number_of_samples_voiced = 4;
    
tempPos = [];    
for counter_Voiced = 1:length(voiced_sepGCI)-1
    segmentGCI = voiced_sepGCI(counter_Voiced):voiced_sepGCI(counter_Voiced+1);
    if(length(segmentGCI)<=number_of_samples_voiced && segmentGCI(end)<length(GCIs))
        tempPos = [tempPos; segmentGCI(2:end)'];
    end
end
GCIs(tempPos) = [];
newGCIs = GCIs;



tempPos = [];    
for counter_Voiced = 1:length(voiced_sepGOI)-1
    segmentGOI = voiced_sepGOI(counter_Voiced):voiced_sepGOI(counter_Voiced+1);
    if(length(segmentGOI)<=number_of_samples_voiced && segmentGOI(end)<length(GOIs))
        tempPos = [tempPos; segmentGOI(2:end)'];
    end
end
GOIs(tempPos) = [];
newGOIs = GOIs;

