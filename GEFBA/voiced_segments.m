function [starts,finishes] = voiced_segments(GCIs, GOIs,fs,MinP)
% This function returns the boundaries of the voiced segments. We made two
% main assumptions:
%  1. A voiced segment starts and ends with a GCI.
%  2. Two voiced segments can have a minimum distance of two Max{PP}. If
%  the distance is shorter then we consider them as one unified voiced
%  segment.
%
% Reference: 
%       1. A.I. Koutrouvelis, G.P. Kafentzis, N.D. Gaubitch, R. Heusdens, 
%          "A Fast Method for High-Resolution Voiced/Unvoiced Detection
%          and Glottal Closure/Opening Instant Estimation of Speech"
%
%
% Last modified: 11/07/2015
%
% Author: Andreas Koutrouvelis
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Copyright 2015: Delft University of Technology. The software is free for
% non-commercial use. This program comes WITHOUT ANY WARRANTY.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


GCIs = GCIs(:);
GOIs = GOIs(:);

if(isempty(GCIs) || isempty(GOIs))
    starts = [];
    finishes = [];
    return;
end

MaxPP = round(fs/MinP);


% right boundaries of voiced segments
finishes_pos = find(diff(GCIs)>2*MaxPP);
finishes = GCIs(finishes_pos);
finishes = [finishes; GCIs(end)];
finishes = finishes(:);

% left barriers of voiced segments
starts = [];
starts = [starts; GCIs(1)];
starts = [starts; GCIs(finishes_pos+1)];
% Npos = length(starts);
% for i=1:Npos
%      temp = GOIs(GOIs<starts(i));
%      if(~isempty(temp) && (starts(i)-temp(end))<MaxPP)
%         starts(i) = temp(end);
%      end
% end
starts = starts(:);