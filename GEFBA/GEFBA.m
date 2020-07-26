function [GOIs, GCIs, Ees, FZCIs, EstGFD, EstRes, A] = GEFBA(signal,fs,params)
% Purpose: 
%       The Glottal closure/opening instants estimation using a 
%       forward-backward algorithm (GEFBA) method estimates the glottal 
%       parameters only in voiced segments.
%
% Arguments:
%       1) signal [Lsig x 1]: The speech signal to be analyzed.
%       2) fs [1 x 1]: The sampling frequency.
%       3) params [struct] (optional): It is a struct array of control 
%          variable vectors of the 3 states: "Strict", "Moderate" and 
%          "Relaxed". Each state gives a different set of values to the 
%          parameters alpha1, alpha2,.... The default values are according 
%          to the published paper [1].
%
%
% Return:
%       1) GOIs: Glottal opening instants.
%       2) GCIs: Glottal closure instants.
%       3) Ees: Negative amplitudes at GCIs.
%       4) FZCIs: first zero crossing instants.
%       5) EstGFD [Shift*Nfr x 1]: It is the estimated glottal flow
%                                  derivative of phase 1. Note that Shift 
%                                  and Nfr are defined inside the code of 
%                                  all_pole_speech_analysis function. 
%                                  Approximetely Shift*Nfr = Lsig.
%       6) EstRes [Shift*Nfr x 1]: It is the linear prediction residual.
%       7) A [p+1 x Nfr]:It is the matrix containing the lpc vectors of
%                        all frames.
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
if(nargin<3)
    % initialization of control variable vectors of the three states: 
    % "Strict", "Moderate" and "Relaxed".
    params = struct('alpha1',   {0.4  0.4  0.65}, ...
                    'alpha2',   {1.6  1.6  1.4 }, ...
                    'beta1',    {0.9  0.85 0.75}, ...
                    'beta2',    {1.1  1.15 1.3 }, ...
                    'gama1',    {0.4  0.4  0.3 }, ...
                    'gama2',    {2    2.5  3.5 }, ...
                    'delta1',   {0    0.6  0.55}, ...
                    'delta2',   {0    1.5  1.6   }, ...
                    'epsilon1', {0.5  0.3  0.25 }, ...
                    'epsilon2', {1    2.6  2.7   }, ...
                    'zeta1',    {0.3  0.35  0.4}, ...
                    'zeta2',    {3    3.1  3.5   });
end

signal = signal(:);


N = 2*round(fs/60);     % Length of analysis frame
p = round(fs/1000)+16;  % Order of LPC analysis.

%%%%%%%%%%%%%%%%%%%%%%%%%%%% Phase 1 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[EstRes,EstGFD,A] = all_pole_speech_analysis(signal,p,N);


%%%%%%%%%%%%%%%%%%%%%%%%%%%% Phase 2 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[GCIs,FZCIs,GOIs] = EstimateGlottalParams(EstGFD,fs,params);
GCIs = unique(GCIs);
FZCIs = unique(FZCIs);
Ees = EstGFD(GCIs);


% The following two commands fix a bag inside the code
GOIs(diff(GOIs)<fs/800) = [];  % At the boundaries of voiced unvoiced detections
                               % sometimes there are two very close gois 
                               % that only one should be kept

FZCIs(diff(FZCIs)<fs/800) = []; % At the boundaries of voiced unvoiced detections
                                % sometimes there are two very close gois 
                                % that only one should be kept                               
end





function [GCIs,FZCIs,GOIs] = EstimateGlottalParams(GFD_estimate,fs,params)
% Purpose: 
%          This function is phase 2 which consists of two main steps: 
%          step1 and step2. Step1 finds a highly voiced frame and its 
%          glottal params. Step2 fills the voiced gaps left and right of 
%          the highly voiced frame.
%
% Arguments:
%          1) GFD_estimate: the glottal flow derivative estimate.
%          2) fs: the sampling frequency.
%          3) params: struct array of the three control variable vectors.
%
% Return:
%          1) GCIs: estimated glottal closure instants.
%          2) FZCIs: estimated first zero crossing instants.
%          3) GOIs: estimated glottal opening instants.

LGFD = length(GFD_estimate);
TemporaryGFD = GFD_estimate;
TemporaryGFD(TemporaryGFD>=0) = 0;
inv_negative_part = -TemporaryGFD;
GCIs = [];
FZCIs = [];
GOIs = [];
pitchesss = [];
frame_length = round(fs/20);  % 4 times the maximum pitch period
Shift = round(frame_length/2);
interv = 1:frame_length;

VU = 0; % 0: not highly voiced, 1: highly voiced, 2: voiced, 3: unvoiced
pitch_period = 0;

%%%%%%%%%%%%%%%%%%%%%%%%%%%% FSM starts %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
while(1)

    if(VU == 1)
        % Step 2: Filling the voiced gaps
        temporary1 = GCIs(GCIs<TempLeftTes(1));
        if(~isempty(temporary1))
            interv_stop = temporary1(end)+round(fs/200);
        else
            interv_stop = 1;
        end
        startPos = TempLeftTes(1);
        [TempLeftTes,TempLeftZeroCross,TempLeftGOIs,TempPitchesLeft] = MB(inv_negative_part,GFD_estimate,startPos,VU,pitch_period,interv_stop,fs,params);
        TempLeftTes(1) = [];
        TempLeftGOIs(1) = [];
        TempLeftTes = sort(TempLeftTes,'ascend');
        TempLeftZeroCross = sort(TempLeftZeroCross,'ascend');
        TempLeftGOIs = sort(TempLeftGOIs,'ascend');
         

        
        startPos = TempRightTes(end);
        interv_stop = length(inv_negative_part);
        [TempRightTes,TempRightZeroCross,TempRightGOIs,TempPitchesRight] = MF(inv_negative_part,GFD_estimate,startPos,VU,pitch_period,interv_stop,fs,params);
        
     
        VU = 3; % unvoiced
        if(~isempty(TempRightTes))
            TempRightTes(1) = [];
            TempRightZeroCross(1) = [];
            TempRightGOIs(1) = [];
        end
        TempTes = [TempLeftTes(:); TempRightTes(:)];    
        FZCIs = [FZCIs;TempLeftZeroCross(:);TempRightZeroCross(:)];
        GOIs = [GOIs;TempLeftGOIs(:);TempRightGOIs(:)];
        pitchesss = [pitchesss; TempPitchesLeft;TempPitchesRight];
        

        GCIs = [GCIs; TempTes];
        if(~isempty(pitchesss))
            pitch_period = mean(pitchesss);
        else
            pitch_period = 0;
        end
    elseif(VU == 0)
        % Step 1: Searching for a highly voiced frame.
        
        
        % the minimum negative peak of the frame is found as a starting 
        % reference candidate GCI.
        [peaks,locs] = findpeaks(inv_negative_part(interv),'SORTSTR','descend');
        while(length(locs)<3)
            interv = interv + Shift;
            if(interv(end)<LGFD)
                [peaks,locs] = findpeaks(inv_negative_part(interv),'SORTSTR','descend');
            else
                return;
            end
        end
        
        startPos = locs(1) + interv(1) - 1;
        interv_stop = interv(1);
        [TempLeftTes,TempLeftZeroCross,TempLeftGOIs,TempPitchesLeft] = MB(inv_negative_part,GFD_estimate,startPos,VU,pitch_period,interv_stop,fs,params);
        TempLeftTes = sort(TempLeftTes,'ascend');
        TempLeftZeroCross = sort(TempLeftZeroCross,'ascend');
        TempLeftGOIs = sort(TempLeftGOIs,'ascend');

        pitch_period = 0;
        interv_stop = interv(end);
        [TempRightTes,TempRightZeroCross,TempRightGOIs,TempPitchesRight] = MF(inv_negative_part,GFD_estimate,startPos,VU,pitch_period,interv_stop,fs,params);
        if(~isempty(TempRightTes))
            TempValueRight = TempRightTes(1);
            TempZcValueRight = TempRightZeroCross(1);
            TempGOIValueRight = TempRightGOIs(1);
            TempRightTes(1) = [];
            TempRightZeroCross(1) = [];
            TempRightGOIs(1) = [];
        end
        
        TempTes = [TempLeftTes(:); TempRightTes(:)];

        TempZerCrosss = [TempLeftZeroCross(:); TempRightZeroCross(:)];
        TempGOIs = [TempLeftGOIs(:); TempRightGOIs(:)];

        TempPitches = [TempPitchesLeft;TempPitchesRight];

        DifferencesTempTes = diff(TempTes);
        DifferencesTempZerCrosss = diff(TempZerCrosss);

        % Threshold checking
        if(length(TempTes)>=4 && (min(DifferencesTempTes) > 0.6*max(DifferencesTempTes)) && (min(DifferencesTempZerCrosss) > 0.4*max(DifferencesTempZerCrosss)))
            VU = 1;
        else
            VU = 0;
        end
        
        if(VU ==0)
            pitchesss = [];
        end

        if(VU == 1)
            if(isempty(TempRightTes))
                TempRightTes = TempValueRight;
                TempRightZeroCross = TempZcValueRight;
                TempRightGOIs = TempGOIValueRight;
            end
            GCIs = [GCIs; TempTes];
            FZCIs = [FZCIs; TempZerCrosss];
            GOIs = [GOIs; TempGOIs];
            pitchesss = [pitchesss; TempPitches];
            pitch_period = mean(pitchesss);
        else
            pitch_period = 0;
            % do nothing
        end
    end
    
    GCIs = sort(GCIs,'ascend');
    FZCIs = sort(FZCIs,'ascend');
    GOIs = sort(GOIs,'ascend');
    
    % Move to next interval
    if(VU == 3)
        interv = GCIs(end)+round(fs/300):GCIs(end)+round(fs/300)+frame_length;
        VU = 0;
    elseif(VU == 0)
        interv = interv + Shift;
    end
    if(interv(end)>LGFD)
        return;
    end

end

end



function [TempLeftTes,TempLeftZeroCross,TempLeftGOIs,TempPitchesLeft] = MB(inv_negative_part,GFD_estimate,startPos,VU,pitch_period,interv_stop,fs,params)
% Move Backward

% Initialization
maximum_possible_pitch = 500; % child
minimum_possible_period = round(fs/maximum_possible_pitch);
minimum_possible_pitch = 80;  % Bass male
maximum_possible_period = round(fs/minimum_possible_pitch);
search_int = round(fs/100);
current_peak_loc = startPos;
TempLeftTes = [];
TempLeftTes = [TempLeftTes; current_peak_loc];
TempLeftZeroCross = [];
TempLeftGOIs = [];
TempPitchesLeft = [];


current_ZC_loc = FindZeroCrossLeft(inv_negative_part,current_peak_loc,search_int);
[current_goi_loc,current_tm] = FindGOILeftFZCI(GFD_estimate,current_ZC_loc,current_peak_loc,0,search_int);
if(isempty(current_ZC_loc) || isempty(current_goi_loc))
    TempLeftZeroCross = [TempLeftZeroCross; current_peak_loc];
    TempLeftGOIs = [TempLeftGOIs; current_peak_loc];
	return;
end
currentDifferenceZC_Te = current_peak_loc - current_ZC_loc;


TempLeftZeroCross = [TempLeftZeroCross; current_ZC_loc];
TempLeftGOIs = [TempLeftGOIs; current_goi_loc];


PeaksPos = [];
pitch_period_temp = 0;
while(current_peak_loc > interv_stop)
    if(VU == 1)
        % State: "relaxed"
        left = params(3).alpha1;
        right = params(3).alpha2;
        leftZC = params(3).beta1;
        rightZC = params(3).beta2;
        leftGOI = params(3).delta1;
        rightGOI = params(3).delta2;
        leftZC_diff = params(3).gama1;
        rightZC_diff = params(3).gama2;
        if(pitch_period_temp)
            Left_Barier = round(left*pitch_period_temp);
            Right_Barier = round(right*pitch_period_temp);
        else
            Left_Barier = round(left*pitch_period);
            Right_Barier = round(right*pitch_period);
        end
        down = params(3).epsilon1;
        up = params(3).epsilon2;
        down_tm = params(3).zeta1;
        up_tm = params(3).zeta2;
    elseif(VU && ~pitch_period)
        error('It should not be here');
    elseif(VU==0 && ~pitch_period)
        % State: "strict"
        leftZC_diff = params(1).gama1;
        rightZC_diff = params(1).gama2;
        Right_Barier = maximum_possible_period;
        Left_Barier = minimum_possible_period;
        down = params(1).epsilon1;
        up = params(1).epsilon2;
        down_tm = params(1).zeta1;
        up_tm = params(1).zeta2;
    elseif(VU==0 && pitch_period)
        % State: "moderate"
        left = params(2).alpha1;
        right = params(2).alpha2;
        leftZC = params(2).beta1;
        rightZC = params(2).beta2;
        leftGOI = params(2).delta1;
        rightGOI = params(2).delta2;
        leftZC_diff = params(2).gama1;
        rightZC_diff = params(2).gama2;
        Left_Barier = round(left*pitch_period);
        Right_Barier = round(right*pitch_period);
        down = params(2).epsilon1;
        up = params(2).epsilon2;
        down_tm = params(1).zeta1;
        up_tm = params(1).zeta2;
    else
        error('It should not be here!!!!')
    end
    
    
    
    if(current_peak_loc - Right_Barier > interv_stop)
        current_area = inv_negative_part(current_peak_loc-Right_Barier : current_peak_loc-Left_Barier);
        offsseet = current_peak_loc-Right_Barier;
    elseif((current_peak_loc - Right_Barier <= interv_stop) && (current_peak_loc-Left_Barier > interv_stop))
        current_area = inv_negative_part(interv_stop : current_peak_loc-Left_Barier);
        offsseet = interv_stop;
    else
        return;
    end
    
    if(length(current_area)<3)
        return;
    end
    
    if(pitch_period)
        [SortedPeaks,SortedPeaksPos] = findpeaks(current_area,'SORTSTR','descend');
        SortedPeaksPos = SortedPeaksPos + offsseet-1;
        if(~isempty(SortedPeaksPos))
                [SortedPeaksPos] = FindCandidateGCIs(SortedPeaksPos,inv_negative_part);
                TeNotFound = 1;
                while(~isempty(SortedPeaksPos) && TeNotFound)
                    % check one by one all candidate glottal parameter set
                    % till you find a valid (i.e., conditions to be satisfied)
                    if(VU)
                        [max_value_intersec, Postemp_peak_intersection] = max(inv_negative_part(SortedPeaksPos));
                        temp_current_peak_loc = SortedPeaksPos(Postemp_peak_intersection);
                    else
                        [minimum_difference,Postemp_peak_intersection] = min( current_peak_loc-SortedPeaksPos );
                        temp_current_peak_loc = SortedPeaksPos(Postemp_peak_intersection);
                    end
                    CorrespZeroCrossSortedPeaksPos = FindZeroCrossLeft(inv_negative_part,temp_current_peak_loc,search_int);
                    if(isempty(CorrespZeroCrossSortedPeaksPos))
                    	return;
                    end                   
                    TempDifferenceZC_Te = temp_current_peak_loc - CorrespZeroCrossSortedPeaksPos;
                    
                 
                    pitch_period_temp = current_peak_loc - temp_current_peak_loc;
                    
                    
                    [Corresp_goi_SortedPeaksPos, corresp_tm] = FindGOILeftFZCI(GFD_estimate,CorrespZeroCrossSortedPeaksPos,temp_current_peak_loc,pitch_period_temp,search_int);
                    if(isempty(Corresp_goi_SortedPeaksPos))
                        return;
                    end
                    
                    FZCI_diff = current_ZC_loc - CorrespZeroCrossSortedPeaksPos;
                    GOI_diff = current_goi_loc - Corresp_goi_SortedPeaksPos;
                    
                    
                    % Conditions
                    if(inv_negative_part(temp_current_peak_loc) >= down*inv_negative_part(current_peak_loc) ...
                       && inv_negative_part(temp_current_peak_loc) <= up*inv_negative_part(current_peak_loc) ...
                       && corresp_tm >= down_tm*current_tm ...
                       && corresp_tm <= up_tm*current_tm ...
                       && pitch_period_temp >= left*pitch_period ...
                       && pitch_period_temp <= right*pitch_period ...
                       && FZCI_diff >= leftZC*pitch_period_temp ...
                       && FZCI_diff <= rightZC*pitch_period_temp ...
                       && TempDifferenceZC_Te>=leftZC_diff*currentDifferenceZC_Te ...
                       && TempDifferenceZC_Te<=rightZC_diff*currentDifferenceZC_Te ...
                       && GOI_diff >= leftGOI*pitch_period_temp ...
                       && GOI_diff <= rightGOI*pitch_period_temp)
                   
                        TempPitchesLeft = [TempPitchesLeft; pitch_period_temp];
                        current_peak_loc = SortedPeaksPos(Postemp_peak_intersection);
                        TempLeftTes = [TempLeftTes; current_peak_loc];
                        current_ZC_loc = CorrespZeroCrossSortedPeaksPos;
                        TempLeftZeroCross = [TempLeftZeroCross; current_ZC_loc];    
                        currentDifferenceZC_Te = TempDifferenceZC_Te;
                        current_goi_loc = Corresp_goi_SortedPeaksPos;
                        TempLeftGOIs = [TempLeftGOIs; current_goi_loc];
                        TeNotFound = 0;
                        pitch_period = pitch_period_temp;
                        current_tm = corresp_tm;
                    else
                        % this candidate glottal parameter set is not valid
                        SortedPeaksPos(Postemp_peak_intersection) = [];
                    end
                end 
                if(TeNotFound)
                    return;
                end
        else
            return;
        end
    else
        [Peaks,PeaksPos] = findpeaks(current_area);
        PeaksPos = PeaksPos + offsseet-1;
        % quantize peaks to local maxima that have zero crossing between
        % them
        PeaksPos = [PeaksPos; current_peak_loc];
        PeaksPos = FindCandidateGCIs(PeaksPos,inv_negative_part);
        PeaksPos(end) = [];
        newPeaksPos = PeaksPos;
        
        % If no pitch period is available then select the closest peak to current_peak_loc with a similar amplitude
        if(~isempty(newPeaksPos))
            TeNotFound = 1;
            while(~isempty(newPeaksPos) && TeNotFound)
            	[minimum_difference,minimum_pos] = min( current_peak_loc - newPeaksPos );
            	temp_current_peak_loc = newPeaksPos(minimum_pos);
                current_ZC_loc = FindZeroCrossLeft(inv_negative_part,temp_current_peak_loc,search_int);
                if(isempty(current_ZC_loc))
                    return;
                end
                TempDifferenceZC_Te = temp_current_peak_loc - current_ZC_loc;
               
                % Conditions
             	if(inv_negative_part(temp_current_peak_loc)>= down*inv_negative_part(current_peak_loc) ...
                   && inv_negative_part(temp_current_peak_loc)<= up*inv_negative_part(current_peak_loc) ...
                   && TempDifferenceZC_Te>=leftZC_diff*currentDifferenceZC_Te ...
                   && TempDifferenceZC_Te<=rightZC_diff*currentDifferenceZC_Te)
                	pitch_period = current_peak_loc - temp_current_peak_loc;
                    [Corresp_goi, corresp_tm] = FindGOILeftFZCI(GFD_estimate,current_ZC_loc,temp_current_peak_loc,pitch_period,search_int);
                    if(isempty(Corresp_goi))
                        return;
                    end
                    if(corresp_tm < down_tm*current_tm || corresp_tm > up_tm*current_tm)
                        newPeaksPos(minimum_pos) = [];
                        continue;
                    end
                    TempLeftTes = [TempLeftTes ; temp_current_peak_loc];
                	TempPitchesLeft = [TempPitchesLeft; pitch_period];
                	current_peak_loc = temp_current_peak_loc;
                    TempLeftZeroCross = [TempLeftZeroCross; current_ZC_loc];
                    TempLeftGOIs = [TempLeftGOIs; Corresp_goi];
                    current_goi_loc = Corresp_goi;
                    currentDifferenceZC_Te = TempDifferenceZC_Te;
                	TeNotFound = 0;
                    current_tm = corresp_tm;
                 else
                	newPeaksPos(minimum_pos) = [];
                 end
            end
            if(TeNotFound)
                return;
            end
        else
           return;
        end 
    end
    
    
end

end





function [TempRightTes,TempRightZeroCross,TempRightGOIs,TempPitchesRight] = MF(inv_negative_part,GFD_estimate,startPos,VU,pitch_period,interv_stop,fs,params)
% Move Forward

% Initialization
maximum_possible_pitch = 500; % child
minimum_possible_period = round(fs/maximum_possible_pitch);
minimum_possible_pitch = 80;  % Bass male
maximum_possible_period = round(fs/minimum_possible_pitch);
search_int = round(fs/100);
current_peak_loc = startPos;
TempRightTes = [];
TempRightTes = [TempRightTes; current_peak_loc];
TempRightZeroCross = [];
TempRightGOIs = [];
TempPitchesRight = [];


current_ZC_loc = FindZeroCrossLeft(inv_negative_part,current_peak_loc,search_int);
[current_goi_loc, current_tm] = FindGOILeftFZCI(GFD_estimate,current_ZC_loc,current_peak_loc,0,search_int);
if(isempty(current_ZC_loc) || isempty(current_goi_loc))
    TempRightZeroCross = [TempRightZeroCross; current_peak_loc];
    TempRightGOIs = [TempRightGOIs; current_peak_loc];
    return;
end
currentDifferenceZC_Te = current_peak_loc - current_ZC_loc;



TempRightZeroCross = [TempRightZeroCross; current_ZC_loc];
TempRightGOIs = [TempRightGOIs; current_goi_loc];




PeaksPos = [];
pitch_period_temp = 0;
while(current_peak_loc < interv_stop)    
    if(VU == 1)
        % State: "relaxed"
        left = params(3).alpha1;
        right = params(3).alpha2;
        leftZC = params(3).beta1;
        rightZC = params(3).beta2;
        leftGOI = params(3).delta1;
        rightGOI = params(3).delta2;
        leftZC_diff = params(3).gama1;
        rightZC_diff = params(3).gama2;
        if(pitch_period_temp)
            Left_Barier = round(left*pitch_period_temp);
            Right_Barier = round(right*pitch_period_temp);
        else
            Left_Barier = round(left*pitch_period);
            Right_Barier = round(right*pitch_period);
        end
        down = params(3).epsilon1;
        up = params(3).epsilon2;
        down_tm = params(3).zeta1;
        up_tm = params(3).zeta2;
    elseif(VU==1 && ~pitch_period)
        error('It should not be here');
    elseif(VU==0 && ~pitch_period)
        % State: "strict"
        leftZC_diff = params(1).gama1;
        rightZC_diff = params(1).gama2;
        Right_Barier = maximum_possible_period;
        Left_Barier = minimum_possible_period;
        down = params(1).epsilon1;
        up = params(1).epsilon2;
        down_tm = params(1).zeta1;
        up_tm = params(1).zeta2;
    elseif(VU==0 && pitch_period)
        % State: "moderate"
        left = params(2).alpha1;
        right = params(2).alpha2;
        leftZC = params(2).beta1;
        rightZC = params(2).beta2;
        leftGOI = params(2).delta1;
        rightGOI = params(2).delta2;
        leftZC_diff = params(2).gama1;
        rightZC_diff = params(2).gama2;
        Left_Barier = round(left*pitch_period);
        Right_Barier = round(right*pitch_period);
        down = params(2).epsilon1;
        up = params(2).epsilon2;
        down_tm = params(1).zeta1;
        up_tm = params(1).zeta2;
    else
        error('It should not be here!!!!')
    end
        
    
    
    if(current_peak_loc + Right_Barier < interv_stop)
        current_area = inv_negative_part(current_peak_loc+Left_Barier : current_peak_loc+Right_Barier);
        offsseet = current_peak_loc+Left_Barier;
    elseif((current_peak_loc + Right_Barier >= interv_stop) && (current_peak_loc+Left_Barier < interv_stop))
        current_area = inv_negative_part(current_peak_loc+Left_Barier:interv_stop);
        offsseet = current_peak_loc+Left_Barier;
    else
        return;
    end
    
    if(length(current_area)<3)
        return;
    end
    
    
    if(pitch_period)
        [SortedPeaks,SortedPeaksPos] = findpeaks(current_area,'SORTSTR','descend');
        SortedPeaksPos = SortedPeaksPos + offsseet-1;        
        if(~isempty(SortedPeaksPos))
                [SortedPeaksPos] = FindCandidateGCIs(SortedPeaksPos,inv_negative_part);
                TeNotFound = 1;
                while(~isempty(SortedPeaksPos) && TeNotFound)
                    % check one by one all candidate glottal parameter set
                    % till you find a valid (i.e., conditions to be satisfied)
                    if(VU)
                        [max_value_intersec, Postemp_peak_intersection] = max(inv_negative_part(SortedPeaksPos));
                        temp_current_peak_loc = SortedPeaksPos(Postemp_peak_intersection);
                    else
                        [minimum_difference,Postemp_peak_intersection] = min( SortedPeaksPos - current_peak_loc );
                        temp_current_peak_loc = SortedPeaksPos(Postemp_peak_intersection);
                    end
                    
                    CorrespZeroCrossSortedPeaksPos = FindZeroCrossLeft(inv_negative_part,temp_current_peak_loc,search_int);
                    if(isempty(CorrespZeroCrossSortedPeaksPos))
                    	return;
                    end
                    TempDifferenceZC_Te = temp_current_peak_loc - CorrespZeroCrossSortedPeaksPos;    
        
                    pitch_period_temp = temp_current_peak_loc-current_peak_loc;
                    
                    [Corresp_goi_SortedPeaksPos, corresp_tm] = FindGOILeftFZCI(GFD_estimate,CorrespZeroCrossSortedPeaksPos,temp_current_peak_loc,pitch_period_temp,search_int);
                    if(isempty(Corresp_goi_SortedPeaksPos))
                        return;
                    end
                    
                    FZCI_diff = CorrespZeroCrossSortedPeaksPos-current_ZC_loc;
                    GOI_diff = Corresp_goi_SortedPeaksPos-current_goi_loc;
                    
                    % Conditions
                    if(inv_negative_part(temp_current_peak_loc) >= down*inv_negative_part(current_peak_loc) ...
                       && inv_negative_part(temp_current_peak_loc) <= up*inv_negative_part(current_peak_loc) ...
                       && corresp_tm >= down_tm*current_tm ...
                       && corresp_tm <= up_tm*current_tm ...
                       && pitch_period_temp >= left*pitch_period ...
                       && pitch_period_temp <= right*pitch_period ...
                       && FZCI_diff >= leftZC*pitch_period_temp ...
                       && FZCI_diff <= rightZC*pitch_period_temp ...
                       && TempDifferenceZC_Te>=leftZC_diff*currentDifferenceZC_Te ...
                       && TempDifferenceZC_Te<=rightZC_diff*currentDifferenceZC_Te ...
                       && GOI_diff >= leftGOI*pitch_period_temp ...
                       && GOI_diff <= rightGOI*pitch_period_temp)
                   
                        TempPitchesRight = [TempPitchesRight; pitch_period_temp];
                        current_peak_loc = SortedPeaksPos(Postemp_peak_intersection);
                        TempRightTes = [TempRightTes; current_peak_loc];
                        current_ZC_loc = CorrespZeroCrossSortedPeaksPos;
                        TempRightZeroCross = [TempRightZeroCross; current_ZC_loc];
                        currentDifferenceZC_Te = TempDifferenceZC_Te;
                        current_goi_loc = Corresp_goi_SortedPeaksPos;
                        TempRightGOIs = [TempRightGOIs; current_goi_loc];
                        TeNotFound = 0;
                        pitch_period = pitch_period_temp;
                        current_tm = corresp_tm;
                    else
                        % this candidate glottal parameter set is not valid
                        SortedPeaksPos(Postemp_peak_intersection) = [];
                    end
                end 
                if(TeNotFound)
                    return;
                end
        else
            return;
        end
    
    else
        [Peaks,PeaksPos] = findpeaks(current_area);
        PeaksPos = PeaksPos + offsseet-1;
        PeaksPos = [current_peak_loc; PeaksPos];
        % quantize peaks to local maxima that have zero crossing between
        % them
        PeaksPos = FindCandidateGCIs(PeaksPos,inv_negative_part);
        PeaksPos(1) = [];
        newPeaksPos = PeaksPos;
        
        
        % If no pitch period is available then select the closest peak to 
        % current_peak_loc if the amplitude is close otherwise select the 
        % next peak with a close amplitude
        if(~isempty(newPeaksPos))
            TeNotFound = 1;
            while(~isempty(newPeaksPos) && TeNotFound)
                [minimum_difference,minimum_pos] = min( newPeaksPos - current_peak_loc );
                temp_current_peak_loc = newPeaksPos(minimum_pos);
                
                current_ZC_loc = FindZeroCrossLeft(inv_negative_part,temp_current_peak_loc,search_int);
                if(isempty(current_ZC_loc))
                    return;
                end
                TempDifferenceZC_Te = temp_current_peak_loc - current_ZC_loc;
                
               
                % Conditions
                if(inv_negative_part(temp_current_peak_loc)>= down*inv_negative_part(current_peak_loc) ...
                   && inv_negative_part(temp_current_peak_loc)<= up*inv_negative_part(current_peak_loc) ...
                   && TempDifferenceZC_Te>=leftZC_diff*currentDifferenceZC_Te ...
                   && TempDifferenceZC_Te<=rightZC_diff*currentDifferenceZC_Te)
                	pitch_period = temp_current_peak_loc-current_peak_loc;
                    [Corresp_goi, corresp_tm] = FindGOILeftFZCI(GFD_estimate,current_ZC_loc,temp_current_peak_loc,pitch_period,search_int);
                    if(isempty(Corresp_goi))
                        return;
                    end
                    if(corresp_tm < down_tm*current_tm || corresp_tm > up_tm*current_tm)
                        newPeaksPos(minimum_pos) = [];
                        continue;
                    end
                    TempRightTes = [TempRightTes ; temp_current_peak_loc];
                	TempPitchesRight = [TempPitchesRight; pitch_period];
                	current_peak_loc = temp_current_peak_loc;
                    TempRightZeroCross = [TempRightZeroCross; current_ZC_loc];
                    TempRightGOIs = [TempRightGOIs; Corresp_goi];
                    current_goi_loc = Corresp_goi;
                    currentDifferenceZC_Te = abs(current_peak_loc - current_ZC_loc);
                    TeNotFound = 0;
                    current_tm = corresp_tm;
                else
                	newPeaksPos(minimum_pos) = [];
                end
            end
            if(TeNotFound)
                return;
            end
        else
           return;
        end 
    end
end
end



function quantized_peaks = FindCandidateGCIs(PeaksPos,inv_negative_part)
% Find the minimum negative peak between each pair of neighbouring zero 
% crossings

i = 1;
while(length(PeaksPos)>=2 && i<length(PeaksPos))
    if(~any( inv_negative_part(PeaksPos(i):PeaksPos(i+1)) == 0))
        if(inv_negative_part(PeaksPos(i+1))>= inv_negative_part(PeaksPos(i)))
        	PeaksPos(i) = [];
        else
            PeaksPos(i+1) = [];
        end
    else
        i = i + 1;
    end
end
quantized_peaks = PeaksPos;

end




function ZeroCrossLeft = FindZeroCrossLeft(inv_negative_part,Te,search_int)
% Find the first zero Crossing left to Te

if(Te-search_int > 1)
    S1 = inv_negative_part(Te-search_int:Te);
else
    S1 = inv_negative_part(1:Te);
    search_int = Te;
end
zeroCrossingsS = find( S1 == 0 );
if(~isempty(zeroCrossingsS))                                   
	ZeroCrossLeft = zeroCrossingsS(end) + Te-search_int-1;
else
	ZeroCrossLeft = [];
end

end






function [goi,max_amp] = FindGOILeftFZCI(GFD_estimateee,tp,te,pitch_period,search_int)
% Find GOI= the first zero Crossing/zero left to tp

% thresh = 2*tp - te; % goi should be less than thresh
max_amp = 0;
thresh = tp - 2;
thresh_begin = round(te-0.8*pitch_period);
Candgois = [];

if(thresh-search_int > 1)
    S1 = GFD_estimateee(thresh-search_int:thresh-1).*GFD_estimateee(thresh-search_int+1:thresh);
else
    S1 = GFD_estimateee(1:thresh-1).*GFD_estimateee(2:thresh);
    search_int = thresh-1;
end

zeroCrossingsS = find( S1 <= 0 );
if(~isempty(zeroCrossingsS))       
    Candgois = zeroCrossingsS(:) + thresh-search_int-1;
    % Discard all positive zero-crossings.
    % Note that only the negative ones are usefull!
    Candgois(GFD_estimateee(Candgois)>0) = [];

       
    if(~isempty(Candgois))
        if(pitch_period~=0)
            if(thresh_begin<tp-1 && thresh_begin>1)
                [max_amp,max_pos] = max(GFD_estimateee(thresh_begin:tp-1));
                max_pos = max_pos+thresh_begin-1;
                Candgois = Candgois(Candgois < max_pos);
                if(~isempty(Candgois))
                   goi = Candgois(end);  % Zero crossing after previous te
                   goi2 = find(GFD_estimateee(goi+1:max_pos-1)<max_amp/2.5);
                   if(~isempty(goi2))
                      goi = goi+goi2(end);
                   end
                else
                   goi = [];
                end
            else
                goi = Candgois(end);
            end
        else
            goi = Candgois(end);
            [max_amp,max_pos] = max(GFD_estimateee(goi:tp-1));
        end
    else
        goi = [];
    end
else
	goi = [];
end

end