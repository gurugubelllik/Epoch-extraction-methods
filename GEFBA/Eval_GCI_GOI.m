function [identsGCI, identsGOI, missesGCI, missesGOI, fasGCI, fasGOI, errorsGCI, errorsGOI, VUDE_samples] = Eval_GCI_GOI(refGOI,refGCI,estGOI,estGCI,MinP,fs)
% Evaluation.
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


identsGCI = 0;
identsGOI = 0;
missesGCI = 0;
missesGOI = 0;
fasGCI = 0;
fasGOI = 0;
errorsGCI = [];
errorsGOI = [];

if(~isempty(estGCI))
    % first cycle
    if(refGCI(1) - (refGCI(2)-refGCI(1))/2 > 1 )
        x = estGCI((estGCI >= (refGCI(1) - (refGCI(2)-refGCI(1))/2 )) & (estGCI < 0.5*(refGCI(1)+refGCI(2))));
        y = estGCI(estGCI < refGCI(1) - (refGCI(2)-refGCI(1))/2);
    else
        x = estGCI((estGCI >= 0.5*(1 + refGCI(1))) & (estGCI < 0.5*(refGCI(1)+refGCI(2))));
    end
    if(~isempty(x))
        if(length(x)==1)
            identsGCI = identsGCI + 1;
            errorsGCI = [errorsGCI; ((refGCI(1)-x)*1000)/fs];
        else
            fasGCI = fasGCI + 1;
        end
    else
        missesGCI = missesGCI + 1;
    end

    % middle cycles
    N = length(refGCI);
    for i=2:N-1
        if((refGCI(i) - refGCI(i-1) < 2*(fs/MinP)) && (refGCI(i+1) - refGCI(i) < 2*(fs/MinP)))
            x = estGCI((estGCI >= 0.5*(refGCI(i-1)+refGCI(i))) & (estGCI < 0.5*(refGCI(i)+refGCI(i+1))));
        elseif((refGCI(i+1) - refGCI(i) >= 2*(fs/MinP)) && (refGCI(i) - refGCI(i-1) >= 2*(fs/MinP)))
            disp('Not possible!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!')
        elseif(refGCI(i) - refGCI(i-1) >= 2*(fs/MinP))
            x = estGCI((estGCI >= (refGCI(i) - (refGCI(i+1)-refGCI(i))/2)) & (estGCI < 0.5*(refGCI(i)+refGCI(i+1))));
        elseif(refGCI(i+1) - refGCI(i) >= 2*(fs/MinP))
            x = estGCI((estGCI >= 0.5*(refGCI(i-1)+refGCI(i))) & (estGCI < refGCI(i)+ (refGCI(i)-refGCI(i-1))/2));
            if(i+2<=length(refGCI))
                y = estGCI((estGCI >= refGCI(i)+ (refGCI(i)-refGCI(i-1))/2) & (estGCI < refGCI(i+1) - (refGCI(i+2)-refGCI(i+1))/2));
            else
                y = estGCI((estGCI >= refGCI(i)+ (refGCI(i)-refGCI(i-1))/2) & (estGCI < refGCI(i+1) - (refGCI(i)-refGCI(i-1))/2));
            end
        end
        if(~isempty(x))
            if(length(x)==1)
                identsGCI = identsGCI + 1;
                errorsGCI = [errorsGCI; ((refGCI(i)-x)*1000)/fs];
            else
                fasGCI = fasGCI + 1;
            end
        else
            missesGCI = missesGCI + 1;
        end 
    end

    % last cycle
    x = estGCI((estGCI >= 0.5*(refGCI(end-1)+refGCI(end))) & (estGCI < (refGCI(end)+ (refGCI(end) - refGCI(end-1))/2)));
    if(~isempty(x))
        if(length(x)==1)
            identsGCI = identsGCI + 1;
            errorsGCI = [errorsGCI; ((refGCI(end)-x)*1000)/fs];
        else
            fasGCI = fasGCI + 1;
        end
    else
        missesGCI = missesGCI + 1;
    end

    y = estGCI(estGCI > (refGCI(end)+ (refGCI(end) - refGCI(end-1))/2));
else
    missesGCI = length(refGCI);
end


if(~isempty(estGOI))
    % first cycle
    if(refGOI(1) - (refGOI(2)-refGOI(1))/2 > 1 )
        x = estGOI((estGOI >= (refGOI(1) - (refGOI(2)-refGOI(1))/2 )) & (estGOI < 0.5*(refGOI(1)+refGOI(2))));
        y = estGOI(estGOI < refGOI(1) - (refGOI(2)-refGOI(1))/2);
    else
        x = estGOI((estGOI >= 0.5*(1 + refGOI(1))) & (estGOI < 0.5*(refGOI(1)+refGOI(2))));
    end
    if(~isempty(x))
        if(length(x)==1)
            identsGOI = identsGOI + 1;
            errorsGOI = [errorsGOI; ((refGOI(1)-x)*1000)/fs];
        else
            fasGOI = fasGOI + 1;
        end
    else
        missesGOI = missesGOI + 1;
    end



    % middle cycles
    N = length(refGOI);
    for i=2:N-1
        if((refGOI(i) - refGOI(i-1) < 2*(fs/MinP)) && (refGOI(i+1) - refGOI(i) < 2*(fs/MinP)))
            x = estGOI((estGOI >= 0.5*(refGOI(i-1)+refGOI(i))) & (estGOI < 0.5*(refGOI(i)+refGOI(i+1))));
        elseif(refGOI(i) - refGOI(i-1) >= 2*(fs/MinP))
            x = estGOI((estGOI >= (refGOI(i) - (refGOI(i+1)-refGOI(i))/2)) & (estGOI < 0.5*(refGOI(i)+refGOI(i+1))));
        elseif(refGOI(i+1) - refGOI(i) >= 2*(fs/MinP))
            x = estGOI((estGOI >= 0.5*(refGOI(i-1)+refGOI(i))) & (estGOI < refGOI(i)+ (refGOI(i)-refGOI(i-1))/2));
            if(i+2<=length(refGOI))
                y = estGOI((estGOI >= refGOI(i)+ (refGOI(i)-refGOI(i-1))/2) & (estGOI < refGOI(i+1) - (refGOI(i+2)-refGOI(i+1))/2));
            else
                y = estGOI((estGOI >= refGOI(i)+ (refGOI(i)-refGOI(i-1))/2) & (estGOI < refGOI(i+1) - (refGOI(i)-refGOI(i-1))/2));
            end
        end
        if(~isempty(x))
            if(length(x)==1)
                identsGOI = identsGOI + 1;
                errorsGOI = [errorsGOI; ((refGOI(i)-x)*1000)/fs];
            else
                fasGOI = fasGOI + 1;
            end
        else
            missesGOI = missesGOI + 1;
        end 
    end



    % last cycle
    x = estGOI((estGOI >= 0.5*(refGOI(end-1)+refGOI(end))) & (estGOI < (refGOI(end)+ (refGOI(end) - refGOI(end-1))/2)));
    if(~isempty(x))
        if(length(x)==1)
            identsGOI = identsGOI + 1;
            errorsGOI = [errorsGOI; ((refGOI(end)-x)*1000)/fs];
        else
            fasGOI = fasGOI + 1;
        end
    else
        missesGOI = missesGOI + 1;
    end
    y = estGOI(estGOI > (refGOI(end)+ (refGOI(end) - refGOI(end-1))/2));
else
    missesGOI = length(refGOI);
end


[Sr,Fr] = voiced_segments(refGCI, refGOI,fs, MinP);

[Sg, Fg] = voiced_segments(estGCI, estGOI, fs, MinP);

VUDE_samples = VoicedUnvoicedError(Sg,Fg,Sr,Fr);