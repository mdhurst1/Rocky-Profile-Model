%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Code developed by Hironori Matsumoto
% Last update : 28 February 2017
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%% CREATE RECTAUGULAR BREAKING WAVE PRESSURE
%%% INPUT : range ... vertical range of wave pressure
%%% INPUT : shift ... vertical shift of wave pressure
%%% INPUT : FH ...vertical range of breaking wave pressure, true...full range, false...half range  
%%% INPUT : UD ...upper or lower part? when fh_brea is false, 1...Upper part, 2...Lower part
%%% OUTPUT : ssrec_p_brea ... pressure shape function
%%% OUTPUT : ssrec_p_brea_sub ... submarine pressure shape function
%%% OUTPUT : ssrec_brea_dist1 ... Upper limit of pressure shape function
%%% OUTPUT : ssrec_brea_dist2 ... Lower limit of pressure shape function

function [ssrec_p_brea, ssrec_p_brea_sub, ssrec_brea_dist1, ssrec_brea_dist2] =  ssrec_break_p_distri(range, shift, FH, UD)

    %%% Clear ...
    clear ssrec_p_brea ssrec_p_brea_sub ssrec_brea_dist1 ssrec_brea_dist2

    range=ceil(range);
    ssrec_p_brea = zeros(size(range,2),max(range)+1);
    
    %%% SET UPPER AND LOWER LIMIT
    if (FH)
        ssrec_brea_dist1 = ceil((range/2)*(1+shift));
        ssrec_brea_dist2 = floor((range/2)*(1-shift));
    else
        if (UD==1)
            ssrec_brea_dist1 = ceil((range/2)*(1+shift))
            ssrec_brea_dist2 = -(ssrec_brea_dist1 - ceil(range/2));
        else
            ssrec_brea_dist2 = floor((range/2)*(1-shift));
            ssrec_brea_dist1 = -ssrec_brea_dist2 + floor(range/2);
        end
    end
        
    %%% CONSTRUCT SHAPE FUNCTION
    for h=1:size(range,2)
        
        q_break = 1 / abs(ssrec_brea_dist1(h)+ssrec_brea_dist2(h));
        for k=1:(ssrec_brea_dist1(h)+ssrec_brea_dist2(h))
            ssrec_p_brea(h,k) = q_break;
        end        
        ssrec_p_brea_sub(h)= q_break;
        
    end
end