%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Code developed by Hironori Matsumoto
% Last update : 28 February 2017
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%% CREATE RECTAUGULAR BREAKING WAVE PRESSURE
%%% INPUT : range ... vertical range of wave pressure
%%% INPUT : FH ...vertical range of breaking wave pressure, true...full range, false...half range  
%%% INPUT : UD ...upper or lower part? when fh_brea is false, 1...Upper part, 2...Lower part
%%% OUTPUT : ssrec_p_bro ... pressure shape function
%%% OUTPUT : ssrec_p_bro_sub ... submarine pressure shape function
%%% OUTPUT : ssrec_bro_dist1 ... Upper limit of pressure shape function
%%% OUTPUT : ssrec_bro_dist2 ... Lower limit of pressure shape function

function [ssrec_p_bro, ssrec_p_bro_sub, ssrec_bro_dist1, ssrec_bro_dist2] = ssrec_brok_p_distri(range, FH, UD)

    %%% Clear ...
    clear ssrec_p_bro ssrec_p_bro_sub ssrec_bro_dist1 ssrec_bro_dist2
    range=ceil(range);
    ssrec_p_bro = zeros(size(range,2),max(range)+1);
    
    %%% SET UPPER AND LOWER LIMIT
    if (FH)
        ssrec_bro_dist1 = ceil(range/2);
        ssrec_bro_dist2 = floor(range/2);
    else
        if(UD==1)
            ssrec_bro_dist1 = ceil(range/2);
            ssrec_bro_dist2 = [0 0 0];
        else
            ssrec_bro_dist1 = [0 0 0];
            ssrec_bro_dist2 = floor(range/2);
        end
    end

    %%% CONSTRUCT SHAPE FUNCTION
    for i=1:size(range,2)		
        q_bro	= 1/(ssrec_bro_dist1(i)+ssrec_bro_dist2(i));
		for k=1:(ssrec_bro_dist1(i)+ssrec_bro_dist2(i))
			ssrec_p_bro(i,k) = q_bro;
        end
		ssrec_p_bro_sub(i) = q_bro;
    end
    
end