%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Code developed by Hironori Matsumoto
% Last update : 28 February 2017
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%% CREATE EXPONENTIAL BREAKING WAVE PRESSURE
%%% INPUT : range ... vertical range of wave pressure
%%% INPUT : shift ... vertical shift of wave pressure
%%% INPUT : FH ...vertical range of breaking wave pressure, true...full range, false...half range  
%%% INPUT : UD ...upper or lower part? when fh_brea is false, 1...Upper part, 2...Lower part
%%% OUTPUT : exp_p_brea ... pressure shape function
%%% OUTPUT : exp_p_brea_sub ... submarine pressure shape function
%%% OUTPUT : exp_brea_dist1 ... Upper limit of pressure shape function
%%% OUTPUT : exp_brea_dist2 ... Lower limit of pressure shape function

function [exp_p_brea, exp_p_brea_sub, exp_brea_dist1, exp_brea_dist2] = ori_break_p_distri(range, shift, FH, UD)

    %%% Clear ...
    clear exp_p_brea exp_p_brea_sub exp_brea_dist1 exp_brea_dist2
    range = ceil(range);
    exp_p_brea = zeros(3,range(3));
    
    %%% SET UPPER AND LOWER LIMIT
    if (FH)
        exp_brea_dist1 = ceil((range/2)*(1+shift));
        exp_brea_dist2 = floor((range/2)*(1-shift));
    else
        if (UD==1)
            exp_brea_dist1 = ceil((range/2)*(1+shift));
            exp_brea_dist2 =  -(exp_brea_dist1 - ceil(range/2));
        else
            exp_brea_dist2 = floor((range/2)*(1-shift));
            exp_brea_dist1 = -ssrec_brea_dist2 + floor(range/2);
        end
    end
    
    %%% CONSTRUCT SHAPE FUNCTION
    for h=1:length(range)
        total = 0;
        if ( mod(range(h),2)==0 )
            i = range(h)/2+1;
            j = range(h)/2;
            for k=range(h)/2:-1:2
                i=i-1;
                j=j+1;
                exp_p_brea(h,i) = exp(-abs(range(h)/2-k)/(0.1*(range(h)/2)));
                exp_p_brea(h,j) = exp(-abs(range(h)/2-k)/(0.1*(range(h)/2)));
            end
            j=j+1;
            exp_p_brea(h,1) = 0;
            exp_p_brea(h,j) = 0;
            
            if (FH)       
                total = sum(exp_p_brea,2);
                exp_p_brea(h,:) = exp_p_brea(h,:) / total(h);
                if(floor((1-shift)*(range(h)/2))==0)    exp_p_brea_sub(h) = exp_p_brea(h,1);
                else                                    exp_p_brea_sub(h) = exp_p_brea(h,floor((1-shift)*(range(h)/2)));
                end
            else
                if (UD==1) 
                    exp_p_brea(h,range(h)/2+1:range(h)) = 0;
                    total = sum(exp_p_brea,2);
                    exp_p_brea(h,:) = exp_p_brea(h,:) / total(h);
                    exp_p_brea_sub(h) = exp_p_brea(h,floor((1-shift)*(range(h)/2)));
                else
                    exp_p_brea(h,1:range(h)/2) = exp_p_brea(h,range(h)/2+1:range(h));
                    exp_p_brea(h,range(h)/2+1:range(h)) = 0;
                    total = sum(exp_p_brea,2);
                    exp_p_brea(h,:) = exp_p_brea(h,:) / total(h);
                    exp_p_brea_sub(h) = exp_p_brea(h,floor((1-shift)*(range(h)/2)));
                end
            end
            
        elseif ( mod(range(h),2)==1 )
            i = ceil(range(h)/2);
            j = i;
            exp_p_brea(h,i)=1;
            for k=floor(range(h)/2):-1:2
                i=i-1;
                j=j+1;
                exp_p_brea(h,i) = exp(-abs(ceil(range(h)/2)-k)/(0.1*floor(range(h)/2)));
                exp_p_brea(h,j) = exp(-abs(ceil(range(h)/2)-k)/(0.1*floor(range(h)/2)));
            end
            j=j+1;
            exp_p_brea(h,1) = 0;
            exp_p_brea(h,j) = 0;
            
            if(FH)
                total = sum(exp_p_brea,2);
                exp_p_brea(h,:) = exp_p_brea(h,:) / total(h);
                exp_p_brea_sub(h) = exp_p_brea(h,ceil((1-shift)*(range(h)/2)));
            else
                if(UD==1)
                    exp_p_brea(h,ceil(range(h)/2)+1:ceil(range(h)/2)+floor(range(h)/2)) = 0;
                    total = sum(exp_p_brea,2);
                    exp_p_brea(h,:) = exp_p_brea(h,:) / total(h);
                    exp_p_brea_sub(h) = exp_p_brea(h,floor((1-shift)*(range(h)/2)));
                else
                    exp_p_brea(h,1:ceil(range(h)/2)) = exp_p_brea(h,ceil(range(h)/2):ceil(range(h)/2)+floor(range(h)/2));
                    exp_p_brea(h,ceil(range(h)/2)+1:ceil(range(h)/2)+floor(range(h)/2)) = 0;
                    total = sum(exp_p_brea,2);
                    exp_p_brea(h,:) = exp_p_brea(h,:) / total(h);
                    exp_p_brea_sub(h) = exp_p_brea(h,floor((1-shift)*(range(h)/2)));
                end
            end
        end
        
    end
    
	
end