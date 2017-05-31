%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Code developed by Hironori Matsumoto
% Last update : 28 February 2017
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%% CREATE TRIANUGLAR UNBROKEN WAVE PRESSURE
%%% INPUT : range ... vertical range of wave pressure
%%% INPUT : FH ...vertical range of breaking wave pressure, true...full range, false...half range  
%%% INPUT : UD ...upper or lower part? when fh_brea is false, 1...Upper part, 2...Lower part
%%% OUTPUT : sstri_p_stan ... pressure shape function
%%% OUTPUT : sstri_p_stan_sub ... submarine pressure shape function
%%% OUTPUT : sstri_stan_dist1 ... Upper limit of pressure shape function
%%% OUTPUT : sstri_stan_dist2 ... Lower limit of pressure shape function

function [sstri_p_stan, sstri_p_stan_sub, sstri_stan_dist1, sstri_stan_dist2] =  sstri_stan_p_distri(range, FH, UD)
   
    %%% Clear ...
    clear sstri_p_stan sstri_p_stan_sub sstri_stan_dist1 sstri_stan_dist2
    range=ceil(range);
    sstri_p_stan = zeros(3,range(3)+1);

    %%% SET UPPER AND LOWER LIMIT
    if(FH)
        sstri_stan_dist1 = ceil(range/2);    
		sstri_stan_dist2 = floor(range/2);
    else
        if (UD==1)
            sstri_stan_dist1 = ceil(range/2);    
    		sstri_stan_dist2 = [0 0 0];
        else
            sstri_stan_dist1 = [0 0 0];    
        	sstri_stan_dist2 = floor(range/2);
        end
    end
 
    %%% CONSTRUCT SHAPE FUNCTION
	for h=1:length(range)
     
        total = 0;
		if (sstri_stan_dist1(h)==sstri_stan_dist2(h))
            i=1;
            j=sstri_stan_dist1(h)+sstri_stan_dist2(h);
            for k=sstri_stan_dist1(h)-1:-1:1
                i=i+1;
                sstri_p_stan(h,i) = ((sstri_stan_dist1(h)-1)-(k-1))/(sstri_stan_dist1(h)-1);
                total = total + sstri_p_stan(h,i);
                j=j-1;
                sstri_p_stan(h,j) = sstri_p_stan(h,i);
                total = total + sstri_p_stan(h,j);
                if (k==1) sstri_p_stan_sub(h) = 1; end
            end
            for k=1:sstri_stan_dist1(h)+sstri_stan_dist2(h)
                sstri_p_stan(h,k) = sstri_p_stan(h,k)/total;
            end
                
        elseif (sstri_stan_dist1(h)~=sstri_stan_dist2(h))
            i=1;
            j=sstri_stan_dist1(h)+sstri_stan_dist2(h);
            for k=sstri_stan_dist1(h)-1:-1:2
                i=i+1;
                sstri_p_stan(h,i) = ((sstri_stan_dist1(h)-1)-(k-1))/(sstri_stan_dist1(h)-1);
                total = total + sstri_p_stan(h,i);
                j=j-1;
                sstri_p_stan(h,j) = sstri_p_stan(h,i);
                total = total + sstri_p_stan(h,j);
            end
            i=i+1;
            sstri_p_stan(h,i) = ((sstri_stan_dist1(h)-1))/(sstri_stan_dist1(h)-1);
            total = total + sstri_p_stan(h,i);
            sstri_p_stan_sub(h) = 1;
            for k=1:sstri_stan_dist1(h)+sstri_stan_dist2(h)
                sstri_p_stan(h,k) = sstri_p_stan(h,k)/total;
            end
   
        end      
	
    end
            
end