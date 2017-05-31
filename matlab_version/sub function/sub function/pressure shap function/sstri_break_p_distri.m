%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Code developed by Hironori Matsumoto
% Last update : 28 February 2017
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%% CREATE TRIANGULAR BREAKING WAVE PRESSURE
%%% INPUT : range ... vertical range of wave pressure
%%% INPUT : shift ... vertical shift of wave pressure
%%% INPUT : FH ...vertical range of breaking wave pressure, true...full range, false...half range  
%%% INPUT : UD ...upper or lower part? when fh_brea is false, 1...Upper part, 2...Lower part
%%% OUTPUT : sstri_p_brea ... pressure shape function
%%% OUTPUT : sstri_p_brea_sub ... submarine pressure shape function
%%% OUTPUT : sstri_brea_dist1 ... Upper limit of pressure shape function
%%% OUTPUT : sstri_brea_dist2 ... Lower limit of pressure shape function

function [sstri_p_brea, sstri_p_brea_sub, sstri_brea_dist1, sstri_brea_dist2] =  sstri_break_p_distri(range, shift, FH, UD)

    %%% Clear ...
    clear sstri_p_brea sstri_p_brea_sub sstri_brea_dist1 sstri_brea_dist2
    range=ceil(range);
    sstri_p_brea = zeros(3,range(3)+1);
    
    %%% SET UPPER AND LOWER LIMIT
    if (FH)
        sstri_brea_dist1 = ceil((range/2)*(1+shift));
		sstri_brea_dist2 = floor((range/2)*(1-shift));
    else
        if (UD==1)
            sstri_brea_dist1 = ceil((range/2)*(1+shift));
            sstri_brea_dist2 = -(sstri_brea_dist1 - ceil(range/2));
        else
            sstri_brea_dist2 = floor((range/2)*(1-shift));
            sstri_brea_dist1 = -sstri_brea_dist2 + floor(range/2);
        end
    end
    
    %%% CONSTRUCT SHAPE FUNCTION
	for h=1:length(range) 
        total=0;
        
        if ( mod(range(h),2)==0 )
            i = range(h)/2+1;
            j = range(h)/2;
            for k=range(h)/2:-1:1
                i=i-1;
                j=j+1;
                sstri_p_brea(h,i) = (k-1)/(range(h)/2-1);
                sstri_p_brea(h,j) = sstri_p_brea(h,i);
            end
            
            if(FH)           
                total = sum(sstri_p_brea,2);
                sstri_p_brea(h,:) = sstri_p_brea(h,:)/total(h);
                if(floor((1-shift)*(range(h)/2))==0)
                    sstri_p_brea_sub(h) = sstri_p_brea(h,ceil((1-shift)*(range(h)/2)));
                else
                    sstri_p_brea_sub(h) = sstri_p_brea(h,floor((1-shift)*(range(h)/2)));
                end
            else
                if(UD==1)
                    sstri_p_brea(h,range(h)/2+1:range(h)) = 0;
                    total = sum(sstri_p_brea,2);
                    sstri_p_brea(h,:) = sstri_p_brea(h,:) / total(h);
                    if(floor((1-shift)*(range(h)/2))==0)
                        sstri_p_brea_sub(h) = sstri_p_brea(h,ceil((1-shift)*(range(h)/2)));
                    else
                        sstri_p_brea_sub(h) = sstri_p_brea(h,floor((1-shift)*(range(h)/2)));
                    end
                elseif (UD==2)
                    sstri_p_brea(h,1:range(h)/2) = sstri_p_brea(h,range(h)/2+1:range(h)); 
                    sstri_p_brea(h,range(h)/2+1:range(h)) = 0;
                    total = sum(sstri_p_brea,2);
                    sstri_p_brea(h,:) = sstri_p_brea(h,:) / total(h);
                    if(floor((1-shift)*(range(h)/2))==0)
                        sstri_p_brea_sub(h) = sstri_p_brea(h,ceil((1-shift)*(range(h)/2)));
                    else
                        sstri_p_brea_sub(h) = sstri_p_brea(h,floor((1-shift)*(range(h)/2)));
                    end
                end
            end     
                    
        elseif ( mod(range(h),2)==1 )
            
            i = ceil(range(h)/2);
            j = i;
            sstri_p_brea(h,i)=i/ceil(range(h)/2);
            for k=1:ceil(range(h)/2)-1
                i=i-1;
                j=j+1;
                sstri_p_brea(h,i) = ((ceil(range(h)/2)-1)-k) / (ceil(range(h)/2)-1);
                sstri_p_brea(h,j) = sstri_p_brea(h,i);
            end
            if(FH)
                total = sum(sstri_p_brea,2);
                sstri_p_brea(h,:) = sstri_p_brea(h,:)/total(h);
                if(floor((1-shift)*(range(h)/2))==0)
                    sstri_p_brea_sub(h) = sstri_p_brea(h,ceil((1-shift)*(range(h)/2)));
                else
                    sstri_p_brea_sub(h) = sstri_p_brea(h,floor((1-shift)*(range(h)/2)));
                end
            else
                if (UD==1)
                    sstri_p_brea(h,ceil(range(h)/2)+1:range(h)) = 0;
                    total = sum(sstri_p_brea,2);
                    sstri_p_brea(h,:) = sstri_p_brea(h,:) / total(h);
                    if(floor((1-shift)*(range(h)/2))==0)
                        sstri_p_brea_sub(h) = sstri_p_brea(h,ceil((1-shift)*(range(h)/2)));
                    else
                        sstri_p_brea_sub(h) = sstri_p_brea(h,floor((1-shift)*(range(h)/2)));
                    end
            elseif (UD==2)
                    sstri_p_brea(h,1:ceil(range(h)/2)) = sstri_p_brea(h,ceil(range(h)/2):range(h));
                    sstri_p_brea(h,ceil(range(h)/2)+1:range(h)) = 0;
                    total = sum(sstri_p_brea,2);
                    sstri_p_brea(h,:) = sstri_p_brea(h,:) / total(h);
                    if(floor((1-shift)*(range(h)/2))==0)
                        sstri_p_brea_sub(h) = sstri_p_brea(h,ceil((1-shift)*(range(h)/2)));
                    else
                        sstri_p_brea_sub(h) = sstri_p_brea(h,floor((1-shift)*(range(h)/2)));
                    end
                end
            end
        end
            
	end
	
end