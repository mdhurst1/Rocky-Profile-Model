%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Code developed by Hironori Matsumoto
% Last update : 28 February 2017
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%% CREAT WEATHERING EFFICACY FUNCTION BY Trenhaile and Kanayay (2005) and Stepehson and Kirk (2000b)
%%% INPUT : num_tidal_range ... tidal range / gridsize
%%% INPUT : start_tide ... starting index of tidal range loop
%%% INPUT : end_tide ... ending index of tidal range loop  
%%% INPUT : current_path ... path to working folder  
%%% INPUT : tidal_values ... tidal range
%%% INPUT : gridsize ... gridsize
%%% OUTPUT : wd1 ... weatehring efficacy function by Trenhaile and Kanayay (2005)
%%% OUTPUT : wd2 ... weatehring efficacy function by Stephenson and Kirk (2000b)
%%% 

function [wd1, wd2]=make_wet_dry_func_v1(num_tidal_range, start_tide, end_tide, current_path, tidal_values, gridsize)

	%%% INITIAL SETTING
    %num =  end_tide -  start_tide + 1;
    num = size(num_tidal_range,2);

    %%% ERROR DISPLAY
    for i=1:num
        tr = num_tidal_range(i);
        if( mod(tr,2)==0 )
        else ('tide/gridsize range must be integer even number')
        end
    end

    
    if( gridsize == 0.1 )
        
        max_wd = max(num_tidal_range);
        wd1 = zeros(max_wd, num);
        wd2 = zeros(max_wd, num);

        for j=1:num
            tide=round(num_tidal_range(j));
            if(tide<=1)
                wd1(1,j)= 1;    wd2(1,j)= 1;
            elseif(tide==2)
                wd1(1,j)= 1/2;  wd1(2,j)= 1/2;
                wd2(1,j)= 1/2;  wd2(2,j)= 1/2;
            elseif(tide==3)
                wd1(1,j)= 0;    wd1(2,j)= 1/2;  wd1(3,j)= 0;
                wd2(1,j)= 0;    wd2(2,j)= 1/2;  wd2(3,j)= 0;
            else
                for i=2:ceil(tide/4)                wd1(i,j)=exp(-((i-ceil(tide/4))^2)/(ceil(tide/2)));            end
                for i=ceil(tide/4)+1:tide-1         wd1(i,j)=exp(-((i-ceil(tide/4))^2)/(tide*ceil(tide/10)));      end
                for i=2:ceil(tide/4)				wd2(i,j)=exp(-((i-ceil(tide/4))^2)/(ceil(tide/2)));            end
                for i=ceil(tide/4)+1:ceil(tide*3/4) wd2(i,j)=1;                                                    end
                for i=ceil(tide*3/4)+1:tide-1		wd2(i,j)=exp(-((i-ceil(tide*3/4))^2)/(ceil(tide/2)));          end
                
            end
        end
        
    %%% WHEN GRIDSIZE < 0.1    
    elseif(gridsize < 0.1)
    
        max_wd  = max(num_tidal_range);
        wd1     = zeros(max_wd, num);
        wd2     = zeros(max_wd, num);

        num_tmp             = num; 
        num_tidal_range_tmp = round(tidal_values/0.1);
        max_wd_tmp          = max(num_tidal_range_tmp);
        wd1_tmp             = zeros(max_wd_tmp, num);
        wd2_tmp             = zeros(max_wd_tmp, num);
        
        for j=1:num_tmp
            tide=round(num_tidal_range_tmp(j));
            if(tide<=1)
                wd1_tmp(1,j)= 1;    wd2_tmp(1,j)= 1;
            elseif(tide==2)
                wd1_tmp(1,j)= 1/2;  wd1_tmp(2,j)= 1/2;
                wd2_tmp(1,j)= 1/2;  wd2_tmp(2,j)= 1/2;
            elseif(tide==3)
                wd1_tmp(1,j)= 0;    wd1_tmp(2,j)= 1/2;  wd1_tmp(3,j)= 0;
                wd2_tmp(1,j)= 0;    wd2_tmp(2,j)= 1/2;  wd2_tmp(3,j)= 0;
            else
                for i=2:ceil(tide/4)                wd1_tmp(i,j)=exp(-((i-ceil(tide/4))^2)/(ceil(tide/2)));            end
                for i=ceil(tide/4)+1:tide-1         wd1_tmp(i,j)=exp(-((i-ceil(tide/4))^2)/(tide*ceil(tide/10)));      end
                %sum_wd1 = sum(wd1(:,j));
                %wd1(:,j) = wd1(:,j)/sum_wd1;
                for i=2:ceil(tide/4)				wd2_tmp(i,j)=exp(-((i-ceil(tide/4))^2)/(ceil(tide/2)));            end
                for i=ceil(tide/4)+1:ceil(tide*3/4) wd2_tmp(i,j)=1;                                                    end
                for i=ceil(tide*3/4)+1:tide-1		wd2_tmp(i,j)=exp(-((i-ceil(tide*3/4))^2)/(ceil(tide/2)));          end
                %sum_wd2 = sum(wd2(:,j));
                %wd2(:,j) = wd2(:,j)/sum_wd2;
        
            end
        end
        
        %%% INTERPOLATION
        if(gridsize==0.05)      add=2;
        elseif(gridsize==0.075) add=1;
        end
            
        for j=1:num 
            wd1(:,j) = interp1([1:max_wd_tmp]',wd1_tmp(:,j),[1:(max_wd_tmp/(max_wd+add)):max_wd_tmp]');
            wd2(:,j) = interp1([1:max_wd_tmp]',wd2_tmp(:,j),[1:(max_wd_tmp/(max_wd+add)):max_wd_tmp]');
        end
            
    end
    
    %%% PRINT
    if(false)
    for j=1:num
        tide=round(num_tidal_range(j));
        if(true)
            figure; 
            x=[1:tide];
            print_wd1 = wd1(:,j);
            print_wd2 = wd2(:,j);
            if (tide < max(num_tidal_range))
                print_wd1(tide+1:max(num_tidal_range)) = [];
                print_wd2(tide+1:max(num_tidal_range)) = [];
            end
            plot(x, print_wd1, x, print_wd2, 'linewidth', 2);	
            filename = ['\', num2str(tidal_values(j)), '-wet&dry.jpg'];
            filename = [current_path, filename];
            saveas(gcf, filename);
            close;
        end
    end
    end
    
end