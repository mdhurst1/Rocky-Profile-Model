%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Code developed by Hironori Matsumoto
% Last update : 28 February 2017
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%% CREATE TIDAL DURATION FUNCTION
%%% INPUT : num_tidal_range ... tidal range / gridsize
%%% INPUT : start_tide ... starting index of tidal range loop
%%% INPUT : end_tide ... ending index of tidal range loop  
%%% INPUT : gridsize ... gridsize
%%% OUTPUT : tr_esf ... tidal duration function

function [tr_esf]=make_tidal_range(num_tidal_range, start_tide, end_tide, gridsize)

	%%% INITIAL SETTING
    total   = 0;
    num =  size(num_tidal_range,2);
    tr_max  = max(num_tidal_range);
    tr_esf  = zeros(tr_max,num);

    %%% CHANGE BOTTOM DEPENDING ON GRIDSIZE
    if ( gridsize == 0.1 )  bottom = 1;
    else                    bottom  = gridsize / 0.1;
    end
    
    
    %%% ERROR DISPLAY
    for i=1:num
        if ( ceil(num_tidal_range(i))==floor(num_tidal_range(i)) )
        else ('tide/gridsize range must be integer number')
            num_tidal_range(i)
            stop
        end
    end
    
    %%% DEVELOP
    for j=1:num
        total=0;
        tr = num_tidal_range(j);
        if(tr <= 2)
            for i=1:tr
                tr_esf(i,j) = 1;
                total = total + tr_esf(i,j)*bottom;
            end
            tr_esf(:,j) = tr_esf(:,j)/total;
        elseif ( ceil(tr)==3 )
            tr_esf(1,j) = 0;    tr_esf(2,j) = 1/bottom;    tr_esf(3,j) = 0;
        elseif ( ceil(tr)==4 )
            tr_esf(1,j) = 0;
            tr_esf(2,j) = 1;    total=total+tr_esf(2,j)*bottom;
            tr_esf(3,j) = 1;    total=total+tr_esf(3,j)*bottom;
            tr_esf(4,j) = 0;
            tr_esf(:,j) = tr_esf(:,j)/total;
        elseif ( tr > 4 & tr < 20)
            if(ceil(0.5*tr)~=floor(0.5*tr))
                for i=1:ceil(0.5*tr)+1      
                    tr_esf(i,j) = sin((i-1)*pi/ceil(tr*0.5));
                    total = total + tr_esf(i,j)*bottom;
                end
                for i=floor(0.5*tr):tr     
                    tr_esf(i,j) = tr_esf(i,j) + ... 
                            sin((i-floor(tr*0.5))*pi/ceil(tr*0.5));
                    total = total + tr_esf(i,j)*bottom;
                end
            else
                for i=1:0.5*tr+1      
                    tr_esf(i,j) = sin((i-1)*pi/(0.5*tr));
                    total = total + tr_esf(i,j)*bottom;
                end
                for i=0.5*tr:tr
                    tr_esf(i,j) = tr_esf(i,j) + ... 
                            sin((i-tr*0.5)*pi/(tr*0.5));
                    total = total + tr_esf(i,j)*bottom;
                end
            end
            tr_esf(:,j) = tr_esf(:,j)/total;
            
        elseif ( tr >= 20)
            for i=1:ceil(0.55*tr)+1     
                tr_esf(i,j) = sin((i-1)*pi/ceil(tr*0.55));
                total = total + tr_esf(i,j)*bottom;
            end
            for i=floor(0.45*tr):tr     
                tr_esf(i,j) = tr_esf(i,j) + ... 
                            sin((i-(tr-floor(tr*0.55)))*pi/ceil(tr*0.55));
                total = total + tr_esf(i,j)*bottom;
            end
            tr_esf(:,j) = tr_esf(:,j)/total;        
        end
            
    end
end