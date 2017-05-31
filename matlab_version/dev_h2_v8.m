%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Code developed by Hironori Matsumoto
% Last update : 28 February 2017
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


clear all
close all


%%%
%%% START OF CALCULATION 
%%%
startT = clock();        % Clock




%%%
%%% GRID LOOP USED TO TEST GRID SIZE SENSITIVITY
%%%
gridsize_values = [0.1];


%%% FOR TESTING GRIDSIZE SENSITIVITY IN GEOMORPH PAPER (Matsumoto et al., 2016)
for gridloop = 1:length(gridsize_values)
    
    
clearvars -except gridloop gridsize_values sec
clf



%%%
%%% FLAG FOR PRINTING
%%%
final_prof_to_txt = false;
print_wvsw_to_fig = false;
print_wvsw_to_txt = false;
print_break_point = false;
print_prof_to_txt = true;
print_weff_to_jpg = false;
print_waty_to_txt = false;
print_pcom_to_jpg = true;
print_tedi_to_jpg = false;
print_snap_to_mov = false;
print_bwer_to_txt = false;
print_nxbw_to_txt = false;



%%%
%%% GRIDSIZE
%%%
gridsize = gridsize_values(gridloop);


%%%  
%%% Process selection :
%%%
bw              = true;		% Back-wearing erosion
dw              = true;		% Down-wearing erosion
weathering      = true;		% Weathering
spr_weathering  = false;    % Supla tidal weathering
rsl_flag        = true;     % Sea level change
tec_flag        = true;     % Tectonic event
pre_platform    = false;    % Occurrence of initial platform



%%% 
%%% SET CONTROLLING VARIABLES : 
%%%
   
maxit		= 8000;                         % MAXIMUM ITERATION
maxcell  	= 500;                          % CALCULATION SPACE IN METRE
x_maxcell	= ceil(maxcell/gridsize);       % HORIZONTAL MAXCELL 
y_maxcell	= ceil(maxcell/20/gridsize);    % VERTICAL MAXCELL

resi_values     = [5e-2 5e-1 1e-0];         % ROCK RESISTANCE VALUE
num_resi_values = resi_values*gridsize;     % CONVERT TO NUMMERICAL RESISTANCE VALUE
start_resi      = 1;                        % STATING IDEX OF RESISTANCE LOOP
end_resi        = 1;                        % ENDING INDEX OF RESISTANCE LOOP

height_values   = [1 2 1.5 2];              % WAVE HEIGHT VALUE
num_height      = round(height_values / gridsize);  % CONVERT TO NUMMERICAL WAVE HEIGHT VALUE
start_height    = 1;                        % STARTING INDEX OF WAVE HEIGHT LOOP
end_height      = 2;                        % ENDING INDEX OF WAVE HEIGHT LOOP

wea_const2      = 5e-2;                     % WEATHERING EFFICACY VALUE 2
wea_const3      = 5e-3;                     % WEATHERING EFFICACY VALUE 3
wea_const4      = 1e-3;                     % WEATHERING EFFICACY VALUE 4

cwea_const2     = 1e-1;                     % CLIFF WEATHERING EFFICACY VALUE 2 
cwea_const3     = 1e-2;                     % CLIFF WEATHERING EFFICACY VALUE 3
cwea_const4     = 5e-3;                     % CLIFF WEATHERING EFFICACY VALUE 4

tidal_values    = [1 2 4 8];                % MEAN TIDAL RANGE
num_tidal_range = round(tidal_values / gridsize);   % CONVERT TO NUMMERICAL MEAN TIDAL RANGE VALUE
for gg=1:size(tidal_values,2)               % IN CASE OF FINER GRIDSIZE ...
    if (num_tidal_range(gg)<1)              
        num_tidal_range(gg)=1; 
    end 
end 
tidal_name  = ['1' '2' '4' '8'];            % STR FOR TIDAL RANGE VALIATION
start_tide  = 1;                            % STARTING INDEX IF TIDAL RANGE LOOP
end_tide    = 1;                            % ENDING INDEX OF TIDAL RANGE LOOP

rsl = [0 2 2 2 2 1.5 1.0 0.5 0];                        % CHANGE RSL FOR 8000 ITERATION
event = [1 1001 2001 3001 4001 5001 6001 7001 8001];    % TIMING OF RSL CHANGE FOR 8000 ITERATION
tectonic    = [2 2 2 2 2];                              % SEA LEVEL CHANGE DUE TO TECTONIC EVENT ( POSITIVE:UPLIFT, NEGATIVE:DOWNLIFT ) 
ttic        = 8000-[4500 3500 1900 1600 250];           % TIMING OF TECTONIC EVENT
store_mwl   = zeros(1,size(tectonic,1));                % VARIABLE FOR STORING MWL



%%% 
%%% NON-CONTROLLING VARIABLES : 
%%% see Matsumoto et al. 2016 for details
%%%

sub_decay_const = 1e-1;         % SUBMARINE DECAY CONST
stan_const  = 0.01;             % CONSTANT FOR UNBROKEN WAVE   
brea_const  = 10;               % CONSTANT FOR BREAKING WAVE
brok_const  = 0.1;              % CONSTANT FOR BROKEN WAVE
swd1        = 1e-1;             % WAVE HEIGHT DECAY FOR BREAKING WAVE
swd2        = 1e-2;             % WAVE HEIGHT DECAY FOR BROKEN WAVE
                                

%%% FLAG FOR MAX EROSION        % WHETHER WAVE ERODE MORE THAN ONE ROCK CELL
max_ero_flag= false;
change    = 2;



%%% variables for multiple erosion
max_bw_erosion  = 1*0.1/gridsize_values(gridloop);      % Flag for dw multiple erosion
max_dw_erosion  = 1*0.1/gridsize_values(gridloop);      % Flag for dw multiple erosion
max_we_erosion  = 1*0.1/gridsize_values(gridloop);      % Flag for dw multiple erosion
max_we_erosion = 1;                                     % TEMPORARY  



%%%  
%%% Other : 
%%%

%%% STR
profile =   'p';	
connect =   '-';	
t_range =   'tide';  
   
%%% ANIMATION TIMING
if(gridsize==0.01)
	print_loop1 = 1;
	print_loop2 = 2;
    print_loop3 = 100;
else
    print_loop1 = 2;
	print_loop2 = 10;
    print_loop3 = 100;
end    
	
%%% CURRENT_PATH(=pwd) DEFINITION
current_path = 'C:\Users\hmat258\Desktop\HIRO\UoA\Research\ROCKY MODEL\Matlab\Source&Result\PhD thesis\Tectonic\';
%current_path = '/home/hmat258/project/matlab/rocky/sealevel/';
%current_path = [current_paith1, num2str(gridsize), '\'];
addpath(genpath(current_path)); % ADD SEARCH PATH


%%% ALLOCATE MEMORY SPACE FOR CALCULATION SPACE 
mat 	= ones(x_maxcell,y_maxcell);

%%% EROSION FLAG AND RESISTANCE
field1  = 'flag';                       % STRUCTURE MEMBER FOR EROSION		 
field2  = 'resistance';                 % STRUCTURE MEMBER FOR MATERIAL
value1  = ones(x_maxcell,y_maxcell);
value2  = ones(x_maxcell,y_maxcell);
s       = struct(field1,value1,field2,value2);
clear value1 value2


%%%   
%%% EROSION SHAPE FUNCTION
%%%
[tr_esf]=make_tidal_range(num_tidal_range,start_tide,end_tide, gridsize);

%%%   
%%% WETTING & DRYING SPATIAL FUNCTION
%%% wd1 ... Trenhaile and Kanyaya, 2005
%%% wd2 ... Stephenson and Kirk, 2000b
[wd1, wd2] = make_wet_dry_func_v1(num_tidal_range, start_tide, end_tide, current_path, tidal_values, gridsize);




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% %%%%%%%%% %%%%%%%%%%%  %%%%%%%%%%%%%%     %%%%%%%%%% %%%%%%%%% %%%%%
%%%%% % %%%%% % %%%%%%%%%% %% %%%%%%%%%%%%%%% %%%%%%%%%%%% % %%%%%%% %%%%%
%%%%% %% %%% %% %%%%%%%%% %%%% %%%%%%%%%%%%%% %%%%%%%%%%%% %% %%%%%% %%%%%
%%%%% %%% % %%% %%%%%%%% %%%%%% %%%%%%%%%%%%% %%%%%%%%%%%% %%% %%%%% %%%%%
%%%%% %%%% %%%% %%%%%%%          %%%%%%%%%%%% %%%%%%%%%%%% %%%% %%%% %%%%%
%%%%% %%%%%%%%% %%%%%% %%%%%%%%%% %%%%%%%%%%% %%%%%%%%%%%% %%%%% %%% %%%%%
%%%%% %%%%%%%%% %%%%% %%%%%%%%%%%% %%%%%%%%%% %%%%%%%%%%%% %%%%%% %% %%%%%
%%%%% %%%%%%%%% %%%% %%%%%%%%%%%%%% %%%%%%%     %%%%%%%%%% %%%%%%% % %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%
%%% PROFILE LOOP
%%% 1=90[deg], 2=85[deg], 3=45[deg], 4=20[deg]
for prof=5:5


    %%% SET INITIAL PROFILE
	if (prof==1)     i_angle = 90;      
	elseif (prof==2) i_angle = 85;	
	elseif (prof==3) i_angle = 75;
	elseif (prof==4) i_angle = 67.5;
	elseif (prof==5) i_angle = 45;
	elseif (prof==6) i_angle = 20;
    end

    
    %%% PROFILE INITIALIZING
    %%% 1: LAND CELL
    %%% 2: SEA/AIR CELL
    if( i_angle~= 90 )
        for y=1:y_maxcell
            x_dist = round(y/tan(i_angle*pi/180));
            for x=1:x_dist mat(x,y)=0; end
        end
    end
    
    
    %%%	
    %%% TIDE LOOP
    %%% 
    for tr=start_tide:end_tide

        %%% TIDE RELATED SETTING
        tide        = num_tidal_range(tr);
        mwl         = round(y_maxcell/2);
        reset_mwl   = mwl;
        
        %%% RESTORE EROSION SHAPE FUNCTION
        clear esf
        esf         = tr_esf(:,tr);
        wetdry_a    = wd1(:,tr);    %%% WEATHERING SHAPE FUNCTION BY TRENHAILE AND KANYAYA (2005) 
        wetdry_b    = wd2(:,tr);    %%% WEATHERING SHAPE FUNCTION BY STEPHENSON AND KIRK (2000b)
        
        %%% CREATE PRE EXISTING PLATFOMR << NOT USED 
        if(pre_platform) 
			for y=mwl-tide/2:y_maxcell
				mat(1:4000,y) = 0;
			end
        end
        
           
		%%% 
		%%% PRESSURE DISTRIBUTION
		%%%
        %%% p_stan ... unbroken wave pressure
		%%% p_breaking ... breaking wave pressure
		%%% p_broken ... broken wave pressure (p_broken=1...Homma & Horikawa, p_broken=2...CERC )
        %%% 
        %%% PRESSURE DISTRIBUTION SENSITIVITY TEST IN GEOMORPH PAPER (Matsumoto et al., 2016)
        %%%
        %%% range...vertidal range of wave pressure
        %%% shift...vertical shift of wave pressure, 0... no shift, 0.5...shift by half the wave height
        %%% fh_brea...vertical range of breaking wave pressure, true...full range of vertical breaking wave pressure, false...half range of vertical breaking wave pressure 
        %%% fh_brok...vertical range of broken wave pressure, true...full range of vertical broken wave pressure, false...half range of vertical broken wave pressure 
        %%% ud_brea.. consideration of upper or lower part of breaking wave pressure when fh_brea is false
        %%% ud_brok.. consideration of upper or lower part of broken wave pressure when fh_brok is false
        
        for pressure=1:1
        
            %%% SET SHAPE OF PRESSURE : 
            clear p_stan p_break p_brok
            clear stan_dist1 stan_dist2 break_dist1 break_dist2 brok_dist1 brok_dist2
            clear p_stan_sub p_break_sub p_brok_sub

            fh_brea = true;
            ud_brea = 1;
            ud_brok = 1;
            
            jump=0;
            
            if ( pressure==1 | pressure==2 | pressure==3 | pressure==4 | pressure==5 | pressure==6 | pressure==7 | pressure==8 ) % RECTANGULAR & RECTANGULAR

                if (pressure == 1)      range = height_values/gridsize;     shift = 0;      fh_brok = true; 
                elseif (pressure == 2)  range = height_values*0.5/gridsize; shift = 0;      fh_brok = true;
                elseif (pressure == 3)  range = height_values/gridsize;     shift = 0.5;    fh_brok = true;
                elseif (pressure == 4)  range = height_values*0.5/gridsize; shift = 0.5;    fh_brok = true;
                elseif (pressure == 5)  range = height_values/gridsize;     shift = 0;      fh_brok = false;  
                elseif (pressure == 6)  range = height_values*0.5/gridsize; shift = 0;      fh_brok = false;
                elseif (pressure == 7)  range = height_values/gridsize;     shift = 0.5;    fh_brok = false;
                elseif (pressure == 8)  range = height_values*0.5/gridsize; shift = 0.5;    fh_brok = false;
                end
               
                [p_stan, p_stan_sub, stan_dist1, stan_dist2]        = ssrec_stan_p_distri(range, true, 1);
                [p_break, p_break_sub, break_dist1, break_dist2]    = ssrec_break_p_distri(range, shift, fh_brea, ud_brea);
                [p_brok, p_brok_sub, brok_dist1, brok_dist2]        = ssrec_brok_p_distri(range, fh_brok, ud_brok);
            

            elseif ( pressure==9 | pressure==10 | pressure==11 | pressure==12 | pressure==13 | pressure==14 | pressure==15 | pressure==16) % RECTANGULAR & TRIANGLE
            
                if (pressure == 9)      range = height_values/gridsize;     shift = 0;      fh_brok = true; 
                elseif (pressure == 10) range = height_values*0.5/gridsize; shift = 0;      fh_brok = true;
                elseif (pressure == 11) range = height_values/gridsize;     shift = 0.5;    fh_brok = true;
                elseif (pressure == 12) range = height_values*0.5/gridsize; shift = 0.5;    fh_brok = true;
                elseif (pressure == 13) range = height_values/gridsize;     shift = 0;      fh_brok = false;  
                elseif (pressure == 14) range = height_values*0.5/gridsize; shift = 0;      fh_brok = false;
                elseif (pressure == 15) range = height_values/gridsize;     shift = 0.5;    fh_brok = false;
                elseif (pressure == 16) range = height_values*0.5/gridsize; shift = 0.5;    fh_brok = false;
                end
                
                [p_stan, p_stan_sub, stan_dist1, stan_dist2]        = ssrec_stan_p_distri(range, true, 1);
                [p_break, p_break_sub, break_dist1, break_dist2]    = ssrec_break_p_distri(range, shift, fh_brea, ud_brea);
                [p_brok, p_brok_sub, brok_dist1, brok_dist2]        = sstri_brok_p_distri(range, fh_brok,ud_brok);

                
            elseif ( pressure==17 | pressure==18 | pressure==19 | pressure==20 | pressure==21 | pressure==22 | pressure==23 | pressure==24 ) % TRIANGLE & RECTANGULAR

                if (pressure == 17)      range = height_values/gridsize;     shift = 0;      fh_brok = true; 
                elseif (pressure == 18)  range = height_values*0.5/gridsize; shift = 0;      fh_brok = true;
                elseif (pressure == 19)  range = height_values/gridsize;     shift = 0.5;    fh_brok = true;
                elseif (pressure == 20)  range = height_values*0.5/gridsize; shift = 0.5;    fh_brok = true;
                elseif (pressure == 21)  range = height_values/gridsize;     shift = 0;      fh_brok = false;  
                elseif (pressure == 22)  range = height_values*0.5/gridsize; shift = 0;      fh_brok = false;
                elseif (pressure == 23)  range = height_values/gridsize;     shift = 0.5;    fh_brok = false;
                elseif (pressure == 24)  range = height_values*0.5/gridsize; shift = 0.5;    fh_brok = false;
                end
                
                [p_stan, p_stan_sub, stan_dist1, stan_dist2]        = ssrec_stan_p_distri(range, true, 1);
                [p_break, p_break_sub, break_dist1, break_dist2]    = sstri_break_p_distri(range, shift, fh_brea, ud_brea);
                [p_brok, p_brok_sub, brok_dist1, brok_dist2]        = ssrec_brok_p_distri(range, fh_brok,ud_brok);
            
            elseif ( pressure==25 | pressure==26 | pressure==27 | pressure==28 | pressure==29 | pressure==30 | pressure==31 | pressure==32) % TRIANGLE & TRIANGLE
            
                if (pressure == 25)      range = height_values/gridsize;    shift = 0;      fh_brok = true; 
                elseif (pressure == 26) range = height_values*0.5/gridsize; shift = 0;      fh_brok = true;
                elseif (pressure == 27) range = height_values/gridsize;     shift = 0.5;    fh_brok = true;
                elseif (pressure == 28) range = height_values*0.5/gridsize; shift = 0.5;    fh_brok = true;
                elseif (pressure == 29) range = height_values/gridsize;     shift = 0;      fh_brok = false;  
                elseif (pressure == 30) range = height_values*0.5/gridsize; shift = 0;      fh_brok = false;
                elseif (pressure == 31) range = height_values/gridsize;     shift = 0.5;    fh_brok = false;
                elseif (pressure == 32) range = height_values*0.5/gridsize; shift = 0.5;    fh_brok = false;
                end
                
                [p_stan, p_stan_sub, stan_dist1, stan_dist2]        = ssrec_stan_p_distri(range, true, 1);
                [p_break, p_break_sub, break_dist1, break_dist2]    = sstri_break_p_distri(range, shift, fh_brea, ud_brea);
                [p_brok, p_brok_sub, brok_dist1, brok_dist2]        = sstri_brok_p_distri(range, fh_brok,ud_brok);

            elseif ( pressure==33 | pressure==34 | pressure==35 | pressure==36 | pressure==37 | pressure==38 | pressure==39 | pressure==40 ) % EXPPONENTIAL & RECTANGULAR

                if (pressure == 33)      range = height_values/gridsize;     shift = 0;      fh_brok = true; 
                elseif (pressure == 34)  range = height_values*0.5/gridsize; shift = 0;      fh_brok = true;
                elseif (pressure == 35)  range = height_values/gridsize;     shift = 0.5;    fh_brok = true;
                elseif (pressure == 36)  range = height_values*0.5/gridsize; shift = 0.5;    fh_brok = true;
                elseif (pressure == 37)  range = height_values/gridsize;     shift = 0;      fh_brok = false;  
                elseif (pressure == 38)  range = height_values*0.5/gridsize; shift = 0;      fh_brok = false;
                elseif (pressure == 39)  range = height_values/gridsize;     shift = 0.5;    fh_brok = false;
                elseif (pressure == 40)  range = height_values*0.5/gridsize; shift = 0.5;    fh_brok = false;
                end
                
                [p_stan, p_stan_sub, stan_dist1, stan_dist2]        = ssrec_stan_p_distri(range, true, 1);
                [p_break, p_break_sub, break_dist1, break_dist2]    = ori_break_p_distri(range, shift, fh_brea, ud_brea);
                [p_brok, p_brok_sub, brok_dist1, brok_dist2]        = ssrec_brok_p_distri(range, fh_brok,ud_brok);
            
            elseif ( pressure==41 | pressure==42 | pressure==43 | pressure==44 | pressure==45 | pressure==46 | pressure==47 | pressure==48) % EXPONENTIAL & TRIANGLE
            
                if (pressure == 41)     range = height_values/gridsize;     shift = 0;      fh_brok = true; 
                elseif (pressure == 42) range = height_values*0.5/gridsize; shift = 0;      fh_brok = true;
                elseif (pressure == 43) range = height_values/gridsize;     shift = 0.5;    fh_brok = true;
                elseif (pressure == 44) range = height_values*0.5/gridsize; shift = 0.5;    fh_brok = true;
                elseif (pressure == 45) range = height_values/gridsize;     shift = 0;      fh_brok = false;  
                elseif (pressure == 46) range = height_values*0.5/gridsize; shift = 0;      fh_brok = false;
                elseif (pressure == 47) range = height_values/gridsize;     shift = 0.5;    fh_brok = false;
                elseif (pressure == 48) range = height_values*0.5/gridsize; shift = 0.5;    fh_brok = false;
                end
                
                [p_stan, p_stan_sub, stan_dist1, stan_dist2]        = ssrec_stan_p_distri(range, true, 1);
                [p_break, p_break_sub, break_dist1, break_dist2]    = ori_break_p_distri(range, shift, fh_brea, ud_brea);
                [p_brok, p_brok_sub, brok_dist1, brok_dist2]        = sstri_brok_p_distri(range, fh_brok,ud_brok); 
            
            elseif ( pressure==49 | pressure==50 ) % NO SPATIAL VARIATION IN WAVE PRESSURE 
                if(pressure==49)    jump = 0;   shift = 0;
                elseif(pressure==50)jump = 0.5; shift = 0;
                end
                range = height_values/gridsize;
                p_stan   = (1./range)';   p_stan_sub = 1./range;  stan_dist1     = [1 1 1] ;     stan_dist2    = [0 0 0];
                p_break  = (1./range)';   p_break_sub= 1./range;  break_dist1    = [1 1 1] ;     break_dist2   = [0 0 0];
                p_brok   = (1./range)';   p_brok_sub = 1./range;  brok_dist1     = [1 1 1] ;     brok_dist2    = [0 0 0];
                
            end
            clear range bd
			
            
            
            %%%
            %%% Weathering LOOP  
            %%% 

            %%% WEATHERING EFFICACY TYPE >> USE WEA 1
            for wea=1:1
				if ( wea == 1 ) weathering = true;  wetdry=wetdry_a;                % 
				elseif (wea==2) weathering = true;  wetdry=wetdry_a;    wea_flag=1; % NO WEATHEIRNG BUT ONLY INITIAL LOOP
                elseif (wea==3) weathering = false; wetdry=wetdry_a;                % NO WEATHEIRNG
                end
            
                %%% WEATHERING EFFICACY COEFFICIENT 
                for weasd=2:2
                if ( weasd == 1 )       wt_const = wea_const1;  % wet & dry weathering const;
                elseif (weasd == 2)     wt_const = wea_const2;  % wet & dry weathering const;
                elseif (weasd == 3)     wt_const = wea_const3;	% wet & dry weathering const;
                elseif (weasd == 4)     wt_const = wea_const4;  % wet & dry weathering const;
                elseif (weasd == 5)     wt_const = wea_const5;  % wet & dry weathering const;
                end
                
                
                %%% TECTONIC EVENT
            	for co=1:size(tectonic,2)
                    
                    
                    %%% FOR CREATION OF FOLDER                
                    dir_path = [current_path, num2str(i_angle), '/tide-', tidal_name(tr), '/pressure', num2str(pressure), ...
							'/wea', num2str(wea), '/weasd', num2str(weasd), '/co', num2str(co)];
							%'\co', num2str(co), '\gs', num2str(gridsize_values(gridloop))];
                    mkdir(dir_path);
                    output_path = [dir_path, '/']; 

            
                    %%%
                    %%% WAVE HEIGHT LOOP
                    %%% 
                    for h=start_height:end_height
            
                        %%% CLEAR WAVE HEIGHT RELTING VARIABLE
                        clear height bd decay_sw1
                        
                        %%% WAVE EROSION FORCE IS NOW PROPORTIONAL TO THE SQUARE OF WAVE HEIGHT 
                        height = [height_values(h)*height_values(h)];
                        wave_height = [num2str(height_values(h)), 'm'];     %%% STR...

                        %%% 
                        break_wave_dist     = num_height(h)/2;                      % HORIZONTAL DISTANCE OF BREAKING WAVE RANGE 
                        bd                  = round(num_height(h)/0.78);            % BREAKING DEPTH : Hb / hb = 0.78
                        decay_sw1           = -log(swd1)/(break_wave_dist);         % WAVE ATTENUATION CONSTANT1
                   
                        %%% WAVE HEIGHT DECAY FROM BREAKING TO BROKEN WAVE
                        bwd = 0.1;
                
                        %%% INITIALIZATION OF FINAL PROF FOR RESULT PROTTING
                        final_prof  =   zeros(y_maxcell,size(resi_values,2));
                                
	
                        %%%  
                        %%% MATERIAL RESISTANCE LOOP  
                        %%%
                        for resi=start_resi:end_resi
				
                            
                            s.flag              = mat;
                            s.resistance        = mat*num_resi_values(resi);
                            effec_resi          = num_resi_values(resi);
                            ix_posi=1;
                            while(s.flag(ix_posi,end)==0) ix_posi=ix_posi+1; end
                            ix_max              = ix_posi;
                    
                            %%% STR OF MATERIAL RESISNTACE
                            material_resistance = ['resi',num2str(resi_values(resi))];
                    
                            
                            
                            %%%
                            %%% VARIABLES FOR PLOTTING
                            %%% 
                            bw_erosion		    = zeros(y_maxcell,1);       % BACK WEARING EROSION
                            total_bw_erosion	= zeros(y_maxcell,maxit);	% TOTAL BACK-WEARING EROSION
                            dw_erosion		    = zeros(x_maxcell,1);       % DOWN WEARING EROSION
                            total_dw_erosion	= zeros(x_maxcell,maxit);	% TOTAL DOWN WEARING ERSION
                            wea_erosion		    = zeros(y_maxcell,1);       % WEATHERING EROSION
                            total_wea_erosion	= zeros(y_maxcell,maxit);	% TOTAL WEATHERING EROSION

                            tidal_posi		    = zeros(tide,1);            % PROFILE POSITION AT TIDAL RANGE
                            posi                = zeros(y_maxcell,1)        % PROFILE POSITION OF WHOLE VERTICAL RANGE
                            save_profile	    = zeros(maxit,y_maxcell);	% SAVE PROFILE ITERATIVELY
                    
                            count_wave_erosion      = zeros(maxit,1);           % COUNT THE WAVE EROSION FORCE IN EACH ITERATION
                            count_wea_erosion       = zeros(maxit,1);           % COUNT THE WEATHERING IN EACH ITERATION
                            count_wave              = zeros(maxit,1);           % COUNT THE NUMBER OF BLOCKS Y WAVE EROSION
                            count_wave_int          = zeros(maxit,1);           % COUNT THE NUMBER OF INTERTIDAL BLOCKS ERODED BY WAVE
                            count_wea               = zeros(maxit,1);           % COUNT THE NUMBER OF BLOCKS ERODED BY WEATHERING
                            count_vert_wave         = zeros(maxit,y_maxcell);   % STORE WAVE EROSION FORCE IN EACH ITERATION
                            count_vert_wea          = zeros(maxit,y_maxcell);   % STORE THE WEATHERING EROSION IN EACH ITERATION
                            
                            
                            %%% RESET OR CLEAR VARIABLES
                            ox_tidal_posi = 0;                              % OUTER PROFILE X POSITION WITHIN TIDAL RANGE
                            ix_tidal_posi = 0;                              % INNER PROFILE X POSITION WITHIN TIDAL RANGE
                            oy_tidal_posi = mwl-ceil(tide/2);               % OUTER PROFILE Y POSITION WITHIN TIDAL RANGE
                            iy_tidal_posi = mwl+ceil(tide/2)-1;             % INNER PROFILE Y POSITION WITHIN TIDAL RANGE

                            ox_posi = 0;                                    % OUTER PROFILE X POSITION 
                            oy_posi	= 1;                                    % OUTER PROFILE Y POSITION 
                            iy_posi = mwl+ceil(tide/2)-1;                   % INNER PROFILE Y POSITION 
					
        					sw2_decay 		= zeros(tide, maxit);           % DECAY CONSTAT FOR BROKEN WAVE
                			surf_width		= zeros(tide, maxit);           % SURF WIDTH
                        	angle_l		= zeros(tide, maxit);               % LOCAL ANGLE
        					angle_g		= zeros(maxit,1);                   % GLOBAL ANGLE
                            angle_g_tidal 	= zeros(maxit,1);               % GLOBAL INTERTIDAL ANGLE
					
                        	break_point_x   = zeros(tide,maxit);            % X BREAK POINT at TIDAL RANGE
                            break_point_y   = zeros(tide,maxit);            % Y BREAK POINT at TIDAL RANGE
                            tidal_x_posi    = zeros(tide,maxit);            % X POSITION AT INTERTIDAL ELEVATION
                            wave_type       = zeros(tide,maxit);            % WAVE TYPE at TIDAL RANGE
                            y_depth         = zeros(x_maxcell,1);				
			
                    
                            %%%  
                            %%% variables for multiple erosion
                            %%%
                            bw_ero_flag     = ones(maxit,1);  % Flag for bw multiple erosion
                            dw_ero_flag     = ones(maxit,1);  % Flag for dw multiple erosion
                            nex_bw_force  = zeros(y_maxcell,maxit);
                            nex_dw_force  = zeros(1,x_maxcell);
                    
                    
                            %%% RESET FOR NON WEATHERING CASE
                            if (wea==2) 
                                weathering = true;
                                wea_flag = 1;
                            end
                    
                            
                            %%% RESET FOR NON WEATHERING CASE
                            failure_flag = 0;
                 
                            mwl = reset_mwl;    % STORE INITIAL MWL
                            tec_flag = true;    % RESET TEC_FLAG
                            rsl_con = 1;        % RSL COUNTER
                 
                 
                 
                 
                    
                            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
                            %%% MAIN EROSION LOOP  
                            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
                            
                            
                             

                            loop = 0;
                            while (loop<maxit) & ((ix_max<(x_maxcell-(1/gridsize))) |(ix_posi<(x_maxcell-(1/gridsize)))|(ix_tidal_posi<(x_maxcell-(1/gridsize))))                
                                loop=loop+1;
                        
                        
                                %%% 
                                %%% RSL CHANGE DUE TO HOLOCENE SEA LEVEL FLUCTUATION
                                %%%                       
                                if(rsl_flag)
                                    if(loop==event(rsl_con))
                                        ver_rsl = rsl(rsl_con+1) - rsl(rsl_con);
                                        hor_rsl = event(rsl_con+1) - event(rsl_con);
                                        mwl = mwl + round(((ver_rsl/hor_rsl)*(loop-event(rsl_con)))/gridsize);
                                        rsl_con = rsl_con + 1;
                                    end
                                end
                        
                                %%% 
                                %%% RSL CHANGE DUE TO TECTONIC EVENTS
                                %%%
                                if(tec_flag)
                                    for i=1:length(ttic)
                                        if( tec_flag & loop==ttic(i) )
                                            mwl = mwl - tectonic(i,co)/gridsize; 
                                            break
                                            if(length(ttic)==i) tec_flag = false; end
                                        end	
                                    end
                                end
                        
                                
                                %%%
                                %%% 
                                %%% CALCULATE BACK WEARING EROSION FORCE
                                %%% 
                                %%%
                                if ( bw )
							
                                    %%% CLEAR PREVIOUS VALUES
                                    i=0; w_type=0; bw_erosion(:) = 0;   
                            
                                    %%% FOR ALL INTERTIDAL ELEVATION
                                    for y=mwl-ceil(tide/2):mwl+ceil(tide/2)-1
                                        i=i+1;
								
                                        
                                        %%% ESF SIZE CHECK
                                        if( i > size(esf,1))
                                            break
                                        end

                                        
                                        %%% FIND X POSITIONS OF SURFICIAL ROCKS
						                x=1;
                                        while (x<=x_maxcell)&(s.flag(x,y)==0)   x = x+1;    end
								
                                        
                                        %%% ESTIMATE HORIZONTAL BREAKING POINT
                                        %%% bp_x ... X POSITION OF BREAKING POINT
                                        y_bd = 0;
                                        if (x==1) bp_x=x; bp_y=0;
                                        else			
                                            bp_x=x;
                                            for xx=x-1:-1:1
                                                for yy=(y-bd):y-1
                                                    if ( s.flag(xx,yy+1)==0 & s.flag(xx,yy)==1)
                                                        bp_y=yy;
                                                        if ( bp_x > xx ) bp_x = xx; end
                                                    end
                                                end
                                            end
                                        end
                                
                                        
                                        %%% SET WAVE TYPE... UNBROKEN? BREAKING? BROKEN?
                                        if (x==1)                                               w_type = 1; %%% UNBROKEN WAVE
                                        elseif ( (x-bp_x)<0 )                                   w_type = 1; %%% UNBROKEN WAVE    
                                        elseif ( ((x-bp_x)<break_wave_dist)&((x-bp_x)>=0) )     w_type = 2; %%% BREAKING WAVE at bp_x
                                        elseif ( (x-bp_x)>=break_wave_dist )                    w_type = 3; %%% BROKEN WAVE at bp_x
                                        end
                                   
                                        
                                        %%% FORCE WEATHEIRNG INACTIVE
                                        if ( wea==2 & wea_flag==1 & w_type==2 )
                                    	    weathering  = false;
                                    	    wea_flag    = 0;
                                        end
                                    
                                        
                                        %%% SET SLOPE ANGLE BETA TO DETERMINE SURF WIDTH
                                        if ( loop==1 | w_type == 1 ) 
                                            angle = i_angle;
                                        elseif ( w_type==2 | w_type==3 )
                                            if ((x-bp_x+1) >= 2)
                                                clear tmp_y_posi
                                                tmp_x_posi = (1:(x-bp_x+1));
                                                tmp_y_posi = zeros(1,(x-bp_x+1));
                                                tmp_y_posi(1)		= bp_y;
                                                tmp_y_posi(x-bp_x+1)= y;
                                                
                                                tmp_i = 1;
                                                for tmp_x=bp_x+1:x-1
                                                    tmp_y = bp_y; 
                                                    tmp_i = tmp_i + 1;
                                                    while (s.flag(tmp_x,tmp_y)==1) tmp_y=tmp_y+1; end
                                                    tmp_y_posi(tmp_i)=tmp_y-1;
                                                end
                                            	pol=polyfit(tmp_x_posi, tmp_y_posi, 1);
                                                angle= atan(pol(1))*180/pi;    
                                            else
                                                angle = i_angle;
                                            end
                                        end
                                    
                                        
                                        %%% SET SURF WIDTH (cell)
                                        %%% SET WAVE ATTENUATION CONSTANT2
                                        if ( angle < 90 )
                                            sw = abs(num_height/tan(angle*pi/180));
                                            decay_sw2 = -log(swd2)/sw;
                                            %decay_sw2 = -log(1e-10)/sw;
                                        else
                                            sw = 0;
                                            decay_sw2 = Inf;
                                        end
									
                                        
                                        %%% PRIMARY EROSION
                                        %%% IN CASE OF STANDING WAVE
                                        if ( w_type == 1 )
                                            k=0;
                                            for yy=y+(stan_dist1-1):-1:y-stan_dist2
                                                k=k+1;
                                                xx=1;
                                                while (xx<=(x_maxcell-1)) & (s.flag(xx,yy)==0)    xx=xx+1;    end
                                                force = stan_const*height*esf(i)*p_stan(k); 
                                                bw_erosion(yy) = bw_erosion(yy) + force;
                                            end
									
                                        %%% IN CASE OF BREAKING WAVE
                                        elseif ( w_type == 2 )			
                                            k=0;    k1=-((break_dist1-1)-(stan_dist1-1));  k2=-((break_dist1-1)-(brok_dist1-1));
                                            
                                            % this loop is the same as the one above because jump=0
                                            for yy=round(jump*(height/gridsize))+y+(break_dist1-1):-1:round(jump*(height/gridsize))+y-break_dist2
                                                k=k+1;  k1=k1+1;    k2=k2+1;

                                                xx=1;
                                                while (xx<=(x_maxcell-1))&(s.flag(xx,yy)==0) xx=xx+1; end
                                        
                                                %%% SPATIALLY BEHAVE AS A STANDING WAVE
                                                if ( xx < bp_x ) & ( (stan_dist1-1) >= (yy-y) ) & ( stan_dist2 >= (y-yy) )                                                    
                                                    force = stan_const*height*esf(i)*p_stan(k1);
                                                    bw_erosion(yy) = bw_erosion(yy) + force;
                                        	  
                                                %%% SPATIALLY BEHAVE AS A BROKEN WAVE                                              
                                                elseif ( ( (xx-bp_x)>= break_wave_dist ) & ( (brok_dist1-1) >= (yy-y) ) & ( brok_dist2 >= (y-yy)) ) 
                                                    force = brok_const*height*bwd*esf(i)*p_brok(k2)*exp(-decay_sw2*(xx-(bp_x+break_wave_dist)));
                                                    bw_erosion(yy) = bw_erosion(yy) + force;
                            				    
                                                %%% OR...    
        										elseif (( xx >= bp_x ) & ( (xx-bp_x) < break_wave_dist ) & ( break_dist1-1 >= (yy-y) ) & ( break_dist2 >= (y-yy) ))
                                                    force = brea_const*height*esf(i)*p_break(k)*exp(-decay_sw1*(xx-(bp_x)));
                                                   	bw_erosion(yy) = bw_erosion(yy) + force;
                                            
                                                end    
                                            end
									
                                        %%% IN CASE OF BROKEN WAVE
                                        elseif ( w_type == 3 )
										
                                            k=0;    
                                            k1=-((brok_dist1-1)-(stan_dist1-1));    
                                            k2=-((brok_dist1-1)-(break_dist1-1));
										
                                            for yy=y+brok_dist1-1:-1:y-brok_dist2
                                                k=k+1;     k1=k1+1;    k2=k2+1;
                                                
                                                xx=1;   
                                                while (xx<=(x_maxcell-1)) & (s.flag(xx,yy)==0)    xx=xx+1;    end
                                    
                                                %%% SPATIALLY BEHAVE AS A STANDING WAVE
                                                if ( xx < bp_x ) & ( (stan_dist1(tmp_p)-1) >= yy-y ) & ( stan_dist2(tmp_p) > y-yy )       
                                                    force = stan_const*height*esf(i)*p_stan(k1);
        											bw_erosion(yy) = bw_erosion(yy) + force;
												
                        						%%% SPATIALLY BEHAVE AS A BREAKING WAVE
                                                elseif ( (xx>=bp_x) & ((xx-bp_x)<break_wave_dist) & ( (break_dist1(tmp_p)-1) >= (yy-y) ) & ( break_dist2(tmp_p) >= (y-yy) ) ) 
                                        
                                                    force  = brea_const*height*esf(i)*p_break(k2)*exp(-decay_sw1*(xx-(bp_x)));
                                                	bw_erosion(yy) = bw_erosion(yy) + force;
                                        		
                                                %%% OR...
                                				elseif ( ((xx-bp_x) >= break_wave_dist) & ( (brok_dist1(tmp_p)-1) >= (yy-y) ) & ( brok_dist2(tmp_p) >= (y-yy)) )
                                        			
                                                    force = brok_const*height*esf(i)*p_brok(k)*exp(-decay_sw2*(xx-(bp_x+break_wave_dist)));
                                            		bw_erosion(yy) = bw_erosion(yy) + force;
                                                end
                                            end
                                        end
									
                                        %%% STORE WAVE TYPE & BREAK POINT
                                        break_point_x(i,loop)   = bp_x;
                                        break_point_y(i,loop)   = bp_y; 
        								wave_type(i,loop)       = w_type;
                						surf_width(i,loop)      = sw;
                        				sw2_decay(i,loop)       = decay_sw2;
                                        angle_l(i,loop)			= angle;
                                        tidal_x_posi(i,loop)    = x;
								
        								if(force==inf)  stop; end
                                
                                    end
                                    
                                    %%% IN CASE OF ONE BLOCK EROSION ... CAN IGNORE
                                    if( loop==1 )
                                        total_bw_erosion(:,loop)  	= bw_erosion;
                                    else
                                        if( max_ero_flag )
                                            if( bw_ero_flag(loop-1) )   total_bw_erosion(:,loop)  	= bw_erosion;
                                            else                        total_bw_erosion(:,loop)  	= nex_bw_force(:,loop-1);
                                            end
                                        else
                                            total_bw_erosion(:,loop)  	= bw_erosion;
                                        end
                                    end
                                    
                                end

                                
                                %%%
                                %%%
                                %%% CALCULATE DOWN WEARING EROSION FORCE
                                %%%
                                %%%
                                if (dw)
							
                                    %%% CLEAR PREVIOUS VALUE
                                    dw_erosion(:) = 0;
							
                                    %%% ESTIMATE SURFICIAL ROCKS WITHIN INTERTIDAL RANGE
                                    i=0; x_tidal_max = 0;   y_tidal_max = 0;
                                    for y=mwl-ceil(tide/2):(mwl+ceil(tide/2)-1)
                                        x=1; 
                                        i=i+1;
                                        while( s.flag(x,y)==0 & x<=ix_max-1 )	x=x+1;	end
                                        tidal_posi(i) = x;
                                        if ( x > x_tidal_max )
                                            x_tidal_max = x;
                                            y_tidal_max = y;
                                        end
                                    end
						
        							%%% ESTIMATE VERTICAL POSITION OF SUBMARINE SURFICIAL ROCKS
                					if ( x_tidal_max >= 1 )	
                        				for x=1:x_tidal_max-1
                                			y = y_tidal_max;
                                        	while (s.flag(x,y)==0) & (y > 1)  y=y-1;    end
        									if (y > 1)      y_depth(x) = y;	
                                            else            y_depth(x)=1;
                                            end
                                        end
								
                                        i = 0;
                                        for y=mwl-ceil(tide/2):mwl+ceil(tide/2)-1
                                            i = i + 1;
                                    
                                            %%% ESF SIZE CHECK
                                            if( i > size(esf,1))
                                                break
                                            end
                                    
                                            bp_x = break_point_x(i,loop);   % OBTAIN BREAKING POINT
                							decay_sw2 = sw2_decay(i,loop);  % BROKEN WAVE DECAY COEFFICIENT
                                            
                            				if ( bp_x >= 1 )
                                                
                                                x=1;
                                                while ( x < tidal_posi(i) )
                                                    %%% IN CASE OF STANDING WAVE
                                                    if (x<bp_x)	
                                                        force = stan_const*height*esf(i)*p_stan_sub;
                                                        decay_dc = -log(sub_decay_const)/(height);
                                                    %%% IN CASE OF BREAKING WAVE
                                                    elseif( x>=bp_x & x<(bp_x+break_wave_dist) ) 
                                                        force = brea_const*height*esf(i)*p_break_sub*exp(-decay_sw1*(x-(bp_x)));
                                                        decay_dc = -log(sub_decay_const)/(height*exp(-decay_sw1*(x-(bp_x))));
                                                    %%% IN CASE OF BROKEN WAVE
                                                    else 
                                                        force = brok_const*bwd*esf(i)*p_brok_sub*exp(-decay_sw2*(x-(bp_x+break_wave_dist)));
                                                        decay_dc = -log(sub_decay_const)/(height*bwd*exp(-decay_sw2*(x-(bp_x+break_wave_dist))));
                                                    end

                                                    %%% ESTIMATE FORCE SURFICIAL ROCKS
        											yy=y-y_depth(x);
                									if ( yy > 0 )
                                						dw_erosion(x) = dw_erosion(x) + force * exp(-decay_dc*yy);
                                        				if ( y_depth(x) >= (mwl-ceil(tide/2)) & (y_depth(x) < mwl+ceil(tide/2)-1) )
                                                			dur=y_depth(x)-(mwl-ceil(tide/2))+1;
                                                        end
                                                    end											
                                                    x=x+1;											
                                                end
                                            end									
                                        end
                                
                                    end
                                    
                                    %%% STORE FOR OUTPUTTING
                                    if(dw_ero_flag(loop))   total_dw_erosion(:,loop) = dw_erosion;
                                    else                    total_dw_erosion(:,loop) = nex_dw_force;
                                    end
                                end
						
                                
                            %%%
                        	%%%
    						%%% Judge erosion ( FW >= FR? )
        					%%%
                            %%%
                            
                            %%%
            				%%% BACK WEARING EROSION
                			%%%
                            if(bw)
                                
                                bw_ero_flag(loop) = 1;

                                for y=1:(mwl+ceil(tide/2)-1+break_dist1)
	
                                    x=1; 
                                    while ( (s.flag(x,y)==0) & ( x < (x_maxcell-10) ) ) x=x+1; end
							
                                    bw_force=total_bw_erosion(y,loop);
                                    if ( bw_force >= s.resistance(x,y) )
                                        tmp_resi = s.resistance(x,y);
                                        
                                        %%% IN CASE OF WAVES IS ABLE TO ACHIVE ONLY ONE ROCK CELL EROSION 
                                        if (change == 1) 
                                    
                                            s.flag(x,y)             = 0;	
                                            s.resistance(x,y)       = 0;
                                    		count_vert_wave(loop,y) = count_vert_wave(loop,y) + tmp_resi;
                                            count_wave(loop)        = count_wave(loop) + 1; 
                                    
                                            if ( y >= mwl-ceil(tide/2) & y <= mwl+ceil(tide/2)-1 )
                                                count_wave_erosion(loop)    = count_wave_erosion(loop)    + tmp_resi;
                                                count_wave_int(loop)        = count_wave_int(loop) + 1;
                                            end
                                            if ( ix_max < (x+1) )  ix_max = (x+1); end
                                            if ( ix_max >= x_maxcell ) ix_max = x_maxcell-1; end
								
                                        %%% IN CASE OF WAVE IS ABLE TO ACHIVE MORE THAN ONE ROCK CELL EROSION 
                                        elseif (change ==2)

                                            %%% CAN ERODE MORE THAN ONE CELL!!
                                            if (bw_force >= (s.resistance(x,y)+effec_resi))
                                                tmp_ero = bw_force-s.resistance(x,y); 
                                                q       = floor(tmp_ero/effec_resi); 
                                                r       = mod(tmp_ero,effec_resi);

                                                ero_count   = 0;
                                                for p=0:q
                                                    if ((x+p)<=x_maxcell)
                                                        s.flag(x+p,y)=0; 
                                                        s.resistance(x+p,y)=0;
                                                        count_vert_wave(loop,y) = count_vert_wave(loop,y) + effec_resi;
                                                        count_wave(loop)        = count_wave(loop) + 1;
                                                    
                                                        if ( y >= mwl-ceil(tide/2) & y <= mwl+ceil(tide/2)-1 )
                                                            count_wave_erosion(loop)    = count_wave_erosion(loop)  + effec_resi;
                                                            count_wave_int(loop)        = count_wave_int(loop) + 1;
                                                        end
                                                        if ( ix_max < (x+p+1) )     ix_max = (x+p+1); end
                                                        if ( ix_max >= x_maxcell )  ix_max = x_maxcell-1; end
                                                    end
                                                end
                                        
                                            %%% CAN ERODE ONLY ONE CELL    
                                            else
                                                s.flag(x,y)             = 0;	
                                                s.resistance(x,y)       = 0;
                                                count_vert_wave(loop,y) = count_vert_wave(loop,y) + tmp_resi;
                                                count_wave(loop)        = count_wave(loop) + 1;
                                        
                                                if ( (y >= mwl-ceil(tide/2)) & (y <= mwl+ceil(tide/2)-1) )
                                                    count_wave_erosion(loop)    = count_wave_erosion(loop) + tmp_resi;
                                                    count_wave_int(loop)        = count_wave_int(loop) + 1;
                                                end
                                                if ( ix_max < (x+1) )  ix_max = (x+1); end
                                                if ( ix_max >= x_maxcell ) ix_max = x_maxcell-1; end
                                        
                                            end
                                        end
                                    end
                                end		
                            end     %%% END OF BACK WEARING EROSION
                        

                            %%% 
                            %%% DOWN WEARING  EROSION
                            %%%
                            if(dw)
                                dw_ero_flag(loop) = 1;
                        
                                %%% FOR ALL THE SURFICIAL ROCK CELLS ...
                                for x=1:x_tidal_max-1
                                    y = y_depth(x);	
                                    dw_force = total_dw_erosion(x,loop);
							
                                    if (s.flag(x,y)==1)
							
                                        while ((s.flag(x,y)==0)&y>1) y=y-1; end 
								
                                        if ( dw_force >= s.resistance(x,y) )
                                            tmp_resi = s.resistance(x,y);
									
                                            %%% IN CASE OF WAVES ACHIVING ONE ROCK CELL EROSION 
                                            if (change==1)
                                                s.flag(x,y)=0;	s.resistance(x,y)=0;
                                                count_vert_wave(loop,y) = count_vert_wave(loop,y) + tmp_resi;
                                                count_wave(loop)        = count_wave(loop) + 1;
                                        
                                                if ( (y >= mwl-ceil(tide/2)) & (y <= mwl+ceil(tide/2)-1) )
                                                    dur=y-(mwl-ceil(tide/2))+1;
                                                    count_wave_erosion(loop)    = count_wave_erosion(loop)  + tmp_resi;
                                                    count_wave_int(loop)        = count_wave_int(loop) + 1;
                                                end
									
                                            %%% IN CASE OF WAVE IS ABLE TO ACHIVE MORE THAN ONE ROCK CELL EROSION 
                                            elseif (change==2)
                                        
                                                %%% CAN ERODE MORE THAN ONE CELL!!
                                                if ( dw_force >= s.resistance(x,y)+effec_resi )
                                                    tmp_ero=dw_force-s.resistance(x,y); 
                                                    q=floor(tmp_ero/effec_resi); 
                                                    r=mod(tmp_ero,effec_resi);
										
                                                    ero_count   = 0;
                                                    for p=0:q
                                                        if ( (y-p) >= 1 )
                                                            s.flag(x,y-p)=0; s.resistance(x,y-p)=0;
                                                            count_vert_wave(loop,y-p)   = count_vert_wave(loop,y-p) + effec_resi;
                                                            count_wave(loop)            = count_wave(loop) + 1;
													
                                                            if ( ( (y-p) >= mwl-ceil(tide/2) ) & ( ( (y-p) <= mwl+ceil(tide/2)-1 ) ) )
                                                                count_wave_erosion(loop)    = count_wave_erosion(loop)  + effec_resi;
                                                                count_wave_int(loop)        = count_wave_int(loop) + 1;
                                                            end
                                                        end
                                                    end
                                                    
                                                %%% CAN ERODE ONLY ONE CELL    
                                                else
                                                    s.flag(x,y)=0;	s.resistance(x,y)=0;
                                                    count_vert_wave(loop,y) = count_vert_wave(loop,y) + tmp_resi;
                                                    count_wave(loop)            = count_wave(loop) + 1;
    										
                                                    if ( y >= mwl-ceil(tide/2) & (y <= mwl+ceil(tide/2)-1) )
                                                        count_wave_erosion(loop)    = count_wave_erosion(loop)  + tmp_resi;
                                                        count_wave_int(loop)        = count_wave_int(loop) + 1;
                                                    end
                                            
                                                end
                                            end				
                                        end
                                    end
                				end		
                            end    %%% END OF BACK WEARING EROSION
              
                            
			  
                            %%%
                            %%% CALCULATE MASS FAILURE : version 2
                            %%% 
                            %%% ORIGINALLY A MINIMUM CEL NUMBER WITH WICHIH CANTILEVER MASS COLLAPSE CALCULATED, 
                            %%% BUT NOW ALL THE CANTILEVER BLOCK IS REMOVED WHEN EROSION DEPTH IS 1 M.   
                            %%%  
                            x_dist=0;   y=mwl-ceil(tide/2)-break_dist1-1;
                            while (y <= (mwl+ceil(tide/2)+break_dist1) & y<=y_maxcell-2 )
                                y=y+1; 
							
                                x=1;
            					while (s.flag(x,y)==0) & (x <= (ix_max-1)) x=x+1;  end
        
                                %%% COUNT THE WIGHT and ARM LENGTH
                                count=0; x_dist=0;
                                if (x>1+round(y/tan(i_angle*pi/180)))
                                    if (s.flag(x-1,y+1)==1) 
                                        yy=y+1; x_dist=1;
                                        while (yy<=(y_maxcell-1))&(s.flag(x-1,yy)==1)
        									yy=yy+1;
            								if(yy == mwl+ceil(tide/2)+(break_dist1(tmp_p)-1)+1) count = count+600;
                							else	count=count+1;
                                            end
                                        end
                                        in=1; xx=x-2;
                                        while (xx>0) & (in==1)
                                            yyy=y+1;
                                            count_old=count;
                                            while (yyy <= yy)
                                                if ((s.flag(xx,yyy)==1)&(s.flag(xx+1,yyy)==1))
                                                    if (yy == mwl+ceil(tide/2)+(break_dist1(tmp_p)-1)+1) count=count+600;
    												else count=count+1;
                                                    end
                                                end
                                                yyy=yyy+1;
                                            end
                                            if (count_old == count) in=0;
                                            else    xx=xx-1; x_dist=x_dist+1;
                                            end
                                        end
                                    end
                                end
                                
                                %%% ORIGINAL CRITERIA
                                %if ( count*x_dist > (y_maxcell*y_maxcell)/10 )
                                %%% NOW CHANGES TO ... IF x_dist IS LARGER THAN 1 m... 
                                if ( x_dist > 10 )
                                    for xxx=x-1:-1:(x-x_dist)
                                        for yyyy=y+1:yy s.flag(xxx,yyyy)=0; end
                                    end
                                    failure_flag = 1;
                                end
                            end
                            
                            
                            %%%
                            %%% CALCULATE MASS FAILURE : version 2
                            %%%
                            if(false)
                            clear posi
                            posi = zeros(y_maxcell,1);	% set size of posi
                            f=0;
                            tidal_posi=0;
                            for y=1:y_maxcell
                                x=1;
                                while( s.flag(x,y)==0 & x <= ix_max )	x=x+1;	end
                                posi(y)=x;
                                if( (y>=(mwl-ceil(tide/2))+add_tide) & (y<=(mwl+ceil(tide/2)-1)) )
                                    f=f+1;
                                    tidal_posi(f)=x;
                                end
                            end
                            max_x_value = max(posi((mwl-ceil(tide/2)+add_tide):(mwl+ceil(tide/2)-1)));
                            max_y_ele   = min(find(max_x_value==posi));
                            for y=max_y_ele:y_maxcell
                                if( posi(y) < max_x_value )
                                    s.flag(1:max_x_value,y) = 0;
                                    s.resistance(1:max_x_value,y) = 0;
                                end
                            end
                            end
                           
                            
                            %%%
                            %%% FIND INNER and OUTER POSITIONS
                            %%%
                            
    						%%% find inner and outer position below MLWS
        					ox_posi	= min(posi);		
            				ix_posi = max(posi);
                			oy_temp = find( posi==ox_posi );
                    		iy_temp = find( posi==ix_posi );
                            oy_posi = max(oy_temp);
                        	iy_posi = min(iy_temp);
                        	g_angle = atan ((iy_posi-oy_posi)/(ix_posi-ox_posi))*(180/pi);
                            if ( ox_posi==ix_posi ) 
                                g_angle =90;
                                oy_posi = iy_posi;
                            end
							
                            % find inner position within tidal range
                            ix_tidal_posi 	= max(tidal_posi);
                            ox_tidal_posi	= min(tidal_posi);
                            iy_tidal_temp 	= find( tidal_posi==ix_tidal_posi );
                            oy_tidal_temp 	= find( tidal_posi==ox_tidal_posi );
                            iy_tidal_posi 	= mwl-(ceil(tide/2))-1+min(iy_tidal_temp);
                            oy_tidal_posi 	= mwl-(ceil(tide/2))-1+max(oy_tidal_temp);
                            if ( ox_tidal_posi==ix_tidal_posi ) 
                                oy_tidal_posi = mwl-ceil(tide/2);
                                iy_tidal_posi = mwl+ceil(tide/2)-1;
                            end
                               
                            
						
                            %%%
                            %%%
                            %%% WEATHERING PROCESS
                            %%%					
                            %%%
                            if (weathering)
                				i=0;

                    			% clear previous value
                        		wea_erosion(:) = 0;
                                for y=(mwl+ceil(tide/2)-1):-1:(mwl-ceil(tide/2))
                                    i=i+1;
                                
                                    %%% ESF SIZE CHECK
                                    if( i > size(esf,1))    break;  end
                    
                                    %force = effec_resi*wt_const*wetdry(i);
                                    previous = 0;
                                    for x=1+round(y/tan(i_angle*pi/180)):ix_tidal_posi %%% CHECK ALL THE SURFICIAL ROCKS
                                    
                                        %%% WEATHEIRNG FORCE
                                        force = wt_const*(gridsize/0.1)*wetdry(i);
                                    
                                        %%% FIND ROCK CELL
                                        if( s.flag(x,y)==1 )
                                        
                                            %%% IF INITIAL POSITION, NO NEED TO CHECK EXPOSURE
                                            if ( x==1+round(y/tan(i_angle*pi/180)) ) 
                                                for yy=0:max_we_erosion-1
                                                
                                                    if( y-yy <= mwl+ceil(tide/2)-1 & y-yy >= (mwl-ceil(tide/2)) )
                                                        left = s.resistance(x,y-yy);
                                                        s.resistance(x,y-yy)        = s.resistance(x,y-yy)      - force;
                                                        wea_erosion(y-yy)           = force;
                                                        count_wea_erosion(loop)     = count_wea_erosion(loop)   + force;
                                                        count_vert_wea(loop,y-yy)   = count_vert_wea(loop,y-yy) + force;

                                                        if( s.resistance(x,y-yy) < 0 )
                                                            previous = x;
                                                            s.flag(x,y-yy)              = 0;	
                                                            s.resistance(x,y-yy)        = 0;
                                                            count_wea_erosion(loop)     = count_wea_erosion(loop)   - (force-left);
                                                            count_vert_wea(loop,y-yy)   = count_vert_wea(loop,y-yy) - (force-left);
                                                            count_wea(loop)             = count_wea(loop) + 1;
                                                            force                       = force-left;
                                                        
                                                            if ( ix_max < x )           ix_max          = x; end
                                                            if ( ix_max >= x_maxcell )  ix_max          = x_maxcell-1; end
                                                            if ( iy_tidal_posi < y-yy ) iy_tidal_posi   = y-yy; end
                                                        end
                                                    end
                                                end
                                            
                                            %%% IF X_MAXCELL POSITION
                                            elseif ( x == x_maxcell )
                                            
                                                for yy=0:max_we_erosion-1
                                                    if( y-yy <= mwl+ceil(tide/2)-1 & y-yy >= (mwl-ceil(tide/2)) )
                                                        if ( (s.flag(x-1,y-yy)==0) | (s.flag(x,y+1-yy)==0) | (s.flag(x,y-1-yy)==0) )  
                
                                                            left = s.resistance(x,y-yy);
                                                            s.resistance(x,y-yy)        = s.resistance(x,y-yy)      - force;
                                                            wea_erosion(y-yy)           = force;
                                                            count_wea_erosion(loop)     = count_wea_erosion(loop)   + force;
                                                            count_vert_wea(loop,y-yy)   = count_vert_wea(loop,y-yy) + force;
                               			
                                                            if( s.resistance(x,y-yy) < 0 )
                                                    
                                                                previous = x;
                                                                s.flag(x,y-yy)              = 0;	
                                                                s.resistance(x,y-yy)        = 0;
                                                                count_wea_erosion(loop)     = count_wea_erosion(loop)   - (force-left);
                                                                count_vert_wea(loop,y-yy)   = count_vert_wea(loop,y)     - (force-left);
                                                                count_wea(loop)             = count_wea(loop) + 1;
                                                                force                       = force - left;
                                                       
                                                                if ( ix_max < x )           ix_max = x; end
                                                                if ( ix_max >= x_maxcell )  ix_max = x_maxcell-1; end
                                                                if ( iy_tidal_posi < y-yy ) iy_tidal_posi = y-yy;	end
                                                        
                                                            end
                                                        end
                                                    end
                                                end
                                       
                                            %%% OTHER ....    
                                            else
        
                                                for yy=0:max_we_erosion-1
        
                                                    if( ((s.flag(x-1,y-yy)==0)&(previous~=(x-1))) | (s.flag(x,y+1-yy)==0) | (s.flag(x,y-1-yy)==0) | (s.flag(x+1,y-yy)==0))
                    
                                                        if( y-yy <= mwl+ceil(tide/2)-1 & y-yy >= (mwl-ceil(tide/2)) ) 
                                                    
                                                            left = s.resistance(x,y-yy);
                                                            s.resistance(x,y-yy)        = s.resistance(x,y-yy)      - force;
                                                            wea_erosion(y-yy)           = force;
                                                            count_wea_erosion(loop)     = count_wea_erosion(loop)   + force;
                                                            count_vert_wea(loop,y-yy)   = count_vert_wea(loop,y-yy) + force;
                                                
                                                            if( s.resistance(x,y-yy) < 0 )
                                                                previous                    = x;
                                                                s.flag(x,y-yy)              = 0;	
                                                                s.resistance(x,y-yy)        = 0;
                                                                count_wea_erosion(loop)     = count_wea_erosion(loop)   - (force-left);
                                                                count_vert_wea(loop,y-yy)   = count_vert_wea(loop,y-yy) - (force-left);
                                                                count_wea(loop)             = count_wea(loop) + 1;
                                                                force                       = force - left;
                                                            
                                                                if ( ix_max < x )           ix_max = x; end
                                                                if ( ix_max >= x_maxcell )  ix_max = x_maxcell-1; end
                                                                if ( iy_tidal_posi < y-yy ) iy_tidal_posi = y-yy; end
                                                            end
                                                        end
                                                    end
                                                end
                                            end
                                        end
                                    end
                                end
							
                                total_wea_erosion(:,loop)  = wea_erosion(:);
                            end

                            
                            
                            %%%
                            %%% SUPRATIDAL CLIFF WEATHERING << NOT USED NOW
                            %%%
                            if( spr_weathering & count_wea_erosion(loop)>count_wave_erosion(loop) )
                                elevation = 0;
    							for y=(mwl+ceil(tide/2)):y_maxcell;
        							elevation=elevation+1;
            						if(i_angle==90) x=1;
                					else			x=1+round(y/tan(i_angle*pi/180));
                                    end
                                    while( s.flag(x,y)==0 & x < x_maxcell )	x=x+1;	end
                            		cl_posi(elevation) = x;
                                end
                                cl_ang = atan ((cl_posi(elevation)-cl_posi(1))/elevation)*(180/pi);
                                elevation = 0;
    							for y=(mwl+ceil(tide/2)):y_maxcell;
        							elevation=elevation+1;
            						if(cl_ang+0.5<60)
                						recess = round(tan((cl_ang+0.5)*pi/180)*elevation);
                    					cl_posi(elevation) = cl_posi(1)+recess;
                        				if(i_angle==90) x=1;
                            			else			x=1+round(y/tan(i_angle*pi/180));
                                        end
                                    	while ( x <= cl_posi(elevation)) 
                                        	s.flag(x,y)=0;
                                            x=x+1;
                                        end
                                    end
                                end
                            end
                        
                        
                        
                            %%%    
                            %%% UPDATE X_&Y_PROFILE, INNER_POSISION
                            %%% CALCULATE FINAL ANGLE
                            %%%
                            
                            clear posi                  % as size change due to sea level change
                            posi = zeros(y_maxcell,1);	% set size of posi
                            tidal_posi=0; 
                            f=0;
						
                            %%% FIND INNER POSITION 
                            for y=1:y_maxcell
							    x=1+round(y/tan(i_angle*pi/180));
                                while( s.flag(x,y)==0 & x <= ix_max-1 )	x=x+1;	end
                                posi(y)=x;
                                if( (y>=(mwl-ceil(tide/2)+add_tide)) & (y<=(mwl+ceil(tide/2)-1)) )
                                    f=f+1;
                                    tidal_posi(f)=x;
                                end
                            end
                        
                            
                            
                            %%%
                            %%% ANGLE CALCULATION
                            %%%
                            g_angle = atan ((iy_posi-oy_posi)/(ix_posi-ox_posi))*(180/pi);
                            
                            if ( ox_posi==ix_posi ) 
                                g_angle =90;
                                oy_posi = iy_posi;
                            end
                                                
                            if (failure_flag == 0)  g_tidal_angle = 90;
                            else
                                if(wea == 1)
                                    if (ix_tidal_posi==ox_tidal_posi)
                                        oy_tidal_posi = iy_tidal_posi;
                                        g_tidal_angle = 90;
                                    elseif (iy_tidal_posi == (mwl-ceil(tide/2)))
                                        ox_tidal_posi = ix_tidal_posi;
                                        oy_tidal_posi = iy_tidal_posi;
                                        g_tidal_angle = 90;
                                    else
                                        g_tidal_angle   = atan((iy_tidal_posi-(mwl-ceil(tide/2)))/(ix_tidal_posi-tidal_posi(1)))*(180/pi);
                                    end
                                
                                elseif(wea == 2)
                                    oy_posi = 1;
                                    while (s.flag(1,oy_posi)==1) oy_posi = oy_posi + 1; end
                                    g_tidal_angle = atan((iy_posi-oy_posi)/(ix_posi))*(180/pi);
                                end
                            end
                        
                        
                        
                            %%%    
                            %%% PLOT FIGURE
                            %%%
                            if (print_snap_to_mov)
                            
                                if ( (loop==1) | ((loop<=10) & (mod(loop,print_loop1)==0)) | ((loop<=500) & (mod(loop,print_loop2)==0)) | ((loop<=ttic) & (mod(loop,print_loop3)==0)) | ...
                                    ( (loop<=ttic+10) & (loop>ttic) & (mod(loop,print_loop1)==0)) | ((loop<=ttic+500) & (loop>ttic+10) & (mod(loop,print_loop2)==0)) | ((loop>ttic+500)&(mod(loop,print_loop3)==0)) )
                                
                                    if ( ix_max+10 > x_maxcell )	ix_max = x_maxcell-10;	end
                                    x_temp1	=1:ix_max+10;
                                    x_temp2	=1:ix_max+10;
                                    x_line	=1:ix_max+10;
                                    for i=2:y_maxcell     x_temp1=horzcat(x_temp1,x_temp2);   end
                            
                                    y_temp1	= s.flag(1:ix_max+10,1)';
                                    for i=2:y_maxcell     
                                        y_temp2	= i*s.flag(1:ix_max+10,i)';
                                        y_temp1	= horzcat(y_temp1,y_temp2);
                                    end
                
                                    y_high(1:size(x_temp2))=mwl+ceil(tide/2)-1;
                                    y_low(1:size(x_temp2))=mwl-ceil(tide/2);
                                    winsize = [ 0 0 1920 1080];

                                    figure('visible', 'off', 'Position', winsize); hold on;  axis([1 x_maxcell 1 y_maxcell]);
                                    plot(x_temp1,y_temp1, '.');
                                    plot(x_line, y_high, '--r', 'linewidth', 10);
                                    plot(x_line, y_low, '--r','linewidth',10);
                                    text(2,y_maxcell-(0*y_maxcell/10+1),['Iteration=',num2str(loop)],'Fontsize',10,'Backgroundcolor',[.7 .9 .7]);
                                    text(2,y_maxcell-(0.5*y_maxcell/10+1),['Rock Hardness=',material_resistance],'Fontsize',10,'Backgroundcolor',[.7 .9 .7]);
                                    text(2,y_maxcell-(1.5*y_maxcell/10+1),['Tidal angle=',num2str(g_tidal_angle)],'Fontsize',10,'Backgroundcolor',[.7 .9 .7]);
							
                                    %%% MAKE AVI MOVIE
                                    if (loop==1)
                                        extension = '.avi';
                                        filename = [output_path, profile, num2str(i_angle), connect, wave_height, connect, material_resistance, extension];
                                        aviobj=avifile(filename, 'fps', 4);
                                        aviobj=addframe(aviobj,gcf);
                                    else
                                        aviobj=addframe(aviobj,gcf);
                                    end
                                    close;
                                end
                            end
						
                            %%%
                            %%% STORE VALUABLES
                            %%%
                            angle_g(loop)			= g_angle;
                            angle_g_tidal(loop)		= g_tidal_angle;
                            	
                            
                            %%%
                            %%% SAVE PROFILE
                            %%%
                            for y=1:y_maxcell
                                x=1;
                                while ( s.flag(x,y)==0 & x < x_maxcell ) x=x+1; end
                                save_profile(loop,y) = x;
                            end
						
                            
                        end
                        
                        
                        
                        
                        
                        
                        
                        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                        %%% END OF EROSION MAIN LOOP
                        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                        
                        
                        
                        
                        
                        
                        
                        
                        
                        %%% CLOSE AVI FILE
                        if(print_snap_to_mov)
                            aviobj=close(aviobj);   
                        end
                    
                        %%% COPY TOTAL RESULT
        				final_prof(:,resi)	= save_profile(loop,:);
            			final_it(resi)		= loop;
                        store_mwl(co)       = mwl;
                                          
                    
                        
                        %%%
                        %%% PRINT FIGURES 
                        %%%           
                
                        %%% OUTPUT final_PROF to txt
                        if(print_bwer_to_txt)
                            extension = '-total_bwero.txt';
                            filename=[output_path, profile, num2str(i_angle), connect, wave_height, connect, material_resistance, extension];
                            txt = total_bw_erosion;
                            %t_txt = flipud(txt.');
                            dlmwrite(filename, txt);
                            clear txt
                        end
                
                        %%% OUTPUT final_PROF to txt
                        if(print_nxbw_to_txt)
                            extension = '-nex_bwero.txt';
                            filename=[output_path, profile, num2str(i_angle), connect, wave_height, connect, material_resistance, extension];
                            txt = nex_bw_force;
                            %t_txt = flipud(txt.');
                            dlmwrite(filename, txt);
                            clear txt
                        end
                    
                        %%% OUTPUT final_PROF to txt
                        if(final_prof_to_txt)
                            extension = '-final_prof.txt';
                            filename=[output_path, profile, num2str(i_angle), connect, wave_height, connect, material_resistance, extension];
                            txt = s.flag;
                            t_txt = flipud(txt.');
                            dlmwrite(filename, t_txt);
                            clear txt t_txt
                        end
                    
                        %%% OUTPUT break_point_x to txt
                        if(print_break_point)
                            extension = '-break_point_x.txt';
                            filename=[output_path, profile, num2str(i_angle), connect, wave_height, connect, material_resistance, extension];
                            dlmwrite(filename, break_point_x(:,:));
                        end                    
                    
                    
                        %%% OUTPUT OF WAVE vs. WEATHEIRNG TO FIG
                        if (loop ~= maxit)	
                            x=[1:loop];
                            x_mean=[5:loop-5];
                            count_wave_erosion(loop+1:maxit)= [];
                            count_wea_erosion(loop+1:maxit) = [];
                            count_vert_wave(loop+1:maxit,:) = [];
                            count_vert_wea(loop+1:maxit,:)  = [];
                            count_wave(loop+1:maxit)        = [];
                            count_wave_int(loop+1:maxit)    = [];
                            count_wea(loop+1:maxit)        = [];
                        
                        else
                            x=[1:maxit];
                            x_mean=[5:maxit-5];
                        end
                    
                    
                    
                        %%% OUTPUT OF WAVE vs. WEATHEIRNG TO TXT
                        if(print_wvsw_to_txt)
                    
                            %%% Calculate mean
                            clear mean_count_wave_erosion mean_count_wea_erosion 
                            mean_count_wave_erosion = zeros(loop-9,1);
                            mean_count_wea_erosion = zeros(loop-9,1);
                            for me=1:loop-9
                                mean_wave = mean(count_wave_erosion(me:me+8));
                                mean_wea  = mean(count_wea_erosion(me:me+8));
                            end
                    
                            extension = '-wave_vs_wea.txt';
                            filename=[output_path, profile, num2str(i_angle), connect, wave_height, connect, material_resistance, extension];
                            txt = x;
                            txt = vertcat(txt,count_wave_erosion');
                            txt = vertcat(txt,count_wea_erosion');
                        
                            loop_tmp = loop;
                            if(loop<=8) 
                                loop_tmp=loop;
                                loop=100;
                            end
                        
                            if(~isempty(mean_count_wave_erosion))   
                                mean_count_wave_erosion(loop-8:loop)=0;
                                txt = vertcat(txt,mean_count_wave_erosion');    
                            end
                            if(~isempty(mean_count_wea_erosion))    
                                mean_count_wea_erosion(loop-8:loop)=0;
                                txt = vertcat(txt,mean_count_wea_erosion');    
                            end
                        
                            count_tmp = sum(count_vert_wave,2);
                            txt = vertcat(txt,count_tmp');
                        
                        
                            %%% NUMBER OF ERODED BLOCKS
                            txt = vertcat(txt,count_wave');
                            txt = vertcat(txt,count_wave_int');
                            txt = vertcat(txt,count_wea');
                        
                            dlmwrite(filename, txt);
                            clear txt x
                            loop=loop_tmp;
                    
                        end
                    
                        %%% OUTPUT SAVE_PROFILE TO txt
                        if(print_prof_to_txt)
                            extension = '-save_prof.txt';
                            filename=[output_path, profile, num2str(i_angle), connect, wave_height, connect, material_resistance, extension];
                            txt = save_profile;
                            dlmwrite(filename, txt);
                            clear txt
                        end
                    
                        %%% WHOLE EFFECT EROSION
                        if(print_weff_to_jpg)
                            figure; hold on; 
                            subplot(2,1,1);
                            barh(effec_whole_erosion, 1.0);
                            legend('EffectTotalErosion','location','southeast');
                            xmax=max(effec_whole_erosion); xmax=ceil(xmax);
                            if (xmax < 1) xmax = 1; end
                            xlim([0 xmax]); ylim([50 200]);
					
                            subplot(2,1,2);
                            xlim([0 xmax]); ylim([50 200]);
					
                            extension = '-whole_effec_erosion.jpg';
                            filename=[output_path, profile, num2str(i_angle), connect, wave_height, connect, material_resistance, extension];
                            saveas(gcf, filename);
                            close;
                        end
                    
                    					
                        %%% OUTPUT WAVE-TYPE TO JPG
                        if(print_waty_to_txt)
                            xaxis = loop;
                            clear x_w_type1 x_w_type2 x_w_type3 y_w_type1 y_w_type2 y_w_type3
                            totalx_type1 = zeros(1,maxit*tide);
                            totalx_type2 = zeros(1,maxit*tide);
                            totalx_type3 = zeros(1,maxit*tide);
                            totaly_type1 = zeros(1,maxit*tide);
                            totaly_type2 = zeros(1,maxit*tide);
                            totaly_type3 = zeros(1,maxit*tide);
                            x_w_type1     = zeros(1,maxit*tide);
                            x_w_type2     = zeros(1,maxit*tide);
                            x_w_type3     = zeros(1,maxit*tide);
                            y_w_type1     = zeros(1,maxit*tide);
                            y_w_type2     = zeros(1,maxit*tide);
                            y_w_type3     = zeros(1,maxit*tide);
                            count_type1   = zeros(1,tide);
                            count_type2   = zeros(1,tide);
                            count_type3   = zeros(1,tide);
                    
                            extension = '-wave_type.jpg';
                            filename=[output_path, profile, num2str(i_angle), connect, wave_height, connect, material_resistance, extension];
                    
                            wave_type_tmp 	= zeros(tide,maxit);
                            for wt=1:tide
                                tmp = wave_type(wt,1);
                                tmp_num = 2;
                               count_num = 1;
                                move_num = 1;
                                while ( tmp_num < maxit )
                                    while( tmp_num < maxit && tmp == wave_type(wt, tmp_num))
                                        tmp_num = tmp_num + 1;
                                        count_num = count_num + 1;
                                    end
                                    wave_type_tmp(wt,move_num)=count_num;
                                    move_num = move_num + 1;
                                    count_num = 1;
                                    tmp_num = tmp_num + 1;
                                    if(tmp_num<maxit) tmp = wave_type(wt,tmp_num); end
                                end
                            end
                    
                            %%% TXT
                            %total_type1 = x_w_type1;    total_type1 = vertcat(total_type1, y_w_type1);  
                            %total_type2 = x_w_type2;    total_type2 = vertcat(total_type2, y_w_type2);  
                            %total_type3 = x_w_type3;    total_type3 = vertcat(total_type3, y_w_type3);  
    
        					total_type1 = horzcat(totalx_type1',totaly_type1'); 
            				total_type2 = horzcat(totalx_type2',totaly_type2'); 
                			total_type3 = horzcat(totalx_type3',totaly_type3'); 
                    
                            extension1 = '-wave_type1.txt';
                            extension2 = '-wave_type2.txt';
                            extension3 = '-wave_type3.txt';
                            filename=[output_path, profile, num2str(i_angle), connect, wave_height, connect, material_resistance, extension1];
                            dlmwrite(filename,total_type1);
                            filename=[output_path, profile, num2str(i_angle), connect, wave_height, connect, material_resistance, extension2];
                            dlmwrite(filename,total_type2);
                            filename=[output_path, profile, num2str(i_angle), connect, wave_height, connect, material_resistance, extension3];
                            dlmwrite(filename,total_type3);
                    
                            extension1 = '-wave_type11.txt';
                            extension2 = '-wave_type22.txt';
                            extension3 = '-wave_type33.txt';
                            filename=[output_path, profile, num2str(i_angle), connect, wave_height, connect, material_resistance, extension1];
                            dlmwrite(filename,count_type1);
                            filename=[output_path, profile, num2str(i_angle), connect, wave_height, connect, material_resistance, extension2];
                            dlmwrite(filename,count_type2);
                            filename=[output_path, profile, num2str(i_angle), connect, wave_height, connect, material_resistance, extension3];
                            dlmwrite(filename,count_type3);
    
            				extension1 = 'wave_type.txt';
                			filename=[output_path, profile, num2str(i_angle), connect, wave_height, connect, material_resistance, extension1];
                            dlmwrite(filename,wave_type(:,:,1));
                            extension1 = 'wave_type_tmp.txt';
                            filename=[output_path, profile, num2str(i_angle), connect, wave_height, connect, material_resistance, extension1];
                            dlmwrite(filename,wave_type_tmp);
                    
                        end
                    
                    
                        %%% OUTPUT TOTAL EROSION DISTRE
                        if(print_tedi_to_jpg)
                        
                            if (loop ~= maxit)	x=[1:loop];
                            else                x=[1:maxit];
                            end

                			%%% FOR Total Erosion distri
                    		if ( loop ~= maxit )
                        		%total_wave_distri(:,loop+1:maxit)       = [];
                            	%total_submarine_distri(:,loop+1:maxit)	= [];
                                %total_weathering_distri(:,loop+1:maxit)	= [];
                                %total_sum_distri(:,loop+1:maxit)        = [];
                            end
                            x=[1:y_maxcell];			
            
                            figure;
                    		%%% Total MIX distri1
                        	subplot(3,1,1);	hold on;
                            %area(%total_sum_distri(:,loop),x, 'facecolor', 'r');
                            %area(%total_wave_distri(:,loop),x, 'facecolor', 'b');
        					%area(%total_submarine_distri(:,loop),x, 'facecolor', 'm');
            				%area(%total_weathering_distri(:,loop),x, 'facecolor', 'g');
                			legend('total', 'wave', 'submarine', 'weathering', 'location', 'southeast');
                    		hold off;
                    
                            %%% Total MIX distri2
                            subplot(3,1,2); hold on;
                            %p=pplot(%total_sum_distri(:,loop),x);
                            set(p,'EdgeAlpha',0.4);     set(p,'EdgeColor','r');     set(p,'LineWidth',2);
                            %p=pplot(%total_wave_distri(:,loop),x);	
                            set(p,'EdgeAlpha',0.4);     set(p,'EdgeColor','b');     set(p,'LineWidth',2); 
                            %p=pplot(%total_submarine_distri(:,loop),x);
                            set(p,'EdgeAlpha',0.4);     set(p,'EdgeColor','m');     set(p,'LineWidth',2);
                            %p=pplot(%total_weathering_distri(:,loop),x);	
                            set(p,'EdgeAlpha',0.4);     set(p,'EdgeColor','g');     set(p,'LineWidth',2);
                            hold off;
                    
                            %%% Total MIX distri
                            subplot(3,1,3); hold on; 
                            for i=1:round(loop/10):loop	
                                %p=pplot(%total_sum_distri(:,i),x);
                                set(p,'EdgeAlpha',0.4);     set(p,'EdgeColor','r');     set(p,'LineWidth',1);
                                %p=pplot(%total_wave_distri(:,i),x);
                                set(p,'EdgeAlpha',0.4);     set(p,'EdgeColor','b');     set(p,'LineWidth',1);
                                %p=pplot(%total_submarine_distri(:,i),x);
                                set(p,'EdgeAlpha',0.4);     set(p,'EdgeColor','m');     set(p,'LineWidth',1);
                                %p=pplot(%total_weathering_distri(:,i),x);
                                set(p,'EdgeAlpha',0.4);     set(p,'EdgeColor','g');     set(p,'LineWidth',1);
                            end
                            hold off;
					
                            extension = '-total_area_distri.jpg';
                            filename=[output_path, profile, num2str(i_angle), connect, wave_height, connect, material_resistance, extension];
                            saveas(gcf, filename);
                            close;
                    
                            %%% Total Wave force distribution
                            extension = '-total_wave2_distri.jpg';
                            filename=[output_path, profile, num2str(i_angle), connect, wave_height, connect, material_resistance, extension];
                			figure; hold on; 
                    		for i=1:round(loop/10):loop	
                        		%p=pplot(%total_wave_distri(:,i),x);	
                            	set(p,'EdgeAlpha',0.4);     set(p,'EdgeColor','b');     set(p,'LineWidth',1);
                            end
                            saveas(gcf, filename);
                            hold off;	close;
                                                
                        end
					
                                        
                    end %%% RESISTANCE LOOP
				
                    
    				%%%
        			%%% PROFILE COMPARISON
            		%%%
                	if(print_pcom_to_jpg)
					
                        y_line_high(1:x_maxcell)=mwl+ceil(tide/2)-1;
                        y_line_low(1:x_maxcell)=mwl-ceil(tide/2);
                        x_temp=1:x_maxcell;
                    
                        figure();
                        num_count=0;
                        
                        for resi=start_resi:end_resi
                        
                            num_count   = num_count + 1;
                        
                            if(false)
                                x_total     = 1:x_maxcell;    
                                x_total_tmp = 1:x_maxcell;
                                for i=2:y_maxcell
                                    x_total=horzcat(x_total,x_total_tmp);
                                end
                                y_total=final_prof(1:x_maxcell,1,resi)';
                                for i=2:y_maxcell
                                    y_total_tmp=i*final_prof(1:x_maxcell,i,resi)';
                                    y_total=horzcat(y_total,y_total_tmp);
                                end
                            end
                        
                           plot_num = end_resi-start_resi+1;
                            subplot(plot_num,1,num_count); hold on;
                            plot(final_prof(:,resi),[1:y_maxcell] , 'b', 'linewidth',2); hold on; 
                            plot(x_temp, y_line_high, '--r', x_temp, y_line_low, '--r');
                            text(2,y_maxcell-(0.5*y_maxcell/10+1),['Rock Hardness=',num2str(resi_values(resi))],'Fontsize',10,'Backgroundcolor',[.7 .9 .7]);
                            text(2,y_maxcell-(2.5*y_maxcell/10+1),['Iteration=',num2str(final_it(resi))],'Fontsize',10,'Backgroundcolor',[.7 .9 .7]);
                            axis([1 x_maxcell 1 y_maxcell]);
										
                        end
            
                        if (h==1)		extension = [num2str(height_values(1)),'m-',num2str(resi_values(1)),'.jpg'];
                        elseif (h==2)	extension = [num2str(height_values(2)),'m-',num2str(resi_values(1)),'.jpg'];
                        elseif (h==3) 	extension = [num2str(height_values(3)),'m-wave.jpg'];
                        elseif (h==4)	extension = [num2str(height_values(4)),'m-',num2str(height_values(5)), 'm-wave1.jpg'];
                        elseif (h==5) 	extension = [num2str(height_values(4)),'m-',num2str(height_values(5)), 'm-wave2.jpg'];
                        elseif (h==6)	extension = [num2str(height_values(4)),'m-',num2str(height_values(5)), 'm-wave3.jpg'];
                        end
				    
                        filename=[output_path, tidal_name(tr), connect, 'P', num2str(pressure), connect, 'wea-', num2str(wea), connect, 'weasd-', num2str(wt_const), '-co-', num2str(co), connect, extension];
                        saveas(gcf, filename);
                        output_path2 = [current_path, num2str(i_angle), '/profile_comp', '-grid-', num2str(gridsize_values(gridloop)), '/'];
                        filename2=[output_path2, tidal_name(tr), connect, 'P', num2str(pressure), connect, 'wea-', num2str(wea), '-weasd-', num2str(wt_const), '-co-', num2str(co),'-', extension];
                        saveas(gcf, filename2);
                        hold off;
                        close;
                    end
                   
                end % WAVE HEIHGT LOOP
            end % CO LOOP
            
            %%% STORE MWL
            filename=[current_path, 'MTR-', tidal_name(tr), 'm-mwl.csv'];
            dlmwrite(filename,store_mwl);
            
            end % WEA SPEED LOOP
            end % WEA LOOP
            
           
        end % PRESSURE DISTRIBUTION LOOP
    
        %%% Pint final info
    	final_output_path = [current_path, num2str(i_angle),'/'];
        filename = [final_output_path, 'tide_', num2str(tidal_values(tr)), '_grid_', num2str(gridsize_values(gridloop)), '_final.txt'];
        dlmwrite(filename, final_info);
        
                        
    end % TIDAL RANGE LOOP
    
end % PROFILE LOOP


end

endT = clock();
ntime = etime(endT,startT);
nhour = floor(ntime/60/60);
nmin = floor((ntime-nhour*3600)/60);
nsec = ntime-nhour*3600-nmin*60;
disp(sprintf('%s%s', 'Start time  =  ',datestr(startT,31)));
disp(sprintf('%s%s', 'Finish time  =  ',datestr(endT,31)));
disp(sprintf('%s%d%s%02d%s%04.1f%s', 'Elapsed time  =  ',nhour,' hour  ',nmin,' min  ',nsec,' sec  '));

sec(gridloop) = nhour*3600+nmin*60+nsec;
filename = [final_output_path, 'tide-', num2str(tidal_values(tr)),'grid-', num2str(gridsize_values(gridloop)), '-num_type.txt'];
out = [current_path, num2str(gridloop)];
dlmwrite(filename,sec);
