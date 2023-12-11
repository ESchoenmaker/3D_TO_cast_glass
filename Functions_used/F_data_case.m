function [p_ym_L, p_ymmin_L, p_nu_L, p_den_L, d_ln_L, d_wd_L, d_ht_L, d_FE_L, vole_L,...
    nodx_L, nody_L, nodz_L, nelx_L, nely_L, nelz_L, tot_elements_L,dmin_L, dmax_L, pen_L, ac_L, vol_frac_L, initial_v_L, name_run_L, beta_L, q_L, tens_limit_L, comp_limit_L, max_u_L, loadcase_1_L, av_L] = F_data_case()
%data
p_ym_L = 70 ;                   % [kN/mm2]  Young's modulus - 70 GPa = 70 (kN/mm2) = 70 *10^6 [kN/m2]
p_ymmin_L = 0.001;              % min Young's modulus to avoid singularity - GPa [kN/mm2]  
p_nu_L = 0.2;                   % Poisson's ratio
p_den_L = 2500 /10^9;           % Density - [kg/mm3]
%Finite Element Dimensions
d_FE_L =  20;                   % [mm] 
vole_L = d_FE_L^3;              % [mm3]
%input dimensions
d_ln_L = 2100;%4200;                  % [mm]
d_wd_L = d_FE_L; %800;             % [mm]
d_ht_L = 300;                   % [mm]
% nodes
nodx_L = 0 : d_FE_L : d_ln_L;   % [mm]
nody_L = 0 : d_FE_L : d_wd_L;   % [mm]
nodz_L = 0 : d_FE_L : d_ht_L;   % [mm]
% elements
nelx_L = numel(nodx_L)-1;
nely_L = numel(nody_L)-1;
nelz_L = numel(nodz_L)-1;
tot_elements_L = nelx_L*nely_L*nelz_L;
% Final size
d_ln_L=nelx_L*d_FE_L;           % [mm]
d_wd_L=nely_L*d_FE_L;           % [mm]
d_ht_L=nelz_L*d_FE_L;           % [mm]
%loads
Fsw = p_den_L*d_ht_L*0.0098;    % kN/mm2 permanent load (selfweight of the slab) 
Ffl = p_den_L*20*0.0098 ;       % kN/mm2 permanent load (2 10 mm float glass sheets
Fsp = 5 /10^6;                  % kN/mm2 short term load (people) 
Fso = 0.4 /10^6 ;               % kN/mm2 short term load (maintenance)
Sp = 1.2;                       % Safety factor permanent
Ss = 1.5;                       % Safety factor short term
loadcase_1_L = (Sp*(Fsw+Ffl)*d_ln_L*d_wd_L + Ss*(Fso + Fsp)*d_ln_L*d_wd_L); % [kN] 
%filtering 
dmin_L = 3*d_FE_L;              % [mm]
dmin_L = dmin_L*1.25;           % [mm]
%annealing
max_anneal_T_L = 5*24*60*60;    % [sec]
d_ann_L = ann_thickness(max_anneal_T_L, p_den_L*10^9, p_nu_L, p_ym_L)*1000; %Result is in meters without addition of *1000
dmax_L = min(d_ann_L, 2*dmin_L);
%Watch out chaning this - the value is calculated backwards to a radius instead of diameter in some of the filtering functions
%missing
pen_L = 3 ;                     % penalization value
ac_L = 3.5;                     % compliance constraint
q_L = 2.8 ;                     % relaxation parameter
vol_frac_L = 0.3;               % volume constraint fraction
initial_v_L = 0.4;              % Initial volume fraction
av_L = 0.9 ;                    % allowable percentage volume for annealing constraint
%Heaviside function
beta_L = 0;
%limits 
tens_limit_L = 6.4*1000 /10^6;   % [kN/mm2]
comp_limit_L = 500*1000 /10^6;   % [kN/mm2]
max_u_L = d_ln_L*2/500 ;         % [mm]
%setting up data saving
name_run_L = '231107_Case_gif_cmpl';
name_folder = 'Results\';
dir = append(name_folder, name_run_L);
if not(isfolder(dir))
   mkdir("Results",name_run_L)
end

