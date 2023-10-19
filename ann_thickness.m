% Function to calculate thickness based on maximum annealing time
function d_ann = ann_thickness(anneal_t_L, p_den_L, p_nu_L, p_ym_L)
anneal_range = 530-460; % - Borosilicate (degrees Celsius or Kelvin)
p_ym_L = p_ym_L*(10^3); % - MPa
b = 0.3; % - shape factor - (-)
l = 1.15; % - thermal conductivity (borosilicate glass) - (W/m*K)
cp = 800; % - specific heat capacity (borosilicate glass) - (J/kg*K)
max_perm_stress = 1; % - MPa
th_coeff = 3.25*(10^(-6)); % - Borosilicate (K^-1)
d_ann = sqrt((anneal_t_L*max_perm_stress)/((anneal_range*p_ym_L*th_coeff*p_den_L*cp*b)/((1-p_nu_L)*l))); % - m
end