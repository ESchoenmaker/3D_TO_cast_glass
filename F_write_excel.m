function F_write_excel(xi_L, p_ym_L, p_ymmin_L, p_nu_L, p_den_L, d_Ln_L, d_wd_L, d_ht_L, d_FE_L, vole_L, nelx_L, nely_L, nelz_L, dmin_L, dmax_L, pen_L, ac_L, vol_frac_L, initial_v_L, name_run_L, beta_L,...
    itt_L, f_count_L, fval_L, feas_L, first_ord_L, tend_L, loadcase_L, q_L)
    %writes results to excel file
%change 3D matrix into 2D
xi_cell = num2cell(xi_L, [1 2]);
xi_cell = reshape(xi_cell, size(xi_L,3), 1);
desired = cell2mat(xi_cell);

%writing it to the excel
filename = name_run_L;
folder = 'Results';
fullFileName = fullfile(folder,filename, filename);
fullFileName = append(fullFileName, '.xlsx');
sheet1 = 1;  
start_t = "A"+1; 
start_m = "A"+3; 
input_arg = table(p_ym_L, p_ymmin_L, p_nu_L, p_den_L, d_Ln_L, d_wd_L, d_ht_L, d_FE_L, vole_L, nelx_L, nely_L, nelz_L, dmin_L, dmax_L, pen_L, ac_L, vol_frac_L,initial_v_L, beta_L, loadcase_L, q_L);

writetable(input_arg, fullFileName, 'Sheet', sheet1, 'Range', start_t);
writematrix(desired, fullFileName, 'Sheet', sheet1, 'Range', start_m); 

%write information optimization loop
sheet2 = 2;
input_arg_2 = table(itt_L, f_count_L, fval_L, feas_L, first_ord_L, tend_L);
start_itt_1 = "A"+1;
start_itt_n = "A"+(itt_L+1);
if itt_L == 1
    writetable(input_arg_2, fullFileName, 'Sheet', sheet2, 'Range', start_itt_1)
else
    writetable(input_arg_2, fullFileName, 'Sheet', sheet2, 'Range', start_itt_n, 'WriteVariableNames', false)
end
end