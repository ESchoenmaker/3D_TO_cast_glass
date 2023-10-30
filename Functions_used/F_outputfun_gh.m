function stop = F_outputfun_gh(x, optimValues, ~)
%Output function that is integrated into the optimization
%  Includes plotting and Writing data to excel  

stop = false;   %if true optimization stops, if false optimization continues

% Get data for plotting
itt = optimValues.iteration; 
f_count = optimValues.funccount;
fval = optimValues.fval;
feas = optimValues.constrviolation;
first_ord = optimValues.firstorderopt;
tend = toc;

[p_ym, p_ymmin, p_nu, p_den, d_Ln, d_wd, d_ht, d_FE, vole,...
    ~, ~, ~, nelx, nely, nelz, ~, dmin, dmax, pen, ac, vol_frac,initial_v, name_run, beta, q, tens_limit, comp_limit, max_u, loadcase_1, ~] = F_data_gh; 
[neigh_dis, neigh_sum] = F_filtering(nelx, nely, nelz, dmin, d_FE); 

% In case of starting the algorithm from a previous iteration
% itt = itt + 45;

% Activate in case of DPHC - numbers get writen to output slurm file
display = ['itt: ', num2str(itt), ' F_count: ', num2str(f_count), ' f(x): ', num2str(fval), ' feasibility: ', num2str(feas), ' first_order ', num2str(first_ord), ' time ', num2str(tend)];
disp(display);

% Saving final density matrix
folder = 'Results';
% name_overwrite = 'Density_final.mat';
% save(fullfile(folder, name_run, name_overwrite), 'x', '-mat')
% Saving density matrix each iteration
formatSpec = 'Density_it_%d.mat';
m_name = sprintf(formatSpec, itt);
save(fullfile(folder, name_run, m_name), 'x', '-mat')

% Plotting
%F_plot_results(name_run, nelx, nely, nelz, d_FE, x, itt, neigh_dis, neigh_sum) 
%F_plot_results_mirror(name_run, nelx, nely, nelz, d_FE, x, itt, neigh_dis, neigh_sum)

% Writing to excel
F_write_excel(x, p_ym, p_ymmin, p_nu, p_den, d_Ln, d_wd, d_ht, d_FE, vole, nelx, nely, nelz, dmin, dmax, pen, ac, vol_frac,initial_v, name_run, beta,...
    itt, f_count, fval, feas, first_ord, tend, loadcase_1, q)
end