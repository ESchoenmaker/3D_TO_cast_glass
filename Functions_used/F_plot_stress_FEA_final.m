function F_plot_stress_FEA_final(x_L, neigh_dis_L, neigh_sum_L,Kel_L, dofs_free_L, F_L, edofs_L)

[p_ym_L, p_ymmin_L, p_nu_L, ~, ~, ~, ~, d_FE_L, ~,...
    ~, ~, ~, nelx_L, nely_L, nelz_L, tot_elements_L, ~, ~, pen_L, ~, ~, ~, name_run_L, ~, q_L, ~, ~, ~, ~, ~] = F_data_case;

xPhys_L = x_L;
xPhys_L(:) = (neigh_dis_L*x_L(:))./neigh_sum_L;

% Stiffness Matrix Assembly
dofs_all_L = 1:3*(nelx_L+1)*(nely_L+1)*(nelz_L+1);
nodegrd = reshape(1:(nelz_L+1)*(nelx_L+1),nelx_L+1,nelz_L+1);
nodegrd = nodegrd';
nodeids2 = reshape(nodegrd(2:end,1:end-1),nelz_L*nelx_L,1);
nodeidz = 0:(nelz_L+1)*(nelx_L+1):(nely_L)*(nelz_L+1)*(nelx_L+1);
nodeidz = nodeidz(:, 2:end);
nodeids = repmat(nodeids2,size(nodeidz))+repmat(nodeidz,size(nodeids2));
edofVec = 3*(nodeids(:)-1)+1;
edofMat = repmat(edofVec,1,24)+ ...
    repmat([0 1 2 -3*(nelx_L+1)*(nelz_L+1) + [0 1 2 3 4 5] 3 4 5 ...
    -3*(nelx_L+1) + [0 1 2] -3*((nelx_L+1)*(nelz_L+1)+(nelx_L+1)) + [ 0 1 2 3 4 5] -3*(nelx_L+1) + [3 4 5]],tot_elements_L,1);
iK = reshape(kron(edofMat,ones(24,1))',24*24*tot_elements_L,1);
jK = reshape(kron(edofMat,ones(1,24))',24*24*tot_elements_L,1);
sK = reshape(Kel_L(:)*(p_ymmin_L+xPhys_L(:)'.^pen_L*(p_ym_L-p_ymmin_L)),24*24*tot_elements_L,1); % - [kN/mm2]
Ktot = sparse(iK,jK,sK); Ktot = (Ktot+Ktot')/2;

% Finite Element Analysis 
U_L = zeros(numel(dofs_all_L),1); 
tolit = 1e-8;
maxit = 8000;
if tot_elements_L >= 10000
    M = diag(diag(Ktot(dofs_free_L, dofs_free_L)));
    U_L(dofs_free_L, :) = pcg(Ktot(dofs_free_L, dofs_free_L),F_L(dofs_free_L, :), tolit, maxit, M);
else
    U_L(dofs_free_L, :) = Ktot(dofs_free_L, dofs_free_L)\ F_L(dofs_free_L, :);              % - [mm2] 
end

[max_tens_stress_L, max_comp_stress_L,princ_stress_tens, princ_stress_comp] = F_stress_calc_gauss(tot_elements_L,  p_ym_L, p_ymmin_L, p_nu_L, edofs_L, U_L, d_FE_L, pen_L, xPhys_L);
weighted_stress_tens = reshape(xPhys_L, [numel(xPhys_L), 1]).^(pen_L - q_L).*max_tens_stress_L; %kN/mm2
max_tens_stress = max(weighted_stress_tens)
max_tens_stress_non_w = max(max_tens_stress_L)
max_comp = min(max_comp_stress_L) %kN/mm2
%plotting
folder = 'Results';

%plot tensile stress
max_princ = figure('Visible', 'off');
set(gcf, 'Position', [0,0,1920,1080]);
view(-30, 30);
F_plot_stress_new(nelx_L, nely_L, nelz_L, d_FE_L, princ_stress_tens, tot_elements_L, xPhys_L) %when plotting tensile blue is lowest, for compressive blue is highest
axis equal; 
axis padded;
axis off;
fname = "figure_final_max_princ.jpg";
exportgraphics(max_princ, fullfile(folder,name_run_L, fname), Resolution=300)

%plot compression stress
min_princ = figure('Visible', 'off');
set(gcf, 'Position', [0,0,1920,1080]);
view(-30, 30);
F_plot_stress_new(nelx_L, nely_L, nelz_L, d_FE_L, princ_stress_comp, tot_elements_L, xPhys_L) %when plotting tensile blue is lowest, for compressive blue is highest
axis equal; 
axis padded;
axis off;
fname = "figure_final_min_princ.jpg";
exportgraphics(min_princ, fullfile(folder,name_run_L, fname), Resolution=300)

