%Compliance Optimization with boundary conditions read from grasshopper
[p_ym, p_ymmin, p_nu, p_den, d_ln, d_wd, d_ht, d_FE, vole,...
    nodx, nody, nodz, nelx, nely, nelz, tot_elements, dmin, dmax, pen, ac,vol_frac, initial_v, name_run, beta, q, tens_limit, comp_limit, max_u, loadcase_1, av, name_excel] = F_data_gh;

tic;                                                                            % Starting timer per loop

%%%% Boundary Conditions %%%%
dofs_all = 1:3*numel(nodx)*numel(nody)*numel(nodz);  
%loads
load_element = loadcase_1 / (nelx*nely); 
[loadnid_middle, loaddof_middle] = f_read_nodes_gh(name_excel, 'Loads_middle', d_FE, nodx, nodz, nody, 0);
[loadnid_sides, loaddof_sides] = f_read_nodes_gh(name_excel, 'Loads_sides', d_FE, nodx, nodz, nody, 0);
F = zeros(numel(dofs_all), 1);                                                  
F(loaddof_middle, 1) = -1*load_element;                                         % [kN / per node]
F(loaddof_sides, 1) = -1*load_element/2;      

%supports
[nodes_fixed_x, dofs_fixed_x] = f_read_nodes_gh(name_excel, 'Supports_fixed_x', d_FE, nodx, nodz, nody, -1); % [y x z]
[nodes_fixed_z, dofs_fixed_z] = f_read_nodes_gh(name_excel, 'Supports_fixed_z', d_FE, nodx, nodz, nody, 0); % [y x z]
[nodes_mirror_x, dofs_mirror_x] = f_read_nodes_gh(name_excel, 'Supports_mirror_x', d_FE, nodx, nodz, nody, -1);
[nodes_mirror_y, dofs_mirror_y] = f_read_nodes_gh(name_excel, 'Supports_mirror_y', d_FE, nodx, nodz, nody, -2);

dof_fixed = sort(cat(1, dofs_fixed_x, dofs_fixed_z));       
dof_fixed = sort(cat(1, dof_fixed, dofs_mirror_x));                               % Activate this in case there is mirroring in the x direction
dof_fixed = sort(cat(1, dof_fixed, dofs_mirror_y));                             % Activate this in case there is mirroring in the Z direction
dofs_free = setdiff(dofs_all, dof_fixed);                                   % Set free nodes

%%%% Preprocessing Step %%%%
%construct starting point
xi(1:nelz, 1:nelx, 1:nely) = initial_v;                                         % Starting density 
% xi = deal(load('Density_it_45.mat').x);                                         % Activate in case of starting optimization from previous density map                          

%constructing matrix neighbours
neigh_ann = F_filtering_ann(nelx, nely, nelz, dmax, d_FE);                      % Annealing neighbours
[neigh_dis, neigh_sum] = F_filtering(nelx, nely, nelz, dmin, d_FE);             % Filtering neighbours

%constructing stiffness matrix
Kel = F_StiffnesMatrix (1,p_nu,d_FE,d_FE,d_FE);                                 % E imput is set to 1 so this can be iterated in the optimization

%Constructing connectivity matrix
[connectivity, edofs] = F_connectivity(nelx, nely, nelz, tot_elements, nodx, nodz);

%%%% Optimization Setup %%%%
f = @(x)ObjFcn(x, nelx, nely, nelz, ...
    nodx, nody, nodz, p_ymmin, p_ym, pen, Kel, dofs_free, F, neigh_dis, neigh_sum);
lb(1:nelz, 1:nelx, 1:nely) = 0.001;
ub(1:nelz, 1:nelx, 1:nely) = 1;
fhes = @(x, lambda)HessianFcn(x, lambda, nely, nelx, nelz, pen, p_ym, p_ymmin, neigh_dis, neigh_sum);
constr = @(x)CnstrFcn(x, vol_frac, neigh_dis, neigh_sum, neigh_ann,...
    nelx, nely, nelz, Kel,p_ymmin, pen, p_ym, dofs_free, F, p_nu, edofs,d_FE, tens_limit, comp_limit, q, nodx, max_u, connectivity, vole, av);

options_display = optimoptions('fmincon','SpecifyConstraintGradient', true, 'SpecifyObjectiveGradient',true,'Hessian', 'user-supplied', 'HessFcn',fhes,  'Display', 'iter', 'PlotFcn', ...
   'optimplotfval', 'OutputFcn', @F_outputfun_gh, 'MaxIterations', 1000, 'MaxFunctionEvaluations', 700000, ...
   'EnableFeasibilityMode', true, 'SubproblemAlgorithm', "cg", ...
   'Algorithm', 'interior-point', 'StepTolerance', 1e-19,   'BarrierParamUpdate', 'predictor-corrector' );

options_nodisplay = optimoptions('fmincon','SpecifyConstraintGradient', true, 'SpecifyObjectiveGradient',true,'Hessian', 'user-supplied', 'HessFcn',fhes,  'Display', 'none', ...
   'OutputFcn', @F_outputfun_gh, 'MaxIterations', 10000, 'MaxFunctionEvaluations', 700000, ...
   'EnableFeasibilityMode', true, 'SubproblemAlgorithm', "cg", ...
   'Algorithm', 'interior-point', 'StepTolerance', 1e-9,   'BarrierParamUpdate', 'predictor-corrector' );

[x, fval] = fmincon(f, xi, [],[],[],[],lb,ub,constr, options_nodisplay);

%%%% Post optimization Process %%%%
folder = 'Results';
Optimization_plot = gcf;
fname = "Optimization_plot_" + name_run + ".jpg";
exportgraphics(Optimization_plot, fullfile(folder,name_run, fname))                         % Saving optimization graph

F_plot_results_mirror(name_run, nelx, nely, nelz, d_FE, x, name_run, neigh_dis, neigh_sum)  % Plot resulting denstiy map
F_plot_stress_FEA_final_gh(x, neigh_dis, neigh_sum,Kel, dofs_free, F, edofs)     % Plot resulting principal stress maps

%%%% Setup Global Values needed for hessian function %%%%
function setGlobalval(val)
global x
x = val;
end
function return_global = getGlobal
global x
return_global = x;
end


%%%% Constraint function %%%%%
function [c, ceq, gradc, gradceq] = CnstrFcn(xi_L, vol_frac_L, neigh_dis_L, neigh_sum_L, neigh_ann_L,...
    nelx_L, nely_L, nelz_L, Kel_L,p_ymin_L, pen_L, p_ym_L, dofs_free_L, F_L, nu_L, edofs_L,d_FE_L, tens_limit_L, comp_limit_L, q_L, nodx_L, max_u_L, connectivity_L, vole_L, av_L)
% Apply filter
xPhys = xi_L;
xPhys(:) = (neigh_dis_L*xi_L(:))./neigh_sum_L;

% Volume constraint 
vf_L = sum(xPhys.*vole_L, "all");
c_vol = vf_L/(vol_frac_L*vole_L*numel(xPhys)) - 1;                              % Volume constraint function
dv = (ones(numel(xPhys),1)*vole_L)./(numel(xPhys)*vole_L*vol_frac_L);           % Applying Derivative of filtered density 
grad_vol = neigh_dis_L*(dv(:)./neigh_sum_L);                                    % Volume gradient function

% Annealing constraint 
Hann = neigh_ann_L*reshape(xPhys, [numel(xPhys), 1])*vole_L;
Hann_limit = sum(neigh_ann_L.*vole_L, 2) .*av_L;
c_ann = Hann./Hann_limit -1;                                                    % Annealing constraint function
grad_ann = sparse((neigh_ann_L.*vole_L) ./Hann_limit');                         % Annealing gradient function

% Stiffness Matrix Assembly
tot_elements_L = nelx_L*nely_L*nelz_L;
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
sK = reshape(Kel_L(:)*(p_ymin_L+xPhys(:)'.^pen_L*(p_ym_L-p_ymin_L)),24*24*tot_elements_L,1); % - [kN/mm2]
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

% stress constraint 
[max_tens_stress, max_comp_stress, ~, ~] = F_stress_calc_gauss(tot_elements_L,  p_ym_L, p_ymin_L, nu_L, edofs_L, U_L, d_FE_L, pen_L, xPhys);
max_tens_stress = max(max_tens_stress,0);                                                  % - [kN/mm2]
max_comp_stress = min(max_comp_stress,0);                                                  % - [kN/mm2]
max_comp_stress = abs(max_comp_stress); 
c_tens = reshape(xPhys, [numel(xPhys), 1]).^(pen_L - q_L).*(max_tens_stress./tens_limit_L) -1;                      % Tensile stress constraint
c_comp = reshape(xPhys, [numel(xPhys), 1]).^(pen_L - q_L).*(max_comp_stress./comp_limit_L) -1;                      % Compression stress constraint
% gradient stress constraint 
grad_tens = (pen_L - q_L).* reshape(xPhys, [numel(xPhys), 1]).^(pen_L - q_L-1).*(max_tens_stress./tens_limit_L);    % Tensile stress gradient as vector
grad_tens = spdiags(grad_tens, 0, tot_elements_L, tot_elements_L);                                                  % Translating to right shape
grad_comp = (pen_L - q_L).* reshape(xPhys, [numel(xPhys), 1]).^(pen_L - q_L-1).*(max_comp_stress./comp_limit_L);    % Compression stress gradient as vector
grad_comp = spdiags(grad_comp, 0, tot_elements_L, tot_elements_L);                                                  % Translating to right shape

% displacement constraint
critnode = numel(nodx_L);                                                               % Top right front node
critdof = critnode *3;                                                                  % Deformation in z direction
ucrit = abs(U_L(critdof));  
c_dis = ucrit/max_u_L - 1;                                                              % Displacement constraint function
% displacement gradient
[index_el, ~] = find(connectivity_L == critnode);                                       % Find which elements the critnode is part of
reshaped_dof = reshape(edofs_L, [24, tot_elements_L])';
index_dofs = reshaped_dof(index_el, :);                                                 % Find DOFs of elements
reshaped_xPhys = reshape(xPhys, [numel(xPhys), 1]);
grad_dis_temp = zeros(numel(index_el), 24);
for i=1:numel(index_el)
    bottom = Ktot(index_dofs(i,:), index_dofs(i,:));                                                                 % K(x)^-1 
    top = U_L(index_dofs(i, :))'*((p_ymin_L+pen_L*reshaped_xPhys(index_el(i))^(pen_L-1)*(p_ym_L-p_ymin_L))*Kel_L);   % U(x)[px(p-1)(E-E_min)K]
    grad_dis_temp(i, :) = -1.*(bottom \ top');                                                                       % Calculation gradient per Element individually
end
iU = reshape(repmat(index_el, 1,24)', [24*numel(index_el), 1]);
jU = reshape(index_dofs', [24*numel(index_el), 1]);
sU = reshape(grad_dis_temp',[24*numel(index_el), 1]);
grad_dis_tot = sparse(iU, jU, sU);	                                                     % Combining gradient elements into one matrix
grad_dis = zeros(tot_elements_L, 1);
grad_dis(index_el, 1) = grad_dis_tot(index_el, critdof)./max_u_L;                        % Displacement Gradient function 

% Choose constraint combination
%c =  c_vol;
%gradc = grad_vol;
c = cat(1, c_vol,c_ann, c_tens, c_comp, c_dis);
gradc = cat(2, grad_vol, grad_ann, grad_tens, grad_comp, grad_dis);

ceq = []; 
gradceq = [];

% Print statements to check if constraints are met during the optimization
if c_vol > 0
    disp('volume constraint not met')
end

if c_dis > 0
    disp('Displacement constraint not met')
end

check = all(c_ann(:)<0);
if check == 0
    disp('Annealing constraint not met')
end
check = all(c_tens(:)<0);
if check == 0
    disp('Tensile stress constraint not met')
end
check = all(c_comp(:)<0);
if check == 0
    disp('Compressive stress constraint not met')
end

end %CnstrFcn; returns: c, gradc, ceq, gradceq

%%%% Obective function %%%%%
function [cmpl, gradf] = ObjFcn(xi_L, nelx_L, nely_L, nelz_L,....
nodx_L, nody_L, nodz_L, p_ymin_L, p_ym_L, pen_L, Kel_L, dofs_free_L, F_L,neigh_dis_L, neigh_sum_L) 
% Apply filter
xPhys = xi_L;
xPhys(:) = (neigh_dis_L*xi_L(:))./neigh_sum_L;

% Stiffness Matrix Assembly
tot_elements_L = nelx_L*nely_L*nelz_L;
dofs_all_L = 1:3*numel(nodx_L)*numel(nody_L)*numel(nodz_L);
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
sK = reshape(Kel_L(:)*(p_ymin_L+xPhys(:)'.^pen_L*(p_ym_L-p_ymin_L)),24*24*tot_elements_L,1);
Ktot = sparse(iK,jK,sK); Ktot = (Ktot+Ktot')/2;

% Finite Element Analysis
U_L = zeros(numel(dofs_all_L),1);
tolit = 1e-8;
maxit = 8000;
if tot_elements_L >= 10000
M = diag(diag(Ktot(dofs_free_L, dofs_free_L)));
U_L(dofs_free_L, :) = pcg(Ktot(dofs_free_L, dofs_free_L),F_L(dofs_free_L, :), tolit, maxit, M);
else
U_L(dofs_free_L, :) = Ktot(dofs_free_L, dofs_free_L)\ F_L(dofs_free_L, :);          % [mm]
end

% Calculating Compliance
ce = reshape(sum((U_L(edofMat)*Kel_L).*U_L(edofMat), 2), [nelz_L, nelx_L, nely_L]); % Gives compliance without taking the physicial(filtered) densities into account
setGlobalval(ce)
c = sum(sum(sum((p_ymin_L + xPhys.^pen_L*(p_ym_L-p_ymin_L)).*ce)));  
dc = -pen_L*(p_ym_L-p_ymin_L)*xPhys.^(pen_L-1).*ce;                                 % Gradient function first part is the [p*xPhys^(p-1) * (E0 - Emin)] part ce = [u^t * K * u]
dc(:) = neigh_dis_L*(dc(:)./neigh_sum_L);                                           % Filtering the gradient function 

cmpl = c;                                                                           % Objective function
gradf = dc(:);                                                                      % Objective Gradient
end % ObjFcn; returns: cmpl, gradf

%%%% Hessian Function %%%%
function h = HessianFcn(x, lambda, nely_L, nelx_L, nelz_L, pen_L, p_ym_L, p_ymin_L, neigh_dis_L, neigh_sum_L)
    return_ce = getGlobal; 
    xPhys = reshape(x,nelz_L,nelx_L,nely_L);
    Hessf = 2*(pen_L*(p_ym_L-p_ymin_L)*xPhys.^(pen_L-1)).^2 ./ (p_ym_L + (p_ym_L-p_ymin_L)*xPhys.^pen_L) .* return_ce;     % Compute Hessian of Obj
    Hessf(:) = neigh_dis_L*(Hessf(:)./neigh_sum_L);                                                                        % Compute Hessian of Obj
    Hessc = 0;                                                                                                             % Compute Hessian of constraints - 0 for linear constraint
    % Hessian of Lagrange
    h = diag(Hessf(:)); %+ lambda.ineqnonlin*Hessc; % blocked the second part out as it was giving an error and hessc is 0 
end % Returns: h

