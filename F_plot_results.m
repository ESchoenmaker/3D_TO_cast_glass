function F_plot_results(name_run_L, nelx_L, nely_L, nelz_L, d_FE_L, xi_L, itt_L, neigh_dis_L, neigh_sum_L)

% Apply Filter
xPhys = xi_L;
xPhys(:) = (neigh_dis_L*xi_L(:))./neigh_sum_L;

nele_L = nelx_L*nely_L*nelz_L;                                    % Total amount elements

% Setup face matrix
A = (0:8:(nele_L-1)*8);
B = repmat(reshape(repmat(A, 6, 1), [nele_L*6, 1]), 1, 4);
face_mat = [1 2 3 4; 2 6 7 3; 6 5 8 7; 5 1 4 8; 1 2 6 5; 4 3 7 8]; % Connection table
face_full = B + repmat(face_mat, nele_L, 1);

% Setup vert matrix
vert_mat = [0 0 0; d_FE_L 0 0; d_FE_L 0 -d_FE_L; 0 0 -d_FE_L;...
                    0 -d_FE_L 0;d_FE_L -d_FE_L 0; d_FE_L -d_FE_L -d_FE_L;0 -d_FE_L -d_FE_L];
x_mat = repmat(reshape(repmat((0:(nelx_L - 1))*d_FE_L, nelz_L*8, 1), [], 1), nely_L, 1);
y_mat = reshape(repmat((0:(nely_L-1))*d_FE_L*-1, (nelz_L*nelx_L*8), 1), [], 1);
z_mat = repmat(reshape(repmat(reshape(repmat((0:(nelz_L-1))*d_FE_L*-1, 8, 1), [], 1), 1, nelx_L), [], 1), nely_L, 1);

coords_mat = cat(2, x_mat, y_mat, z_mat);
vert_tot = coords_mat + repmat(vert_mat, nele_L, 1);

% Setup Colours and Oppacity 
colour = repmat(reshape(repmat(reshape((1-xPhys), 1, []), 6,1), [], 1), 1, 3);
face_alpha = zeros(size(colour(:,1)));                  % Setting up transparency of elements 
face_alpha(colour(:,1)<=0.7) =1;                        % Elements not shown if density lower then 0.3

% Setup Window
density = figure('Visible', 'off');                     % Visibilty turned of for DHPC calculations
set(gcf, 'Position', [0,0,1920,1080]);
view(-30, 30);
%view(3);
axis equal;
axis padded;
axis off;

% Plotting
patch('Faces',face_full,'Vertices',vert_tot, 'FaceVertexCData', colour, 'FaceVertexAlphaData', face_alpha,'FaceAlpha', 'flat', 'EdgeAlpha', 'flat', 'FaceColor','flat');  % Activate for calculations where cubes with density lower than 0.3 should not be shown
%patch('Faces',face_full,'Vertices',vert_tot, 'FaceVertexCData', colour, 'FaceColor','flat');  % Activate for plotting where all elements are shown regardless of density

% Setup Legend
colormap(flipud(gray))
c = colorbar('Ticks',[0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1]);
clim([0 1]);
c.Label.String = 'Density';

% Setup Saving
fname = "figure_plt_new" + itt_L + ".jpg";
folder = 'Results';
density = gcf;
exportgraphics(density, fullfile(folder,name_run_L, fname), Resolution=300)
end

