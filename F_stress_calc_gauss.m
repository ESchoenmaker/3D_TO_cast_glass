function [max_tens_stress_L, max_comp_stress_L, princ_tens_L, princ_comp_L] = F_stress_calc_gauss(tot_elements_L,  E_L, Emin_L, nu_L, edofs_L, U_L, d_FE, pen_L, xPhys_L)
% Function to calculate the maximum tensile and compressive stress
% Function uses principle stresses 
% Function is vectorized for speed 

%reshape xPhys 
xPhys_Lx = reshape(xPhys_L, [numel(xPhys_L), 1]);
%empty array where max stress (comp and tens) is going to be stored per element
max_tens_stress_L = zeros(tot_elements_L, 1);
max_comp_stress_L = zeros(tot_elements_L, 1);
princ_tens_L = zeros(tot_elements_L, 8);
princ_comp_L = zeros(tot_elements_L, 8);
%Material stiffness matrix that is needed later on
D = [1-nu_L nu_L nu_L 0 0 0;...   
    nu_L 1-nu_L nu_L 0 0 0;...
    nu_L nu_L 1-nu_L 0 0 0;...
    0 0 0 (1-2*nu_L)/2 0 0;...
    0 0 0 0 (1-2*nu_L)/2 0;...
    0 0 0 0 0 (1-2*nu_L)/2];

Exi = Emin_L + xPhys_Lx.^pen_L*(E_L - Emin_L);  %including penalty [kN/mm2]
Exi_new = reshape(Exi, [1, 1, numel(Exi)])./((1+nu_L)*(1-2*nu_L)); 
D_comp = Exi_new.*D;          %setting up (6*6*nele) size matrix [kN/mm2]

%create 3d matrix holding the 8 different B matrixes corrosponding to each
%node
B = zeros(6, 24, 8);

length_x = d_FE; % [mm]
length_y = d_FE; % [mm]
length_z = d_FE; % [mm]

coordinates = zeros(8,3);
coordinates(1,:) = [-length_x/2 -length_y/2 -length_z/2];
coordinates(2,:) = [length_x/2 -length_y/2 -length_z/2];
coordinates(3,:) = [length_x/2 length_y/2 -length_z/2];
coordinates(4,:) = [-length_x/2 length_y/2 -length_z/2];
coordinates(5,:) = [-length_x/2 -length_y/2 length_z/2];
coordinates(6,:) = [length_x/2 -length_y/2 length_z/2];
coordinates(7,:) = [length_x/2 length_y/2 length_z/2];
coordinates(8,:) = [-length_x/2 length_y/2 length_z/2];

coordinates_gauss = zeros(8,3);
coordinates_gauss(1,:) = [-1/sqrt(3) -1/sqrt(3) -1/sqrt(3)];
coordinates_gauss(2,:) = [1/sqrt(3) -1/sqrt(3) -1/sqrt(3)];
coordinates_gauss(3,:) = [1/sqrt(3) 1/sqrt(3) -1/sqrt(3)];
coordinates_gauss(4,:) = [-1/sqrt(3) 1/sqrt(3) -1/sqrt(3)];
coordinates_gauss(5,:) = [-1/sqrt(3) -1/sqrt(3) 1/sqrt(3)];
coordinates_gauss(6,:) = [1/sqrt(3) -1/sqrt(3) 1/sqrt(3)];
coordinates_gauss(7,:) = [1/sqrt(3) 1/sqrt(3) 1/sqrt(3)];
coordinates_gauss(8,:) = [-1/sqrt(3) 1/sqrt(3) 1/sqrt(3)];

for k = 1:8
    xi1 = coordinates_gauss(k,1);
    xi2 = coordinates_gauss(k,2);
    xi3 = coordinates_gauss(k,3);
    dShape = (1/8)*[-(1-xi2)*(1-xi3),(1-xi2)*(1-xi3),...
                (1+xi2)*(1-xi3),-(1+xi2)*(1-xi3),-(1-xi2)*(1+xi3),...
                (1-xi2)*(1+xi3),(1+xi2)*(1+xi3),-(1+xi2)*(1+xi3);...
                -(1-xi1)*(1-xi3),-(1+xi1)*(1-xi3),(1+xi1)*(1-xi3),...
                (1-xi1)*(1-xi3),-(1-xi1)*(1+xi3),-(1+xi1)*(1+xi3),...
                (1+xi1)*(1+xi3),(1-xi1)*(1+xi3);-(1-xi1)*(1-xi2),...
                -(1+xi1)*(1-xi2),-(1+xi1)*(1+xi2),-(1-xi1)*(1+xi2),...
                (1-xi1)*(1-xi2),(1+xi1)*(1-xi2),(1+xi1)*(1+xi2),...
                (1-xi1)*(1+xi2)];
     JacobianMatrix = dShape*coordinates;
     auxiliar = inv(JacobianMatrix)*dShape;
     for i=1:3
         for j=0:7
             B(i,3*j+1+(i-1), k) = auxiliar(i,j+1);
         end
     end
     % Construct fourth row
     for j=0:7
         B(4,3*j+1, k) = auxiliar(2,j+1);
     end
     for j=0:7
         B(4,3*j+2, k) = auxiliar(1,j+1);
     end
     % Construct fifth row
     for j=0:7
         B(5,3*j+3, k) = auxiliar(2,j+1);
     end
     for j=0:7
         B(5,3*j+2, k) = auxiliar(3,j+1);
     end
     % Construct sixth row
     for j=0:7
         B(6,3*j+1, k) = auxiliar(3,j+1);
     end
     for j=0:7
         B(6,3*j+3, k) = auxiliar(1,j+1);
     end
end

q = reshape(U_L(edofs_L), [24, 1, tot_elements_L]);

for g_point = 1:8
    stress = pagemtimes(pagemtimes(D_comp, B(:, :, g_point)), q); % [mm] * [kN/mm2] * [mm]
    stress_tensor_new = [stress(1, :, :), stress(6, :, :), stress(5, :, :);
                         stress(6, :, :), stress(2, :, :), stress(4, :, :);
                         stress(5, :, :), stress(4, :, :), stress(3, :, :)];
    princ_stress_new = eig3(stress_tensor_new);
    max_comp_stress_L = min(cat(3, max_comp_stress_L, min(princ_stress_new)'),[],3);
    max_tens_stress_L = max(cat(3, max_tens_stress_L, max(princ_stress_new)'),[],3);
    princ_comp_L(:, g_point) = min(princ_stress_new); princ_tens_L(:, g_point) = max(princ_stress_new);
end