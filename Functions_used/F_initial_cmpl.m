function cmpl = F_initial_cmpl(xPhys_L, nelx_L, nely_L, nelz_L, nodx_L, nody_L, nodz_L, Kel_L, dofs_free_L, pen_L, p_ym_L, p_ymin_L, F_L)
%calculates the initial compliance
%   Detailed explanation goes here
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
    -3*(nelx_L+1) + [0 1 2] -3*((nelx_L+1)*(nelz_L+1)+(nelx_L+1)) + [ 0 1 2 3 4 5] -3*(nelx_L+1) + [3 4 5]],tot_elements_L,1);%this setup is correct (changed the way the nodes are counted)
iK = reshape(kron(edofMat,ones(24,1))',24*24*tot_elements_L,1);
jK = reshape(kron(edofMat,ones(1,24))',24*24*tot_elements_L,1);
sK = reshape(Kel_L(:)*(p_ymin_L+xPhys_L(:)'.^pen_L*(p_ym_L-p_ymin_L)),24*24*tot_elements_L,1); %[mm] * [kN/mm2] = kN/mm
Ktot = sparse(iK,jK,sK); Ktot = (Ktot+Ktot')/2;
U_L = zeros(numel(dofs_all_L),1); 
tolit = 1e-8;
maxit = 8000;
if tot_elements_L >= 10000
M = diag(diag(Ktot(dofs_free_L, dofs_free_L)));
U_L(dofs_free_L, :) = pcg(Ktot(dofs_free_L, dofs_free_L),F_L(dofs_free_L, :), tolit, maxit, M);
else
U_L(dofs_free_L, :) = Ktot(dofs_free_L, dofs_free_L)\ F_L(dofs_free_L, :); % [kN] / [kN/mm] = [mm]
end

%compliance constraint
ce = reshape(sum((U_L(edofMat)*Kel_L).*U_L(edofMat), 2), [nelz_L, nelx_L, nely_L]); % [mm]*[mm]*[mm] gives compliance without taking the physicial(filtered) densities into account
%setGlobalval(ce)
cmpl = sum(sum(sum((p_ymin_L + xPhys_L.^pen_L*(p_ym_L-p_ymin_L)).*ce))); % 

% %assemble total stiffnes matrix
% dofs_all_L = 1:3*numel(nodx_L)*numel(nody_L)*numel(nodz_L);
% Ktot = sparse(numel(dofs_all_L), numel(dofs_all_L));
% for y = 1 : nely_L
%     for z = 1 : nelz_L
%         for x = 1 : nelx_L
%             Exi = p_ymin_L + (xi_L(z,x,y)^pen_L)*(p_ym_L-p_ymin_L);
%             %Exi = Exi*(10^6);
%             n1 = (y-1)*numel(nodx_L)*numel(nodz_L) + (z-1)*numel(nodx_L)+x;
%             n2 = (y-1)*numel(nodx_L)*numel(nodz_L) + z*numel(nodx_L)+x;
%             n3 = y*numel(nodx_L)*numel(nodz_L) + (z-1)*numel(nodx_L)+x;
%             n4 = y*numel(nodx_L)*numel(nodz_L) + z*numel(nodx_L)+x;
%             edof = [n4*3-2; n4*3-1; n4*3; n2*3-2; n2*3-1; n2*3;...
%                 (n2+1)*3-2; (n2+1)*3-1; (n2+1)*3; (n4+1)*3-2; (n4+1)*3-1; (n4+1)*3;...
%                 n3*3-2; n3*3-1; n3*3; n1*3-2; n1*3-1; n1*3;...
%                 (n1+1)*3-2; (n1+1)*3-1; (n1+1)*3; (n3+1)*3-2; (n3+1)*3-1; (n3+1)*3];
%             Kel_xi = Exi*Kel_L;
%             Ktot(edof, edof) = Ktot(edof, edof) + Kel_xi; % - kN/M
%         end
%     end
% end
% 
% %FE analysis 
% U = zeros(numel(dofs_all_L),1);
% U(dofs_free_L, :) = Ktot(dofs_free_L, dofs_free_L)\ F_L(dofs_free_L, :);
% 
% %calculate compliance
% cmpl = 0;
% for y = 1 : nely_L
%     for z = 1 : nelz_L
%         for x = 1 : nelx_L
%             Exi = p_ymin_L + (xi_L(z,x,y)^pen_L)*(p_ym_L-p_ymin_L);
%             %Exi = Exi*(10^6); 
%             n1 = (y-1)*numel(nodx_L)*numel(nodz_L) + (z-1)*numel(nodx_L)+x;
%             n2 = (y-1)*numel(nodx_L)*numel(nodz_L) + z*numel(nodx_L)+x;
%             n3 = y*numel(nodx_L)*numel(nodz_L) + (z-1)*numel(nodx_L)+x;
%             n4 = y*numel(nodx_L)*numel(nodz_L) + z*numel(nodx_L)+x;
%             edof = [n4*3-2; n4*3-1; n4*3; n2*3-2; n2*3-1; n2*3;...
%                 (n2+1)*3-2; (n2+1)*3-1; (n2+1)*3; (n4+1)*3-2; (n4+1)*3-1; (n4+1)*3;...
%                 n3*3-2; n3*3-1; n3*3; n1*3-2; n1*3-1; n1*3;...
%                 (n1+1)*3-2; (n1+1)*3-1; (n1+1)*3; (n3+1)*3-2; (n3+1)*3-1; (n3+1)*3];
%             Uel = U(edof, 1);
%             Kel_xi = Exi*Kel_L;
%             cmpl = cmpl + Uel'*Kel_xi*Uel;
%         end
%     end
% end

end