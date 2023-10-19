function [nodes, dofs] = f_read_nodes_gh(name_excel, sheet_name, d_FE_L, nodx_Lc, nodz_Lc, nody_Lc, mat_dofs)


mat = readmatrix(name_excel, 'Sheet', sheet_name, 'Range', 'A2:C2000');
mat = rmmissing(mat);
coords = [1; numel(nodz_Lc)*numel(nodx_Lc); -1*numel(nodx_Lc)]' .* mat ./ d_FE_L;
nodes = sum(coords, 2) + 1;

dofs = nodes * 3 + mat_dofs;
dofs = reshape(dofs, [numel(dofs), 1]);