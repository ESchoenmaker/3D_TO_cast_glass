function [connectivity_L, edofs_L] = F_connectivity(nelx_L, nely_L, nelz_L, tot_elements_L, nodx_L, nodz_L)
%make matrix size [number_of_elements, 8] that contains that relates the
%global nodal numbers of each element to their local numbering

connectivity_L = zeros(tot_elements_L, 8);
edofs_L = zeros(tot_elements_L*8*3, 1);
count = 1;
for y = 1 : nely_L
    for x = 1 : nelx_L
        for z = 1 : nelz_L
            n1_L = (y-1)*numel(nodx_L)*numel(nodz_L) + (z-1)*numel(nodx_L)+x;
            n2_L = (y-1)*numel(nodx_L)*numel(nodz_L) + z*numel(nodx_L)+x;
            n3_L = y*numel(nodx_L)*numel(nodz_L) + (z-1)*numel(nodx_L)+x;
            n4_L = y*numel(nodx_L)*numel(nodz_L) + z*numel(nodx_L)+x;
            connectivity_L(count, :) = [n4_L, n2_L, n2_L+1, n4_L+1, n3_L, n1_L, n1_L+1, n3_L+1];
            count = count + 1;
        end
    end
end

for i = 1:tot_elements_L
    nodes = connectivity_L(i, :);
    for j = 1:8
        edofs_L((i-1)*24 + j*3-2) = nodes(1, j)*3 -2;
        edofs_L((i-1)*24 + j*3-1) = nodes(1, j)*3 -1;
        edofs_L((i-1)*24 + j*3) = nodes(1, j)*3;
    end
end

end