function [neigh_dis_Lc, neigh_sum_Lc ] = F_filtering(nelx_Lc, nely_Lc, nelz_Lc, dmin_Lc, d_FE_Lc)
% filtering function based on the TOP3D filtering implimentation
% Returns:
% neigh_dis - Matrix containing information on the distances between all elements
% neigh_sum - Matrix containing sum of all distances 


rmin_L = (ceil(dmin_Lc/d_FE_Lc)-1)/2;
nele = nelx_Lc*nelz_Lc*nely_Lc;
iH = ones(nele*(2*(ceil(rmin_L)-1)+1)^2,1);                                                %function as in the top3d script
jH = ones(size(iH));
sH = zeros(size(iH));
k = 0;
for k1 = 1:nely_Lc
    for i1 = 1:nelx_Lc
        for j1 = 1:nelz_Lc
            e1 = (k1-1)*nelx_Lc*nelz_Lc + (i1-1)*nelz_Lc+j1;                                %element we are measuring from
            for k2 = max(k1-(ceil(rmin_L)-1),1):min(k1+(ceil(rmin_L)-1),nely_Lc)            %finding within direction z the + and - element you can go in that axis within the rmin boundary
                for i2 = max(i1-(ceil(rmin_L)-1),1):min(i1+(ceil(rmin_L)-1),nelx_Lc)        % "
                    for j2 = max(j1-(ceil(rmin_L)-1),1):min(j1+(ceil(rmin_L)-1),nelz_Lc)    % "
                        e2 = (k2-1)*nelx_Lc*nelz_Lc + (i2-1)*nelz_Lc+j2;                    %neighbouring element
                        if sqrt((i1-i2)^2+(j1-j2)^2+(k1-k2)^2) <= rmin_L
                            k = k+1;
                            iH(k) = e1;                                                     %vector containing xi - element that we want to calculate the density for 
                            jH(k) = e2;                                                     %vector containing all the element numbers of elements that fall within the rmin distance
                            sH(k) = max(0,rmin_L*d_FE_Lc-sqrt(((i1-i2)*d_FE_Lc)^2+((j1-j2)*d_FE_Lc)^2+((k1-k2)*d_FE_Lc)^2));      %distance between element 1 and 2, if distance is larger sH is smaller
                        end
                    end
                end
            end
        end
    end
end
if (all(iH(:) > 0)) == 0
    cut = find(iH == 0);
    iH = iH(1:(cut-1));
    jH = jH(1:(cut-1));
    sH = sH(1:(cut-1));
end
neigh_dis_Lc = sparse(iH,jH,sH);
neigh_sum_Lc = sum(neigh_dis_Lc,2);  

