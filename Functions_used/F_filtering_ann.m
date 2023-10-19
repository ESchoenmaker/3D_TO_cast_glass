function Neigh_ann_Lc = F_filtering_ann(nelx_Lc, nely_Lc, nelz_Lc, dmax_Lc, d_FE_Lc)
% filtering function based on the TOP3D filtering implimentation
% returns: 
% Neigh_ann -  matrix with information on which elements fall within rmax of each other

rmax_L = ceil(dmax_Lc/d_FE_Lc)/2;
nele = nelx_Lc*nelz_Lc*nely_Lc;
iH = zeros(nele*(2*(ceil(rmax_L)-1)+1)^2,1);                                              %function as in the top3d script
jH = zeros(size(iH));                                                                     % empty matrix where the density of the neighbouring elements will be filled in for maximum constraint
limH = zeros(size(iH));
k = 0;
for k1 = 1:nely_Lc
    for i1 = 1:nelx_Lc
        for j1 = 1:nelz_Lc
            e1 = (k1-1)*nelx_Lc*nelz_Lc + (i1-1)*nelz_Lc+j1;                              %element we are measuring from
            for k2 = max(k1-(ceil(rmax_L)-1),1):min(k1+(ceil(rmax_L)-1),nely_Lc)          %finding within direction z the + and - element you can go in that axis within the rmin boundary
                for i2 = max(i1-(ceil(rmax_L)-1),1):min(i1+(ceil(rmax_L)-1),nelx_Lc)   	  % "
                    for j2 = max(j1-(ceil(rmax_L)-1),1):min(j1+(ceil(rmax_L)-1),nelz_Lc)  % "
                        e2 = (k2-1)*nelx_Lc*nelz_Lc + (i2-1)*nelz_Lc+j2;                  %neighbouring element
                        if sqrt((i1-i2)^2+(j1-j2)^2+(k1-k2)^2) <= rmax_L
                           k = k+1;
                           iH(k) = e1;                                                    %vector containing xi - element that we want to calculate the density for 
                           jH(k) = e2;                                                    %vector containing all the element numbers of elements that fall within the rmin distance
                           limH(k) = 1;
                        end
                    end
                end
            end
        end
    end
end

%trimming arrays
if (all(iH(:) > 0)) == 0
    cut = find(iH == 0);
    iH = iH(1:(cut-1));
    jH = jH(1:(cut-1));
    limH = limH(1:(cut-1));
end

Neigh_ann_Lc = sparse(iH, jH, limH);
