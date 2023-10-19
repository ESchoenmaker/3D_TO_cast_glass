function F_plot_stress_new(nelx_L, nely_L, nelz_L, d_FE_L, stress, label, xPhys_LC)

%find gemiddelde displacement of each element by looping over the nodes,
%getting the U in y direction, adding them up together and deviding by
%amount of nodes

maximum = max(max(stress));
minimum = min(min(stress));
%xPhys_LC = permute(xPhys_LC, [2 1 3]);
xPhys_LC = reshape(xPhys_LC, [numel(xPhys_LC), 1]);

face = [1 2 3 4; 2 6 7 3; 6 5 8 7; 5 1 4 8; 1 2 6 5; 4 3 7 8];
element = 1;
for K = 1:nely_L
    y = (K-1)*d_FE_L;
    for I = 1:nelx_L
        x = (I-1)*d_FE_L;
        for J = 1:nelz_L
            z = (J-1)*-1*d_FE_L;
            if xPhys_LC(element, 1) > 0.3
                vert = [x y-d_FE_L z-d_FE_L;...
                     x y z-d_FE_L;...
                     x+d_FE_L y z-d_FE_L;...
                     x+d_FE_L y-d_FE_L z-d_FE_L;...
                     x y-d_FE_L z;...
                     x y z;...
                     x+d_FE_L y z;...
                     x+d_FE_L y-d_FE_L z];
                col = reshape(stress(element, :), [8,1]);
                patch('Faces',face,'Vertices',vert,'FaceVertexCData',col,'FaceColor','interp');
                hold on 
            end
            element = element+1;
        end
    end
end
colormap("jet")
c = colorbar('Ticks',[minimum, 0, maximum]);
c.Label.String = label;






