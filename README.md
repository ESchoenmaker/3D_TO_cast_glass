# 3D_TO_cast_glass
Three-dimensional topology optimization algorithm for the design of cast glass components.

Functions_used folder contains functions used in the algorithm. 
TO_algorithm folder contains different variations of the topology optimization code. 
There are two versions that completely function in matlab, the boundary conditions here have to be defined in code. The two versions are a volume optimization and a compliance optimization. 
There is a version of the volume optimization that uses a grasshopper script to define the boundary conditions. 

How to use the algorithm in combination with grasshopper:
- create a excel file and copy the link of this in the
- alter the domain and boundary conditions to liking (ensure that the sheet name input is rational)
- enable toggle to write coordinates of nodes to the excel file (file needs to be closed)
- save the excel folder in the same folder as the MATLAB algorithm
- open the matlab script
- change the name_excel to the name of the excel folder containing the boundary conditions in both the main script as the F_data_gh function. 
- input material and optimization specific parameters in F_data_gh, additionally give the run a name
- under boundary conditions in the main script, use the funtion f_read_nodes_gh to read the coordinates out of grasshopper and change them to the right nodal numbers and degrees of freedom corrosponding to the algorithm setup. The last imput 'mat_dofs' refers to in which direction the degree of freedom should be calculated [-2 -1 0] for [y x z] multiple can be kept in case of a node that is fixed for multiple directions. 
- define the load on the individual nodes. The differentiation between load_middle and load_sides ensures that the side nodes have less load applied to them than those carrying a larger surface. 
- ensure that under 'choose constraint combination' the optimization uses all the constraints you want to optimize with
- run the code

 Additional changes output:
The output of the algorithm can be a new excel file, a density matrix and or figure of the geometry. Changes can be made in the F_output_case function to specify the preferred output. Overwriting the excel with the output of each iteration is advised. Not plotting the figure at each iteration can save a lot of time, and is advised for larger optimizaitons. Saving the density plot of each iteration is advised for long calculations. 
