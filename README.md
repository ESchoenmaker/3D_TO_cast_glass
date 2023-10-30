# 3D_TO_cast_glass
Three-dimensional topology optimization algorithm for the design of cast glass components.

The’Functions used’ folder contains functions used in the algorithm. 
’TO algorithm’ folder contains different variations of the topology optimization code. 
There are two versions that completely function in matlab, the boundary conditions here have to be defined in code. 
The two versions are a volume optimization and a compliance optimization.
For both versions there is also a version that uses a grasshopper script to define the boundary conditions.

**How to use the algorithm in combination with grasshopper:** 

- Create a excel file and copy the path to the panel labeled 'Excel to write to' in the 'Boundary_to_excel.gh' file.
- Alter the domain and boundary conditions to liking using the sliders connected to the python component, the worksheet name will be used to read the coordinates in matlab.
- Enable toggle to write coordinates of nodes to the excel file (file needs to be closed)
- Save the excel folder in the same folder as the MATLAB algorithm
- Open the matlab script
- Change the 'name_excel' variable to the name of the excel folder containing the boundary conditions in the 'F_data_gh' function. Additionally in this file changes can be made to the material properties as well as naming the run.  
- Under boundary conditions in the main script, use the funtion f_read_nodes_gh to read the coordinates out of grasshopper. This function will automatically change them to the right nodal numbers and degrees of freedom corresponding to the algorithm setup. The last input 'mat_dofs' refers to in which direction the degree of freedom should be calculated [-2 -1 0] for [y x z] multiple can be kept in case of a node that is fixed for multiple directions.
- Define the load on the individual nodes. The differentiation between load_middle and load_sides ensures that the side nodes have less load applied to them than those carrying a larger surface. 
- Ensure that under 'choose constraint combination' the optimization uses all the constraints you want to optimize with
- run the code
- When finished matlab should show the density map
- If desired, the density map can be transferred to rhino using the 'Excel_to_Geometry.gh' script

**Additional changes output:**

The output of the algorithm can be a new excel file, a density matrix and or figure of the geometry. Changes can be made in the F_output_case function to specify the preferred output. Overwriting the excel with the output of each iteration is advised. Not plotting the figure at each iteration can save a lot of time, and is advised for larger optimizations. Saving the density plot of each iteration is advised for long calculations.

There are three plotting function integrated in the algorithm. 'F\_plot\_results' will plot exactly the design domain as entered without mirroring. 'F_plot_results_mirror' offers the possibility to mirror the plot in the x and y direction. 'F_plot_stress\_FEA_final' will plot two graphs displaying the tensile and compressive stresses in the final design. 

For all three plot options it is possible to plot either all elements, or the elements from a certain density threshold. This threshold can be set by the user, the standard value is set to plot any element with a density above 0.3. 

**Possible error messages:**

Errors regarding issues using reshape or arrays not having the correct sizes for matrix multiplication can often times be solved by clearing the workspace. 
