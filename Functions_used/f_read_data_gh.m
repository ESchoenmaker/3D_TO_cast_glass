function [d_FE, Length, Width, Height] = f_read_data_gh(name_excel, sheet_name)

data = readmatrix(name_excel, 'Sheet',sheet_name, 'Range', 'A2:D2');
d_FE = data(1);
Length = data(2);
Width = data(3);
Height = data(4);