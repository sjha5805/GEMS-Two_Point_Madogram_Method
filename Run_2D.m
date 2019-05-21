clear all;  close all;  clc;

addpath('tools');

NumberOfPoints=10000;

type_sample='Exampe_Rock_Xray_2D_data';
fname=[(type_sample),'.png'];
volume=imread(fname,'png');

S=size(volume);



tic
[Results,VLE_Coord]=Hetero_linescanning2_2D_mex(volume,NumberOfPoints);
toc

save(['Example_Output_2D.mat'],'Results','NumberOfPoints','S');

