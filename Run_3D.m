clear all;  close all;  clc;

addpath('tools');

NumberOfPoints=100000;

fname=['Example_Rock_Xray_Data.mat'];
load(fname);

S=size(volume);

tic
[Results]=Hetero_linescanning2_mex(volume,NumberOfPoints);
toc

save(['Example_Output_3D.mat'],'Results','NumberOfPoints','S');
    

