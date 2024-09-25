%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%   High-fidelity structured illumination microscopy by point-spread-function engineering    %%                                             
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Reference
% Daniel Sage, Dimiter Prodanov, Jean-Yves Tinevez and Johannes Schindelin, "MIJ: Making Interoperability Between ImageJ and Matlab Possible", ImageJ User & Developer Conference, 24-26 October 2012, Luxembourg.
% http://bigwww.epfl.ch/publications/sage1205.html

close all;clear;clc;
currPath = fileparts(mfilename('fullpath'));    
cd(currPath);                                   
addpath(genpath( './Main_fun'));

%% install mij.jar and ij.jar in matlab
% 1) Copy MIJ(interfacing imagej and matlab)\mij.jar into the java directory of Matlab (e.g for Window Machine 'D:\Program Files\MATLAB\R2014a\java\').
% 2) Copy MIJ(interfacing imagej and matlab)\ij.jar (ImageJ) into the java directory of Matlab. 
javaaddpath 'D:\Software\MATLAB\R2023a\java\jar\mij.jar';  % Note to replace the installation path
javaaddpath 'D:\Software\MATLAB\R2023a\java\jar\ij.jar';   % Note to replace the installation path
MIJ.start;

%% Run 'HiFi-SIM'
HiFiSIM;                                        % 



