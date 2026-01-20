%% Nettoyage
    clear 
    close all
    clc
    
%% Addpath
    addpath ('D:\Doctorat\Experience\Matlab-traitement\function');
    addpath ('D:\Doctorat\Experience\Matlab-traitement\function\fct_for_fct');
    addpath ('D:\Doctorat\Experience\Matlab-traitement\function\Visualisation_graph');
    addpath ('D:\Doctorat\Experience\Matlab-traitement');
    addpath ('D:\Doctorat\Experience\Matlab-traitement\fct_mvt');
    addpath ('D:\Doctorat\Experience\Matlab-traitement\Calcule_param');
    addpath (genpath('D:\Doctorat\Experience\Matlab-traitement\lib_try'));
    addpath (genpath('D:\Doctorat\Experience\Matlab-traitement\Art_figures'));
    
%% Import Data    
    %Data = importdata('T_DA_V_SE_Matlab_v2.xlsx');  
    Data = importdata('D_A_IS_SE_Matlab_v1.xlsx'); 
    ST.data = Data.data;
    
%% Variable d'interet    
    [ST.Param] = variable_corr_ST(ST);
    
%% Correlation
    ST.Corr.r.T_V = Corr_All_Sujet( ST.Param.Time , ST.Param.V_IPD);
    ST.Corr.r.T_D = Corr_All_Sujet( ST.Param.Time , ST.Param.Dist);
    ST.Corr.r.T_A = Corr_All_Sujet( ST.Param.Time , ST.Param.Angle);
    ST.Corr.r.V_D = Corr_All_Sujet( ST.Param.V_IPD , ST.Param.Dist);
    ST.Corr.r.V_A = Corr_All_Sujet( ST.Param.V_IPD , ST.Param.Angle);
    
    ST.Corr.r.T_V_D = [ST.Corr.r.T_V, ST.Corr.r.T_D, ST.Corr.r.V_D];
    ST.Corr.r.T_V_A = [ST.Corr.r.T_V, ST.Corr.r.T_A, ST.Corr.r.V_A];

%% Transformation r to Z
    ST.Corr.Z.T_V = r_to_z( ST.Corr.r.T_V );
    ST.Corr.Z.T_D = r_to_z( ST.Corr.r.T_D );
    ST.Corr.Z.T_A = r_to_z( ST.Corr.r.T_A );
    ST.Corr.Z.V_D = r_to_z( ST.Corr.r.V_D );
    ST.Corr.Z.V_A = r_to_z( ST.Corr.r.V_A );
    
    ST.Corr.Z.T_V_D = [ST.Corr.Z.T_V, ST.Corr.Z.T_D, ST.Corr.Z.V_D];
    ST.Corr.Z.T_V_A = [ST.Corr.Z.T_V, ST.Corr.Z.T_A, ST.Corr.Z.V_A];   
    
%% Kmean K = 3  
%     k = 3;
%     idx_1 = kmeans(ST.Corr.Z.T_V_D,k);
%     idx_2 = kmeans(ST.Corr.Z.T_V_A,k);
%     [idx_1, idx_2]
    
    
    