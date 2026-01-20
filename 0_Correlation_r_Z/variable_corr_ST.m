function [Param] = variable_corr_ST(ST)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

    Param.Sujet = ST.data(:,1);
    Param.Cible = ST.data(:,4);
    Param.Material = ST.data(:,6);
    Param.Slope = ST.data(:,8);
    Param.Trial = ST.data(:,10);
    Param.Trial_All = ST.data(:,12);
    Param.Time = ST.data(:,13);
    Param.Dist = ST.data(:,14);
    Param.Angle = ST.data(:,15);
    Param.V_IPD = ST.data(:,16);
    Param.SE = ST.data(:,17);

end

