function [Corr_sujet] = Corr_All_Sujet( X , Y)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

max_loop = length(X)/120;

incr_1 = 1;
incr_2 = 0;

for i = 1:max_loop
    
    [rho] = corr( X( (incr_1:(incr_2 + 120)), 1), Y( (incr_1:(incr_2 + 120)), 1));
    
                     %loop_param.d = loop_param.d +1;
                     [Corr_sujet(i,1)] = rho;
    
    incr_1 = incr_1 + 120;
    incr_2 = incr_2 + 120;
    
    i = i+1;
end

end

