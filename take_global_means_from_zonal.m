% returns global annual mean of provided zonal mean of 3D variable
%  input 1: data - monthly zonal mean data of 3D atm variable - 192 x 70 x 12y x n
%  input 2: wt - latitudinal weight, a 192x1 vector
% output 1:     global_annual_mean -   y x n
% output 2:  ea_global_annual_mean -   y x 1  (ensemble averaged)
% output 3:    global_monthly_mean - 12y x n
% output 4: ea_global_monthly_mean - 12y x 1  (ensemble averaged)
function [global_annual_mean, ea_global_annual_mean, global_monthly_mean, ea_global_monthly_mean] = take_global_means_from_zonal(monthly_data, wt)
    y = size(monthly_data, 3)/12; %number of years
    n = size(monthly_data, 4); %number of ensemble members

    zonal_mean = squeeze(mean(monthly_data,2));
    global_monthly_mean = sum(zonal_mean.*wt)/sum(wt);
    global_mean_monthly_temp = reshape(global_monthly_mean,[12 y n]);
    global_annual_mean = squeeze(mean(global_mean_monthly_temp));

    ea_global_annual_mean = mean(global_annual_mean, 2);
    ea_global_monthly_mean = mean(global_monthly_mean, 2);
end