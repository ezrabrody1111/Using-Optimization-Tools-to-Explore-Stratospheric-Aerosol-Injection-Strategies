% returns latitude and timeseries of annually-averaged zonal means of the given monthly data
% only works on 2D atmospheric variables
% (y is the number of years and n is the numer of ensemble members)
%  input 1: monthly_data - monthly data of 2D variable
% output 1:     zonal_annual_mean -   y x n
% output 2:  ea_zonal_annual_mean -   y x 1  (ensemble averaged)
function [zonal_annual_mean, ea_zonal_annual_mean] = take_zonal_mean(monthly_data)
    y = size(monthly_data, 3)/12; %number of years
    n = size(monthly_data, 4); %number of ensemble members

    zonal_mean = squeeze(mean(monthly_data));
    zonal_mean_monthly_temp = reshape(zonal_mean,[192 12 y n]);
    zonal_annual_mean = squeeze(mean(zonal_mean_monthly_temp,2));

    ea_zonal_annual_mean = mean(zonal_annual_mean, 3);
end