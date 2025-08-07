% returns various timeseries of global means of the given monthly data
% only works on 2D atmospheric variables
% (y is the number of years and n is the numer of ensemble members)
%  input 1: monthly_data - monthly data of 2D variable
%  input 2: wt - latitudinal weight, a 192x1 vector
% output 1:     global_annual_mean -   y x n
% output 2:  ea_global_annual_mean -   y x 1  (ensemble averaged)
% output 3:    global_monthly_mean - 12y x n
% output 4: ea_global_monthly_mean - 12y x 1  (ensemble averaged)
function [global_annual_mean, ea_global_annual_mean, global_monthly_mean, ea_global_monthly_mean] = take_global_means(monthly_data, wt)
    y = size(monthly_data, 3)/12; %number of years
    n = size(monthly_data, 4); %number of ensemble members

    zonal_mean = squeeze(mean(monthly_data));
    global_monthly_mean = squeeze(sum(zonal_mean.*wt)/sum(wt));
    if(size(global_monthly_mean,1)) == 1
        global_monthly_mean = global_monthly_mean'; %this is necessary for n=1
    end
    global_mean_monthly_temp = reshape(global_monthly_mean,[12 y n]);
    global_annual_mean = squeeze(mean(global_mean_monthly_temp));
    if(size(global_annual_mean,1)) == 1
        global_annual_mean = global_annual_mean'; %this is necessary for n=1
    end

    ea_global_annual_mean = mean(global_annual_mean, 2);
    ea_global_monthly_mean = mean(global_monthly_mean, 2);
end