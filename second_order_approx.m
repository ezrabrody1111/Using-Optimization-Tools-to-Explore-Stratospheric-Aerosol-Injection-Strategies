% returns annually-averaged coefficients of 2nd-order approximation of a variable wrt
% latitude, like T0, T1, and T2
% only works on 2D atmospheric variables
% (y is the number of years and n is the numer of ensemble members)
%  input 1: monthly_data - monthly data of 2D variable (288 x 192 x 12y x n)
%  input 2: lat - latitude, a 192x1 vector
% output 1:     fit012 -   3 x y x n
%                 first row -   zero-order coeff
%                second row -  first-order coeff
%                 third row - second-order coeff
% output 2:  ea_fit012 -   3 x y  (ensemble averaged)
function [fit012, ea_fit012] = second_order_approx(monthly_data, lat)
    y = size(monthly_data, 3)/12; %number of years
    n = size(monthly_data, 4); %number of ensemble members

    zonal_mean = squeeze(mean(monthly_data)); %(192 x 12y x n)
    zonal_annual_mean = reshape(zonal_mean, [192 12 y n]); %(192 x 12 x y x n)
    zonal_annual_mean = squeeze(mean(zonal_annual_mean,2)); %(192 x y x n)

    a0 = sum(zonal_annual_mean.*cosd(lat))/sum(cosd(lat)); %(y x n)
    a1 = sum(zonal_annual_mean.*cosd(lat).*sind(lat))/sum(cosd(lat));
    a2 = sum(.5*(3*(sind(lat)).^2-1).*zonal_annual_mean.*cosd(lat))/sum(cosd(lat));

    fit012 = zeros(3,y,n);
    fit012(1,:,:) = a0;
    fit012(2,:,:) = a1;
    fit012(3,:,:) = a2;

    ea_fit012 = mean(fit012, 3);
end