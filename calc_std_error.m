% Calculates the standard error of the supplied variable in each grid box
% over the supplied time interval
% input  - annual_data - nlon x nlat x y x n  or  y x n
% output - std_err     - nlon x nlat
function [std_err] = calc_std_error(annual_data)
%
if ndims(annual_data) > 2
    nlon = size(annual_data,1);
    nlat = size(annual_data,2);
    y = size(annual_data,3); % # of years in the data
    if ndims(annual_data) == 3
        n = 1;
    else
        n = size(annual_data,4);
    end
    
    std_err = zeros(nlon,nlat);
    ann_detrended = zeros(y,n);
    n_eff = zeros(n,1);
    for i =1:nlon
        for j=1:nlat
            for k = 1:n
                ann_detrended(:,k)=detrend(squeeze(annual_data(i,j,:,k)));
                % remove autocorrelation: AR(1)
                n_eff(k)=funcEffectiveSampleSize(ann_detrended(:,k),y);
            end
            std_err(i,j)=sqrt(sum(ann_detrended.^2,"all")/(sum(n_eff-2)*sum(n_eff))); %-2 because 2 years are used to calculate mean and slope
            %std_err(i,j)=std_err(i,j)*sqrt(n); % multiply by sqrt(n) to show...?
        end 
    end
else
    y = size(annual_data,1); % # of years in the data
    if ndims(annual_data) == 1
        n = 1;
    else
        n = size(annual_data,2);
    end

    ann_detrended = zeros(y,n);
    n_eff = zeros(n,1);
    for k = 1:n
        ann_detrended(:,k)=detrend(squeeze(annual_data(:,k)));
        % remove autocorrelation: AR(1)
        n_eff(k)=funcEffectiveSampleSize(ann_detrended(:,k),y);
    end
    std_err=sqrt(sum(ann_detrended.^2,"all")/(sum(n_eff-2)*sum(n_eff)));
    %std_err=std_err*sqrt(n);
end
%
end