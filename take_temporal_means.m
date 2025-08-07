%returns temporal means of 2D variables over a given time frame
% (y is the number of years and n is the numer of ensemble members)
%  input 1: data - monthly (default) or yearly data of 2D variable
%  input 2: data_start - year in which the given data starts
%  input 3: avg_start - year in which the averaging is to start.
%                       Must be >= data_start
%  input 4: avg_end - year in which averaging is to end.
%                     Must be <= data_start + size(monthly_data(3))/12
%  input 5: resolution - OPTIONAL. Temporal resolution of provided data.
%                        'monthly' or 'yearly'. If not specified,
%                        defaults to monthly.
% output 1:    temporal_mean - 289 x 192 x n
% output 2: ea_temporal_mean - 289 x 192      (ensemble averaged)
% output 3: std_temporal_mean - 289 x 192  (standard deviation of the
% temporal mean of each grid cell over the ensemble members)
function [temporal_mean, ea_temporal_mean, std_temporal_mean] = take_temporal_means(data, data_start, avg_start, avg_end, resolution)
    if(nargin < 4 || nargin > 5)
        error('Invalid number of inputs')
    elseif(data_start>avg_start || avg_start>avg_end)
       error('data_start must not be greater than avg_start, and avg_start must not be greater than avg_end')
    elseif(nargin == 4 || (nargin == 5 && strcmp(resolution, 'monthly')))
        data_temp = data(:,:,12*(avg_start-data_start)+1:12*(avg_end-data_start)+12,:);
    elseif(nargin == 5 && strcmp(resolution, 'yearly'))
        data_temp = data(:,:,(avg_start-data_start)+1:(avg_end-data_start)+1,:);
    else
        error('Invalid resolution argument. Must be ''monthly'' or ''yearly''')
    end
    temporal_mean = squeeze(mean(data_temp,3));
    if(size(temporal_mean,1)==288)
        temporal_mean(289,:,:) = temporal_mean(1,:,:);
    end
    ea_temporal_mean = squeeze(mean(temporal_mean,3));
    std_temporal_mean = std(temporal_mean,0,3);
end