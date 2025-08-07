% Returns an ITC timeseries using the precipitation centroid method
%  input 1: PRECT - total precipitation from climate model - nlon x nlat x t x n
%  input 2: lat - latitude vector - lat x 1
% output 1: zonmean_ITCZ - a timeseries of zonal mean ITCZ - t x n
% output 2: ITCZ - a timeseries of ITCZ at each longitude level - nlon x t x n
function [zonmean_ITCZ, ea_zonmean_ITCZ, ITCZ, ea_ITCZ] = funcITCZusingPrecipCentroid(PRECT,lat)
    if ndims(PRECT) < 3 || ndims(PRECT) >4
        error('PRECT must have 3 dimensions (single data set) or 4 dimensions (ensemble).')
    elseif ~isvector(lat)
        error('lat must be a vector.')
    elseif size(PRECT,2) ~= length(lat)
        error('PRECT must have same number of latitude levels as lat')
    end
    nlon = size(PRECT,1); % number of longitude levels in the data
    t = size(PRECT,3); % number of timesteps in the data
    if ndims(PRECT) == 3
        n = 1;
    else
        n = size(PRECT,4); % number of ensemble members in the data
    end
    ITCZ = zeros(nlon,t,n); % preallocate ITCZ array, nlon x t x n
    for k = 1:n
        for j = 1:t
            for i = 1:nlon
                ITCZ(i,j,k) = sum(PRECT(i,lat >= -20 & lat <= 20,j,k)'.*lat(lat >= -20 & lat <= 20).*cosd(lat(lat >= -20 & lat <= 20)))/sum(PRECT(i,lat >= -20 & lat <= 20,j,k)'.*cosd(lat(lat >= -20 & lat <= 20)));
            end % calculate ITCZ for each combination of longitude, timestep, and ensemble member individually
        end
    end
    zonmean_ITCZ = squeeze(mean(ITCZ)); % take zonal mean of ITCZ
    if n == 1
        ea_ITCZ = ITCZ;
        ea_zonmean_ITCZ = zonmean_ITCZ;
    else
        ea_ITCZ = mean(ITCZ,2);
        ea_zonmean_ITCZ = mean(zonmean_ITCZ,2);
    end
end