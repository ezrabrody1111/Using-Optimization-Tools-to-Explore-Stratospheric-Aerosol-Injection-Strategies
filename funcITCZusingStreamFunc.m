% function for calculating ITCZ position using stream function
% author: Yan Zhang
% date: 08/08/2022
% edited: Ezra Brody 2023

% input 1: V_ann_zon - annual zonal mean V data  192 x 70 x y
% input 2: lat - latitude         192 x 1
% input 3: lev - altitude          70 x 1
%output 1: ITCZ - ITCZ timeseries   y x 1
function ITCZ = funcITCZusingStreamFunc(V_ann_zon,lat,lev)
    y = size(V_ann_zon,3); %number of years
    
    % calculate stream function
    stream=zeros(192,y);
    radius=6378000; % m
    g=9.8; %m/s^2
    [~,lev_idx]=min(abs(lev-500));
    for k=1:y
        for i=1:192
            tmp_int=0;
            p=lev_idx;
            while p>1
                tmp_int=tmp_int+100*(lev(p)-lev(p-1))*mean(V_ann_zon(i,p-1:p,k));% change presure from hPa to Pa
                p=p-1;
            end
            lat_rad=pi*lat(i)/180; % latitude in radian
            stream(i,k)=1e-10*2*pi*radius*cos(lat_rad)*tmp_int/g;
        end
    end
    ITCZ=zeros(y,1);
    idx_all=zeros(y,1);
    lat_idx=96;
    for k=1:y
        if stream(lat_idx,k)<=0
            % find the idx of lat where the streamfunction changes to negative
            idx_all(k)=find(stream(lat_idx+1:192,k)>0,1)+lat_idx-1;
        else
            idx_all(k)=find(stream(1:lat_idx,k)<0,1,'last');
        end
    end
%     alpha_arr=zeros(1,20);
    for k=1:y
        % interpolation version
        % stream(idx_all(k))<0, stream(1+idx_all(k))>0
        % interpolate
        alpha=abs(stream(idx_all(k),k)/(stream(1+idx_all(k),k)-stream(idx_all(k),k)));
%         alpha_arr(k)=alpha;
        % approximate latitude that has streamfunction equal to 0
        ITCZ(k)=alpha*lat(idx_all(k))+(1-alpha)*lat(1+idx_all(k));
%         % previous version (by Ewa)
%         lat_all(k)=0.5*lat(idx_all(k))+0.5*lat(1+idx_all(k));
    end
    % 20-yr average
    %ITCZ_mean=mean(ITCZ);
%     plot(alpha_arr)
%     plot(lat_all)
%     n=18;
%     ITCZ_detrended=detrend(lat_all);
%     n_eff=funcEffectiveSampleSize(ITCZ_detrended,n);
%     ITCZ_std=sqrt(sum(ITCZ_detrended.^2)/n_eff);
end