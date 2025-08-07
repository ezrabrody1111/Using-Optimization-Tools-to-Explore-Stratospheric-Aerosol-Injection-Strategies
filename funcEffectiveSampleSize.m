% calculate effective sample size
% input: time-series data and nominal sample size (n)
function n_eff=funcEffectiveSampleSize(data,n)
    % calculate autocorrelations
    [acf,~]=autocorr(data);
    % AR(1) correlation coefficient \rho_1
    rho1=acf(2);
    % effective sample size
    n_eff=n*(1-rho1)/(1+rho1);
end