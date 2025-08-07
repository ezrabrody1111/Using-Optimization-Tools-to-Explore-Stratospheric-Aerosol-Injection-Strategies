%% Simulation information
n_base = 3; %Number of baseline (SSP2-4.5) ensemble members
n_def = 3;
n_EQ = 3;
n_15N_15S = 3;
n_30N_30S = 3;
n_60N_60S = 3;
n_30N = 3;
n_60N = 3;

start_base = 2000; %Number of baseline (SSP2-4.5) ensemble members
start_def = 2035;
start_EQ = 2035;
start_15N_15S = 2035;
start_30N_30S = 2035;
start_60N_60S = 2035;
start_30N = 2035;
start_60N = 2035;

SAI_end = 2069; %end year of all simulations (Dec 31)

years = start_def:SAI_end;
base_years = start_base:SAI_end;
months = 1:12*length(years);
big_months = 1:(1980+12*length(years));

y = length(years);

comp_start = 2050; %year to start comparing steady-state behavior
comp_end   = 2069;
ref_start = 2008; % reference period that has the same GMST as the target in the SAI simulations
ref_end = 2027;

% for plotting the different strategies/scenarios
colors=[0 0 0;
        0.8500 0.3250 0.0980;
        0.9290, 0.6940, 0.1250;
        0.4660 0.6740 0.1880;
        0 0.4470 0.7410;
        0 1 1;
        0.7176 0.2745 1.0000;
        0 1 0];

%% load data
% PRECT and QFLX have already been converted to mm/day!
load('GAUSS.mat');
load('GAUSS_V_zm.mat');
scenarios = {'base','EQ','15N_15S','30N_30S','60N_60S','60N','def','30N'};

%% calculate standard error and standard deviation due to natural variability

PE_se = calc_std_error(squeeze(mean(reshape(PRECT_base(:,:,97:336,:),288,192,12,20,[]),3))-squeeze(mean(reshape(QFLX_base(:,:,97:336,:),288,192,12,20,[]),3)));
P_se = calc_std_error(squeeze(mean(reshape(PRECT_base(:,:,97:336,:),288,192,12,20,[]),3)));
T_se = calc_std_error(squeeze(mean(reshape(TREFHT_base(:,:,97:336,:),288,192,12,20,[]),3)));

PE_std = calc_std_dev(squeeze(mean(reshape(PRECT_base(:,:,97:336,:),288,192,12,20,[]),3))-squeeze(mean(reshape(QFLX_base(:,:,97:336,:),288,192,12,20,[]),3)));
P_std = calc_std_dev(squeeze(mean(reshape(PRECT_base(:,:,97:336,:),288,192,12,20,[]),3)));
T_std = calc_std_dev(squeeze(mean(reshape(TREFHT_base(:,:,97:336,:),288,192,12,20,[]),3)));


%% take global means
varnames = {'PRECT', 'QFLX', 'TREFHT'};
for i = 1:length(scenarios)
    for k=1:length(varnames)
        eval(['[global_annual_mean_' varnames{k} '_' scenarios{i} ', ea_global_annual_mean_' varnames{k} '_' scenarios{i} ',~,~] = take_global_means(' varnames{k} '_' scenarios{i} ', wt);']);
    end
end

%calculate deltas for global mean variables
%delta is defined as the change in a variable due to GHG warming in a
%certain scenario
delta_T0 = mean(global_annual_mean_TREFHT_base(comp_start-start_base+1:comp_end-start_base+1,:)) - mean(global_annual_mean_TREFHT_base(ref_start-start_base+1:ref_end-start_base+1,:));
delta_P0 = mean(global_annual_mean_PRECT_base(comp_start-start_base+1:comp_end-start_base+1,:)) - mean(global_annual_mean_PRECT_base(ref_start-start_base+1:ref_end-start_base+1,:));
P0_std_err = calc_std_error(global_annual_mean_PRECT_def(ref_start-start_base+1:ref_end-start_base+1,:));

%% take temporal means
varnames = {'PRECT', 'QFLX', 'TREFHT'}; % Edit this line to choose which variable to process

for i = 1:length(scenarios)
    for k=1:length(varnames)
        eval(['[temporal_mean_' varnames{k} '_' scenarios{i} ', ea_temporal_mean_' varnames{k} '_' scenarios{i} ',~] = take_temporal_means(' varnames{k} '_' scenarios{i} ', start_' scenarios{i} ', comp_start, comp_end);']);
    end
end

%for reference period
for k=1:length(varnames)
    eval(['[temporal_mean_' varnames{k} '_ref, ea_temporal_mean_' varnames{k} '_ref,~] = take_temporal_means(' varnames{k} '_base, start_base, ref_start, ref_end);']);
end

%calculate deltas for 2D variables
%delta is defined as the change in a variable due to GHG warming in a
%certain scenario
for k=1:length(varnames)
    eval(['delta_' varnames{k} ' = temporal_mean_' varnames{k} '_base - temporal_mean_' varnames{k} '_ref;'])
end

%% Get SSI
SSI_comp_start = 2030; % choose an earlier comparison period for baseline scenario because SSI disappears and thus becomes nonlinear later on
SSI_comp_end = 2049;
SSI_delta_T0 = mean(global_annual_mean_TREFHT_base(SSI_comp_start-start_base+1:SSI_comp_end-start_base+1,:)) - mean(global_annual_mean_TREFHT_base(ref_start-start_base+1:ref_end-start_base+1,:));

for i=1:length(scenarios)
    eval(['S_ICEFRAC_' scenarios{i} '= ICEFRAC_' scenarios{i} '(:,:,9:12:size(ICEFRAC_' scenarios{i} ',3),:);']);
    eval(['SSI_' scenarios{i} '= S_ICEFRAC_' scenarios{i} '(:,97:192,:,:).*AREA(:,97:192);']);
    if i>1
        eval(['total_arctic_SSI_' scenarios{i} ' = squeeze(sum(SSI_' scenarios{i} ',[1 2]));']);
        eval(['ea_total_arctic_SSI_' scenarios{i} '_combined = cat(1,mean(total_arctic_SSI_base(1:35,:),2),mean(total_arctic_SSI_' scenarios{i} ',2));'])   %2000-2069
        eval(['ea_temporal_mean_SSI_' scenarios{i} ' = mean(ea_total_arctic_SSI_' scenarios{i} '_combined(comp_start-start_base+1:comp_end-start_base+1));'])
    else
        total_arctic_SSI_base = squeeze(sum(SSI_base,[1 2]));
        ea_total_arctic_SSI_base = squeeze(mean(total_arctic_SSI_base,2));
        ea_temporal_mean_SSI_base = mean(ea_total_arctic_SSI_base(SSI_comp_start-start_base+1:SSI_comp_end-start_base+1));
    end
end
ea_temporal_mean_SSI_ref = mean(ea_total_arctic_SSI_base(ref_start-start_base+1:ref_end-start_base+1));

%% T1 and T2
for k = 1:length(scenarios)
    eval(['[T012_' scenarios{k} ', ea_T012_' scenarios{k} '] = second_order_approx(TREFHT_' scenarios{k} ', lat);'])
    eval(['temporal_mean_T1_' scenarios{k} ' = squeeze(mean(T012_' scenarios{k} '(2,comp_start-start_' scenarios{k} '+1:comp_end-start_' scenarios{k} '+1,:),2))'';'])
    eval(['temporal_mean_T2_' scenarios{k} ' = squeeze(mean(T012_' scenarios{k} '(3,comp_start-start_' scenarios{k} '+1:comp_end-start_' scenarios{k} '+1,:),2))'';'])
end
temporal_mean_T1_ref = squeeze(mean(T012_base(2,ref_start-start_base+1:ref_end-start_base+1,:),2))';
temporal_mean_T2_ref = squeeze(mean(T012_base(3,ref_start-start_base+1:ref_end-start_base+1,:),2))';

%% ITCZ (streamfunction)
for k = 1:4
    eval(['n = n_' scenarios{k} ';'])
    for i = 1:n
        eval(['ITCZ_' scenarios{k} '(:,i) = funcITCZusingStreamFunc(squeeze(mean(reshape(V_zm_' scenarios{k} '(:,:,:,i),192,70,12,[]),3)),lat,lev);'])
    end
end
for k = 5:8
    eval(['n = n_' scenarios{k} ';'])
    eval(['ITCZ_' scenarios{k} ' = zeros(35,3);'])
    for i = 1:n
        eval(['ITCZ_' scenarios{k} '(11:35,i) = funcITCZusingStreamFunc(squeeze(mean(reshape(V_zm_' scenarios{k} '(:,:,10*12+1:35*12,i),192,70,12,[]),3)),lat,lev);'])
    end
end

for k = 1:length(scenarios)
    eval(['temporal_mean_ITCZ_' scenarios{k} ' = squeeze(mean(ITCZ_' scenarios{k} '(comp_start-start_' scenarios{k} '+1:comp_end-start_' scenarios{k} '+1,:)));'])
end
temporal_mean_ITCZ_ref = squeeze(mean(ITCZ_base(ref_start-start_base+1:ref_end-start_base+1,:)));

%% calculate std error for scalar vars
SSI_std_err = calc_std_error(total_arctic_SSI_base(ref_start-start_base+1:ref_end-start_base+1,:))
T1_std_err = calc_std_error(squeeze(T012_base(2,ref_start-start_base+1:ref_end-start_base+1,:)))
T2_std_err = calc_std_error(squeeze(T012_base(3,ref_start-start_base+1:ref_end-start_base+1,:)))
ITCZ_std_err = calc_std_error(ITCZ_base(ref_start-start_base+1:ref_end-start_base+1,:))

%save('std_err.mat', "ITCZ_std_err", "T2_std_err", "T1_std_err", "SSI_std_err", "T_se", "P0_std_err", "P_se", "PE_se", "-append")

%% calculate alphas
%alpha is defined as the change in a variable due to GHG warming
%per degree of GHG warming
for k = 1:length(varnames)
    eval(['alpha_' varnames{k} ' = zeros(289,192,n_base);'])
    for i = 1:n_base
        eval(['alpha_' varnames{k} '(:,:,i) = (temporal_mean_' varnames{k} '_base(:,:,i) - temporal_mean_' varnames{k} '_ref(:,:,i))/delta_T0(i);'])
        eval(['ea_alpha_' varnames{k} ' = mean(alpha_' varnames{k} ',3);'])
    end
end

alpha_SSI = (ea_temporal_mean_SSI_base-ea_temporal_mean_SSI_ref)/mean(SSI_delta_T0);
alpha_T1 = (temporal_mean_T1_base-temporal_mean_T1_ref)./delta_T0;
alpha_T2 = (temporal_mean_T2_base-temporal_mean_T2_ref)./delta_T0;
alpha_ITCZ = (temporal_mean_ITCZ_base-temporal_mean_ITCZ_ref)./delta_T0;
alpha_P0 = (mean(global_annual_mean_PRECT_base(comp_start-start_base+1:comp_end-start_base+1,:))-mean(global_annual_mean_PRECT_base(ref_start-start_base+1:ref_end-start_base+1,:)))./delta_T0;
ea_alpha_T1 = mean(alpha_T1);
ea_alpha_T2 = mean(alpha_T2);
ea_alpha_ITCZ = mean(alpha_ITCZ);
ea_alpha_P0 = mean(alpha_P0);
%save("alpha.mat", "ea_alpha_P0", "ea_alpha_ITCZ", "ea_alpha_T2", "ea_alpha_T1", "alpha_SSI", "ea_alpha_TREFHT", "ea_alpha_PRECT", "ea_alpha_QFLX")
%% calculate mus
%mu is defined as the change in a variable due to SAI cooling
%per degree of SAI cooling
for k = 1:length(varnames)
    for j = 2:length(scenarios)
        eval(['n = n_' scenarios{j} ';'])
        eval(['cooling_' scenarios{j} ' = zeros(n,1);'])
        eval(['mu_' varnames{k} '_' scenarios{j} ' = zeros(289,192,n);'])
        for i = 1:n
            eval(['cooling_' scenarios{j} '(i) = mean(global_annual_mean_TREFHT_base(comp_start-start_base+1:comp_end-start_base+1,i)) - mean(global_annual_mean_TREFHT_' scenarios{j} '(comp_start-start_' scenarios{j} '+1:comp_end-start_' scenarios{j} '+1,i));'])
            eval(['mu_' varnames{k} '_' scenarios{j} '(:,:,i) = (temporal_mean_' varnames{k} '_' scenarios{j} '(:,:,i) - temporal_mean_' varnames{k} '_base(:,:,i))/cooling_' scenarios{j} '(i);'])
            eval(['ea_mu_' varnames{k} '_' scenarios{j} ' = mean(mu_' varnames{k} '_' scenarios{j} ',3);'])
            %eval(['save("mu.mat", "ea_mu_' varnames{k} '_' scenarios{j} '", "-append")'])
        end
    end
end
for j = 2:length(scenarios)
eval(['n = n_' scenarios{j} ';'])
    eval(['cooling_' scenarios{j} ' = zeros(n,1);'])
    eval(['mu_T1_' scenarios{j} ' = zeros(n,1);'])
    eval(['mu_T2_' scenarios{j} ' = zeros(n,1);'])
    eval(['mu_ITCZ_' scenarios{j} ' = zeros(n,1);'])
    eval(['mu_P0_' scenarios{j} ' = zeros(n,1);'])
    for i = 1:n
        eval(['cooling_' scenarios{j} '(i) = mean(global_annual_mean_TREFHT_base(comp_start-start_base+1:comp_end-start_base+1,i)) - mean(global_annual_mean_TREFHT_' scenarios{j} '(comp_start-start_' scenarios{j} '+1:comp_end-start_' scenarios{j} '+1,i));'])
        eval(['mu_T1_' scenarios{j} '(i) = (temporal_mean_T1_' scenarios{j} '(i) - temporal_mean_T1_base(i))/cooling_' scenarios{j} '(i);'])
        eval(['mu_T2_' scenarios{j} '(i) = (temporal_mean_T2_' scenarios{j} '(i) - temporal_mean_T2_base(i))/cooling_' scenarios{j} '(i);'])
        eval(['mu_ITCZ_' scenarios{j} '(i) = (temporal_mean_ITCZ_' scenarios{j} '(i) - temporal_mean_ITCZ_base(i))/cooling_' scenarios{j} '(i);'])
        eval(['mu_P0_' scenarios{j} '(i) = (mean(global_annual_mean_PRECT_' scenarios{j} '(comp_start-start_' scenarios{j} '+1:comp_end-start_' scenarios{j} '+1,i))-mean(global_annual_mean_PRECT_base(comp_start-start_base+1:comp_end-start_base+1,i)))/cooling_' scenarios{j} '(i);'])
    end
    eval(['ea_mu_SSI_' scenarios{j} ' = ((ea_temporal_mean_SSI_' scenarios{j} ' - ea_temporal_mean_SSI_ref) - alpha_SSI*mean(delta_T0))/mean(cooling_' scenarios{j} ');'])
    eval(['mu_T1_arr(j-1) = mean(mu_T1_' scenarios{j} ');'])
    eval(['mu_T2_arr(j-1) = mean(mu_T2_' scenarios{j} ');'])
    eval(['mu_ITCZ_arr(j-1) = mean(mu_ITCZ_' scenarios{j} ');'])
    eval(['mu_P0_arr(j-1) = mean(mu_P0_' scenarios{j} ');'])
    eval(['mu_SSI_arr(j-1) = ea_mu_SSI_' scenarios{j} ';'])
end
load('mu_predicted.mat')
%save('mu.mat', "mu_T2_arr", "mu_T1_arr", "mu_ITCZ_arr", "mu_SSI_arr", "mu_P0_arr", "-append")
%save('mu.mat', "cooling_EQ", "cooling_15N_15S", "cooling_30N_30S", "cooling_60N_60S", "cooling_60N", "cooling_def", "cooling_30N", "-append")

%% get injection rates
%inj_rate_arr is 35(years) by 5(strategies)
% the order of strategies is:
% 0N_1.0, 15N_15S, 30N_30S, 60N_60S_1.0, Global+1C
load('yearly_injection_rate.mat')

%30S, 15S, 15N, 30N
injections_def_001 = [0	0.602	0.602	0;
0.10002879	1.08381145	1.19767915	0.15063666;
0.97254312	0.9706783	1.06703946	1.01537031;
1.71057871	0.92421212	1.18659987	1.82719549;
0.92702809	2.44699299	2.23098342	0.83102384;
1.36615414	2.94916733	2.696849	1.25401266;
0.25477223	4.43297164	4.16264958	0.1346291;
1.3906934	4.23386891	3.35237533	0.99891848;
0.96693728	5.1052305	3.91710783	0.43888276;
2.07152038	4.98148196	3.32617895	1.33583016;
1.81334111	6.1968299	4.09054657	0.87721518;
1.06315715	7.21680646	4.82470287	0;
1.69575573	6.74951196	4.71712437	0.79247236;
1.49254054	7.55183112	4.97226662	0.34606743;
1.32094568	8.14737202	5.17524423	0;
1.25905812	8.01631559	5.18343483	0;
1.42929716	8.36269676	5.14677816	0;
1.52727699	8.59617595	5.15980273	0;
1.42059678	8.39867815	5.2023354	0;
1.59795402	8.47173814	4.87634161	0;
1.62420636	8.65059089	4.99612659	0;
1.81975184	9.09865523	5.00421359	0;
1.58880441	8.96030355	5.38549362	0;
1.59907653	9.09069513	5.49277293	0;
1.73547488	9.16011083	5.25529235	0;
1.53100956	9.39335413	5.94858263	0;
1.56582884	9.45587915	5.93276425	0;
1.95027785	9.01865085	5.63425968	0.446104;
1.73964936	9.85633586	5.9421248	0;
1.72489942	10.19873114	6.31770745	0;
1.82009081	10.35751113	6.2623068	0;
1.95737456	10.56733309	6.16324033	0;
1.92597264	10.64490365	6.3114652	0;
1.94336039	10.74405499	6.37149412	0;
1.97561914	10.7681852	6.32304214	0];

injections_def_002 = [0	0.602	0.602	0;
0.00765396	1.19994001	1.18271859	0;
0.07052848	1.76667419	1.60798512	0;
0.0710203	2.23909432	2.07929865	0;
0.12368635	2.80824364	2.52994936	0;
0.33868513	3.60657815	2.84453661	0;
0.42979577	4.30985446	3.34281397	0;
0.42804177	4.73632746	3.77323348	0;
0.85856472	5.69355767	3.76178706	0;
0.70995839	6.17311487	4.57570849	0;
0.74065348	6.65868421	4.99221389	0;
0.79820958	6.86461364	5.06864208	0;
0.98945483	7.27625464	5.04998127	0;
0.91265661	7.30804662	5.25456925	0;
0.90926635	7.61049399	5.56464469	0;
1.06869819	7.89902063	5.49444969	0;
1.3383247	8.26832036	5.25708979	0;
1.17464724	8.22839077	5.58543448	0;
1.24495535	8.42602778	5.62487824	0;
1.49655364	8.8867881	5.51954242	0;
1.69418578	8.97484628	5.16292828	0;
1.67753223	9.01474069	5.24029317	0;
1.69948225	9.25703446	5.4331994	0;
1.6072268	9.28058944	5.66432913	0;
1.88307752	9.46773285	5.23080843	0;
1.88331054	9.69657008	5.45912137	0;
2.02904116	9.90587418	5.34053156	0;
1.92681423	9.70986086	5.37452884	0;
2.02223177	9.88881927	5.33879778	0;
2.1029179	10.01691144	5.28534615	0;
2.24868339	10.26309946	5.20356184	0;
2.37833691	10.38710113	5.03584309	0;
2.51174332	10.30643428	4.65501183	0;
2.77591967	10.65810598	4.41228672	0;
2.68514298	10.91969206	4.87812035	0];

injections_def_003 = [0	0.602	0.602	0;
0	1.11215765	1.22263609	0.04910153;
0	1.60114471	1.71316347	0.04978612;
0	2.19607756	2.2667932	0.03142917;
0.76918435	2.26178958	2.42549111	0.84194059;
0.43146276	3.49939678	3.3462087	0.36337917;
1.3875995	3.38615949	3.23083702	1.31856729;
2.08149595	3.6171052	3.10471346	1.85376629;
1.40138668	4.94143553	4.23409782	1.08701436;
1.99256867	5.13513096	4.25014408	1.59924116;
0.53646404	7.09364906	5.88660498	0;
0.55013692	7.41622424	6.17841616	0;
0.65784011	7.48226793	6.00212768	0;
1.61044478	7.05528093	5.35555585	0.85501141;
2.38983488	6.94317353	4.75915892	1.41916173;
3.97891118	5.94176501	3.72595255	2.99410564;
3.48349762	6.50373651	4.19932639	2.45931534;
4.36866325	5.95593563	3.91095613	3.45978347;
3.68388677	6.69773552	4.39268807	2.65942124;
3.83380107	6.59853807	4.24718285	2.78875431;
4.04242878	6.741613	3.9698319	2.81052607;
3.97726411	7.06527235	4.21930699	2.71239062;
4.2336638	7.02267227	4.07949478	2.92558491;
5.77144431	5.69018379	2.89363895	4.5285355;
4.37205932	6.69733224	4.13427131	3.23292113;
3.30767236	8.13487445	4.8950426	1.86774709;
4.11734615	7.30348019	4.4240668	2.83760687;
4.0960642	7.23711262	4.39397093	2.83244568;
2.45655567	8.71844289	5.4631284	1.00974923;
2.13384171	9.30967238	5.76322293	0.55764196;
2.60184815	8.97160577	5.18387646	0.9184129;
2.79131696	8.89434688	5.00595173	1.06314133;
3.56784392	8.31189648	4.58366977	1.91085427;
2.08637772	9.86186246	5.63872816	0.20942914;
1.76654249	10.22209748	6.24737688	0];

% Arctic High (60N) (2035-2071)
% injection_Arctic_high = [9.1782; 10.49668034; 11.0667278718; 11.0164038371; 11.0626947726; 11.2648447319; 11.0949464367; 11.3684169759;...
%     11.3997395915; 11.3095318408; 11.1778255661; 11.1573853842; 11.1734064663; 11.6374206082; 11.5546664967; 11.4426735541; 11.0355683292;...
%     10.6229658139; 10.301512415; 10.0687804386; 10.0078037632; 9.45483087301; 9.04461198403; 9.05299517759; 9.32301248653; 9.82637477143; ...
%     9.95212912789; 9.78926432157; 10.1539789566; 10.2887091035; 10.7002643962; 11.4796016888; 12.2556812013; 13.0135932335; 13.3354220392; ...
%     13.1184875126; 13.1069430668];

ea_injections_def = (injections_def_001 + injections_def_002 + injections_def_003)/3;
ea_inj_def_end = mean(ea_injections_def(16:35,:));
ea_inj_EQ_end = mean(inj_rate_arr(16:35,1));
ea_inj_15N_15S_end = mean(inj_rate_arr(16:35,2));
ea_inj_30N_30S_end = mean(inj_rate_arr(16:35,3));
ea_inj_60N_60S_end = mean(inj_rate_arr(16:35,4));
ea_inj_60N_end = 12;

ControlLog_EQ_001 = readmatrix(strcat('ControlLog_b.e21.BWSSP245.f09_g17.release-cesm2.1.3.WACCM-MA-1deg.feedback_0N_1.0.001.txt'));
ControlLog_EQ_002 = readmatrix(strcat('ControlLog_b.e21.BWSSP245.f09_g17.release-cesm2.1.3.WACCM-MA-1deg.feedback_0N_1.0.002.txt'));
ControlLog_EQ_003 = readmatrix(strcat('ControlLog_b.e21.BWSSP245.f09_g17.release-cesm2.1.3.WACCM-MA-1deg.feedback_0N_1.0.003.txt'));
ControlLog_15N_15S_001 = readmatrix(strcat('ControlLog_b.e21.BWSSP245.f09_g17.release-cesm2.1.3.WACCM-MA-1deg.feedback_15N_15S.001.txt'));
ControlLog_15N_15S_002 = readmatrix(strcat('ControlLog_b.e21.BWSSP245.f09_g17.release-cesm2.1.3.WACCM-MA-1deg.feedback_15N_15S.002.txt'));
ControlLog_15N_15S_003 = readmatrix(strcat('ControlLog_b.e21.BWSSP245.f09_g17.release-cesm2.1.3.WACCM-MA-1deg.feedback_15N_15S.003.txt'));
ControlLog_30N_30S_001 = readmatrix(strcat('ControlLog_b.e21.BWSSP245.f09_g17.release-cesm2.1.3.WACCM-MA-1deg.feedback_30N_30S.001.txt'));
ControlLog_30N_30S_002 = readmatrix(strcat('ControlLog_b.e21.BWSSP245.f09_g17.release-cesm2.1.3.WACCM-MA-1deg.feedback_30N_30S.002.txt'));
ControlLog_30N_30S_003 = readmatrix(strcat('ControlLog_b.e21.BWSSP245.f09_g17.release-cesm2.1.3.WACCM-MA-1deg.feedback_30N_30S.003.txt'));
ControlLog_60N_60S_001 = readmatrix(strcat('ControlLog_b.e21.BWSSP245.f09_g17.release-cesm2.1.3.WACCM-MA-1deg.feedback_60N_60S_1.0.001.txt'));
ControlLog_60N_60S_002 = readmatrix(strcat('ControlLog_b.e21.BWSSP245.f09_g17.release-cesm2.1.3.WACCM-MA-1deg.feedback_60N_60S_1.0.002.txt'));
ControlLog_60N_60S_003 = readmatrix(strcat('ControlLog_b.e21.BWSSP245.f09_g17.release-cesm2.1.3.WACCM-MA-1deg.feedback_60N_60S_1.0.003.txt'));

inj_def_end = [mean(injections_def_001(16:35,:)); mean(injections_def_002(16:35,:)); mean(injections_def_003(16:35,:))];
inj_EQ_end = [mean(ControlLog_EQ_001(16:35,9)); mean(ControlLog_EQ_002(16:35,9)); mean(ControlLog_EQ_003(16:35,9))];
inj_15N_15S_end = 2*[mean(ControlLog_15N_15S_001(16:35,9)); mean(ControlLog_15N_15S_002(16:35,9)); mean(ControlLog_15N_15S_003(16:35,9))];
inj_30N_30S_end = 2*[mean(ControlLog_30N_30S_001(16:35,8)); mean(ControlLog_30N_30S_002(16:35,8)); mean(ControlLog_30N_30S_003(16:35,8))];
inj_60N_60S_end = 2*[mean(ControlLog_60N_60S_001(16:35,8)); mean(ControlLog_60N_60S_002(16:35,8)); mean(ControlLog_60N_60S_003(16:35,8))];

%% plot global mean TREFHT timeseries
figure('Renderer', 'painters', 'Position', [100 200 1000 450])
clf
ax = axes;
hold on
for i = 1:n_base
    plot(.5+base_years,global_annual_mean_TREFHT_base(:,i),':','Linewidth',1,'HandleVisibility','off','Color',colors(1,:));
end
for j = 2:length(scenarios)
    eval(['n = n_' scenarios{j} ';'])
    for i = 1:n
        eval(['plot(.5+years,global_annual_mean_TREFHT_' scenarios{j} '(:,i),'':'',''Linewidth'',1,''HandleVisibility'',''off'',''Color'',colors(j,:));'])
    end
end
plot(.5+base_years,ea_global_annual_mean_TREFHT_base,'k','Linewidth',3);
for j = 2:length(scenarios)
    eval(['n = n_' scenarios{j} ';'])
    eval(['plot(.5+years,ea_global_annual_mean_TREFHT_' scenarios{j} ',''LineWidth'',n,''Color'',colors(j,:))'])
end
hold off
set(gca,'FontSize',14);
axis([2000 2069 287.5 290]);
title('Global Mean Surface Temperature');
ylabel('Temperature (K)')
xlabel('Year');
legend('SSP2-4.5', 'Multi-Objective','EQ','15N-15S','30N-30S','60N-60S','30N','60N','Location','eastoutside');
grid on

%% plot TREFHT and PRECT anomaly maps for multi-objective
Cmin = -3;
Cmax = 3;
[contours, color_map] = custom_color_map(Cmin,Cmax,25,'*RdBu');

% TREFHT error check
Z=abs((ea_temporal_mean_TREFHT_def-ea_temporal_mean_TREFHT_ref)./(2*[T_se; T_se(1,:)]));
[lon_null,lat_null]=find(Z<1);
lon_null=(lon_null-1)/287*360;
lat_null=(lat_null-1)/191*180-90;

figure('Renderer', 'painters', 'Position', [60 200 1200 400])
clf
ax1=subplot(1,2,1); 
hold on
worldmap([-90 90],[-180 180])
setm(ax1, 'meridianlabel', 'off', 'parallellabel', 'off')
box on
contourfm(lat,lon2,ea_temporal_mean_TREFHT_def'-ea_temporal_mean_TREFHT_ref',contours,'LineStyle','none')
scatterm(lat_null,lon_null,.2,'k','filled');
load coastlines
geoshow(coastlat,coastlon,'Color','k')
caxis([Cmin Cmax])
cb=colorbar;
ylabel(cb,'\DeltaK')
hold off
title('\DeltaT from Reference Period 2050-2069 Multi-Objective');
colormap(ax1, color_map);

Cmin = -3;
Cmax = 3;
[contours, color_map] = custom_color_map(Cmin,Cmax,11,'BrBG');

% TREFHT error check
Z=abs((ea_temporal_mean_PRECT_def-ea_temporal_mean_PRECT_ref)./(2*[P_se; P_se(1,:)]));
[lon_null,lat_null]=find(Z<1);
lon_null=(lon_null-1)/287*360;
lat_null=(lat_null-1)/191*180-90;

ax2=subplot(1,2,2);
hold on
worldmap([-90 90],[-180 180])
setm(ax2, 'meridianlabel', 'off', 'parallellabel', 'off')
box on
contourfm(lat,lon2,double(ea_temporal_mean_PRECT_def')-ea_temporal_mean_PRECT_ref',contours,'LineStyle','none')
scatterm(lat_null,lon_null,.2,'k','filled');
load coastlines
geoshow(coastlat,coastlon,'Color','k')
caxis([Cmin Cmax])
cb=colorbar;
ylabel(cb,'\Deltamm/day')
hold off
title('\DeltaP from Reference Period 2050-2069 Multi-Objective');
colormap(ax2, color_map);

%% plot TREFHT and PRECT anomaly maps for multi-objective and baseline
Cmin = -3;
Cmax = 3;
[T_contours, T_color_map] = custom_color_map(Cmin,Cmax,25,'*RdBu');

figure('Renderer', 'painters', 'Position', [80 100 1200 700])
clf

% multi-objective TREFHT error check
Z=abs((ea_temporal_mean_TREFHT_def-ea_temporal_mean_TREFHT_ref)./(2*[T_se; T_se(1,:)]));
[lon_null,lat_null]=find(Z<1);
lon_null=(lon_null-1)/287*360;
lat_null=(lat_null-1)/191*180-90;

ax1=subplot(2,2,1); 
hold on
worldmap([-90 90],[-180 180])
setm(ax1, 'meridianlabel', 'off', 'parallellabel', 'off')
box on
contourfm(lat,lon2,ea_temporal_mean_TREFHT_def'-ea_temporal_mean_TREFHT_ref',T_contours,'LineStyle','none')
scatterm(lat_null,lon_null,.2,'k','filled');
load coastlines
geoshow(coastlat,coastlon,'Color','k')
caxis([Cmin Cmax])
cb=colorbar;
ylabel(cb,'\DeltaK')
hold off
title('\DeltaT from Reference Period 2050-2069 Multi-Objective');
colormap(ax1, T_color_map);

Cmin = -3;
Cmax = 3;
[P_contours, P_color_map] = custom_color_map(Cmin,Cmax,11,'BrBG');

% multi-objective PREFHT error check
Z=abs((ea_temporal_mean_PRECT_def-ea_temporal_mean_PRECT_ref)./(2*[P_se; P_se(1,:)]));
[lon_null,lat_null]=find(Z<1);
lon_null=(lon_null-1)/287*360;
lat_null=(lat_null-1)/191*180-90;

ax2=subplot(2,2,2);
hold on
worldmap([-90 90],[-180 180])
setm(ax2, 'meridianlabel', 'off', 'parallellabel', 'off')
box on
contourfm(lat,lon2,double(ea_temporal_mean_PRECT_def')-ea_temporal_mean_PRECT_ref',P_contours,'LineStyle','none')
scatterm(lat_null,lon_null,.2,'k','filled');
load coastlines
geoshow(coastlat,coastlon,'Color','k')
caxis([Cmin Cmax])
cb=colorbar;
ylabel(cb,'\Deltamm/day')
hold off
title('\DeltaP from Reference Period 2050-2069 Multi-Objective');
colormap(ax2, P_color_map);

% baseline TREFHT error check
Z=abs((ea_temporal_mean_TREFHT_base-ea_temporal_mean_TREFHT_ref)./(2*[T_se; T_se(1,:)]));
[lon_null,lat_null]=find(Z<1);
lon_null=(lon_null-1)/287*360;
lat_null=(lat_null-1)/191*180-90;

ax3=subplot(2,2,3); 
hold on
worldmap([-90 90],[-180 180])
setm(ax3, 'meridianlabel', 'off', 'parallellabel', 'off')
box on
contourfm(lat,lon2,ea_temporal_mean_TREFHT_base'-ea_temporal_mean_TREFHT_ref',T_contours,'LineStyle','none')
scatterm(lat_null,lon_null,.2,'k','filled');
load coastlines
geoshow(coastlat,coastlon,'Color','k')
caxis([Cmin Cmax])
cb=colorbar;
ylabel(cb,'\DeltaK')
hold off
title('\DeltaT from Reference Period 2050-2069 SSP2-4.5');
colormap(ax3, T_color_map);

Cmin = -3;
Cmax = 3;
% baseline PREFHT error check
Z=abs((ea_temporal_mean_PRECT_base-ea_temporal_mean_PRECT_ref)./(2*[P_se; P_se(1,:)]));
[lon_null,lat_null]=find(Z<1);
lon_null=(lon_null-1)/287*360;
lat_null=(lat_null-1)/191*180-90;

ax4=subplot(2,2,4);
hold on
worldmap([-90 90],[-180 180])
setm(ax4, 'meridianlabel', 'off', 'parallellabel', 'off')
box on
contourfm(lat,lon2,double(ea_temporal_mean_PRECT_base')-ea_temporal_mean_PRECT_ref',P_contours,'LineStyle','none')
scatterm(lat_null,lon_null,.2,'k','filled');
load coastlines
geoshow(coastlat,coastlon,'Color','k')
caxis([Cmin Cmax])
cb=colorbar;
ylabel(cb,'\Deltamm/day')
hold off
title('\DeltaP from Reference Period 2050-2069 SSP2-4.5');
colormap(ax4, P_color_map);

%% plot TREFHT anomaly maps
Cmin = -3;
Cmax = 3;
[contours, color_map] = custom_color_map(Cmin,Cmax,25,'*RdBu');

figure('Renderer', 'painters', 'Position', [20 60 1500 800])
clf

Z=abs((ea_temporal_mean_TREFHT_def-ea_temporal_mean_TREFHT_ref)./(2*[T_se; T_se(1,:)]));
[lon_null,lat_null]=find(Z<1);
lon_null=(lon_null-1)/287*360;
lat_null=(lat_null-1)/191*180-90;

ax1=subplot(3,3,1); %multi-objective
ax1.Position = [0.03 0.67 0.28 0.3];
hold on
worldmap([-90 90],[-180 180])
setm(ax1, 'meridianlabel', 'off', 'parallellabel', 'off')
box on
contourfm(lat,lon2,double(ea_temporal_mean_TREFHT_def')-ea_temporal_mean_TREFHT_ref',contours,'LineStyle','none')
scatterm(lat_null,lon_null,.2,'k','filled');
load coastlines
geoshow(coastlat,coastlon,'Color','k')
caxis([Cmin Cmax])
hold off
title('\DeltaT from Reference Period 2050-2069 Multi-Objective');
colormap(ax1, color_map);

Z=abs((ea_temporal_mean_TREFHT_EQ-ea_temporal_mean_TREFHT_ref)./(2*[T_se; T_se(1,:)]));
[lon_null,lat_null]=find(Z<1);
lon_null=(lon_null-1)/287*360;
lat_null=(lat_null-1)/191*180-90;

ax2=subplot(3,3,2); %EQ
ax2.Position = [0.36 0.67 0.28 0.3];
hold on
worldmap([-90 90],[-180 180])
setm(ax2, 'meridianlabel', 'off', 'parallellabel', 'off')
box on
contourfm(lat,lon2,double(ea_temporal_mean_TREFHT_EQ')-ea_temporal_mean_TREFHT_ref',contours,'LineStyle','none')
scatterm(lat_null,lon_null,.2,'k','filled');
load coastlines
geoshow(coastlat,coastlon,'Color','k')
caxis([Cmin Cmax])
hold off
title('\DeltaT from Reference Period 2050-2069 EQ');
colormap(ax2, color_map);

Z=abs((ea_temporal_mean_TREFHT_15N_15S-ea_temporal_mean_TREFHT_ref)./(2*[T_se; T_se(1,:)]));
[lon_null,lat_null]=find(Z<1);
lon_null=(lon_null-1)/287*360;
lat_null=(lat_null-1)/191*180-90;

ax3=subplot(3,3,3); %15N-15S
ax3.Position = [0.69 0.67 0.28 0.3];
hold on
worldmap([-90 90],[-180 180])
setm(ax3, 'meridianlabel', 'off', 'parallellabel', 'off')
box on
contourfm(lat,lon2,double(ea_temporal_mean_TREFHT_15N_15S')-ea_temporal_mean_TREFHT_ref',contours,'LineStyle','none')
scatterm(lat_null,lon_null,.2,'k','filled');
load coastlines
geoshow(coastlat,coastlon,'Color','k')
caxis([Cmin Cmax])
hold off
title('\DeltaT from Reference Period 2050-2069 15N-15S');
colormap(ax3, color_map);

Z=abs((ea_temporal_mean_TREFHT_30N_30S-ea_temporal_mean_TREFHT_ref)./(2*[T_se; T_se(1,:)]));
[lon_null,lat_null]=find(Z<1);
lon_null=(lon_null-1)/287*360;
lat_null=(lat_null-1)/191*180-90;

ax4=subplot(3,3,4); %30N-30S
ax4.Position = [0.03 0.33 0.28 0.3];
hold on
worldmap([-90 90],[-180 180])
setm(ax4, 'meridianlabel', 'off', 'parallellabel', 'off')
box on
contourfm(lat,lon2,double(ea_temporal_mean_TREFHT_30N_30S')-ea_temporal_mean_TREFHT_ref',contours,'LineStyle','none')
scatterm(lat_null,lon_null,.2,'k','filled');
load coastlines
geoshow(coastlat,coastlon,'Color','k')
caxis([Cmin Cmax])
hold off
title('\DeltaT from Reference Period 2050-2069 30N-30S');
colormap(ax4, color_map);

Z=abs((ea_temporal_mean_TREFHT_60N_60S-ea_temporal_mean_TREFHT_ref)./(2*[T_se; T_se(1,:)]));
[lon_null,lat_null]=find(Z<1);
lon_null=(lon_null-1)/287*360;
lat_null=(lat_null-1)/191*180-90;

ax5=subplot(3,3,5); %60N-60S
ax5.Position = [0.36 0.33 0.28 0.3];
hold on
worldmap([-90 90],[-180 180])
setm(ax5, 'meridianlabel', 'off', 'parallellabel', 'off')
box on
contourfm(lat,lon2,double(ea_temporal_mean_TREFHT_60N_60S')-ea_temporal_mean_TREFHT_ref',contours,'LineStyle','none')
scatterm(lat_null,lon_null,.2,'k','filled');
load coastlines
geoshow(coastlat,coastlon,'Color','k')
caxis([Cmin Cmax])
hold off
title('\DeltaT from Reference Period 2050-2069 60N-60S');
colormap(ax5, color_map);

Z=abs((ea_temporal_mean_TREFHT_30N-ea_temporal_mean_TREFHT_ref)./(2*[T_se; T_se(1,:)]));
[lon_null,lat_null]=find(Z<1);
lon_null=(lon_null-1)/287*360;
lat_null=(lat_null-1)/191*180-90;

ax6=subplot(3,3,6); %30N
ax6.Position = [0.69 0.33 0.28 0.3];
hold on
worldmap([-90 90],[-180 180])
setm(ax6, 'meridianlabel', 'off', 'parallellabel', 'off')
box on
contourfm(lat,lon2,double(ea_temporal_mean_TREFHT_30N')-ea_temporal_mean_TREFHT_ref',contours,'LineStyle','none')
scatterm(lat_null,lon_null,.2,'k','filled');
load coastlines
geoshow(coastlat,coastlon,'Color','k')
caxis([Cmin Cmax])
hold off
title('\DeltaT from Reference Period 2050-2069 30N');
colormap(ax6, color_map);

Z=abs((ea_temporal_mean_TREFHT_60N-ea_temporal_mean_TREFHT_ref)./(2*[T_se; T_se(1,:)]));
[lon_null,lat_null]=find(Z<1);
lon_null=(lon_null-1)/287*360;
lat_null=(lat_null-1)/191*180-90;

ax7=subplot(3,3,7); %60N
ax7.Position = [0.03 0 0.32 0.3];
hold on
worldmap([-90 90],[-180 180])
setm(ax7, 'meridianlabel', 'off', 'parallellabel', 'off')
box on
contourfm(lat,lon2,double(ea_temporal_mean_TREFHT_60N')-ea_temporal_mean_TREFHT_ref',contours,'LineStyle','none')
scatterm(lat_null,lon_null,.2,'k','filled');
load coastlines
geoshow(coastlat,coastlon,'Color','k')
caxis([Cmin Cmax])
cb=colorbar;
ylabel(cb,'\DeltaK')
hold off
title('\DeltaT from Reference Period 2050-2069 60N');
colormap(ax7, color_map);

%% plot PRECT anomaly maps
Cmin = -3;
Cmax = 3;
[contours, color_map] = custom_color_map(Cmin,Cmax,11,'BrBG');

figure('Renderer', 'painters', 'Position', [20 60 1500 800])
clf

Z=abs((ea_temporal_mean_PRECT_def-ea_temporal_mean_PRECT_ref)./(2*[P_se; P_se(1,:)]));
[lon_null,lat_null]=find(Z<1);
lon_null=(lon_null-1)/287*360;
lat_null=(lat_null-1)/191*180-90;

ax1=subplot(3,3,1); %multi-objective
ax1.Position = [0.03 0.67 0.28 0.3];
hold on
worldmap([-90 90],[-180 180])
setm(ax1, 'meridianlabel', 'off', 'parallellabel', 'off')
box on
contourfm(lat,lon2,double(ea_temporal_mean_PRECT_def')-ea_temporal_mean_PRECT_ref',contours,'LineStyle','none')
scatterm(lat_null,lon_null,.2,'k','filled');
load coastlines
geoshow(coastlat,coastlon,'Color','k')
caxis([Cmin Cmax])
hold off
title('\Deltaprecip from Reference Period 2050-2069 Multi-Objective');
colormap(ax1, color_map);

Z=abs((ea_temporal_mean_PRECT_EQ-ea_temporal_mean_PRECT_ref)./(2*[P_se; P_se(1,:)]));
[lon_null,lat_null]=find(Z<1);
lon_null=(lon_null-1)/287*360;
lat_null=(lat_null-1)/191*180-90;

ax2=subplot(3,3,2); %EQ
ax2.Position = [0.36 0.67 0.28 0.3];
hold on
worldmap([-90 90],[-180 180])
setm(ax2, 'meridianlabel', 'off', 'parallellabel', 'off')
box on
contourfm(lat,lon2,double(ea_temporal_mean_PRECT_EQ')-ea_temporal_mean_PRECT_ref',contours,'LineStyle','none')
scatterm(lat_null,lon_null,.2,'k','filled');
load coastlines
geoshow(coastlat,coastlon,'Color','k')
caxis([Cmin Cmax])
hold off
title('\Deltaprecip from Reference Period 2050-2069 EQ');
colormap(ax2, color_map);

Z=abs((ea_temporal_mean_PRECT_15N_15S-ea_temporal_mean_PRECT_ref)./(2*[P_se; P_se(1,:)]));
[lon_null,lat_null]=find(Z<1);
lon_null=(lon_null-1)/287*360;
lat_null=(lat_null-1)/191*180-90;

ax3=subplot(3,3,3); %15N-15S
ax3.Position = [0.69 0.67 0.28 0.3];
hold on
worldmap([-90 90],[-180 180])
setm(ax3, 'meridianlabel', 'off', 'parallellabel', 'off')
box on
contourfm(lat,lon2,double(ea_temporal_mean_PRECT_15N_15S')-ea_temporal_mean_PRECT_ref',contours,'LineStyle','none')
scatterm(lat_null,lon_null,.2,'k','filled');
load coastlines
geoshow(coastlat,coastlon,'Color','k')
caxis([Cmin Cmax])
hold off
title('\Deltaprecip from Reference Period 2050-2069 15N-15S');
colormap(ax3, color_map);

Z=abs((ea_temporal_mean_PRECT_30N_30S-ea_temporal_mean_PRECT_ref)./(2*[P_se; P_se(1,:)]));
[lon_null,lat_null]=find(Z<1);
lon_null=(lon_null-1)/287*360;
lat_null=(lat_null-1)/191*180-90;

ax4=subplot(3,3,4); %30N-30S
ax4.Position = [0.03 0.33 0.28 0.3];
hold on
worldmap([-90 90],[-180 180])
setm(ax4, 'meridianlabel', 'off', 'parallellabel', 'off')
box on
contourfm(lat,lon2,double(ea_temporal_mean_PRECT_30N_30S')-ea_temporal_mean_PRECT_ref',contours,'LineStyle','none')
scatterm(lat_null,lon_null,.2,'k','filled');
load coastlines
geoshow(coastlat,coastlon,'Color','k')
caxis([Cmin Cmax])
hold off
title('\Deltaprecip from Reference Period 2050-2069 30N-30S');
colormap(ax4, color_map);

Z=abs((ea_temporal_mean_PRECT_60N_60S-ea_temporal_mean_PRECT_ref)./(2*[P_se; P_se(1,:)]));
[lon_null,lat_null]=find(Z<1);
lon_null=(lon_null-1)/287*360;
lat_null=(lat_null-1)/191*180-90;

ax5=subplot(3,3,5); %60N-60S
ax5.Position = [0.36 0.33 0.28 0.3];
hold on
worldmap([-90 90],[-180 180])
setm(ax5, 'meridianlabel', 'off', 'parallellabel', 'off')
box on
contourfm(lat,lon2,double(ea_temporal_mean_PRECT_60N_60S')-ea_temporal_mean_PRECT_ref',contours,'LineStyle','none')
scatterm(lat_null,lon_null,.2,'k','filled');
load coastlines
geoshow(coastlat,coastlon,'Color','k')
caxis([Cmin Cmax])
hold off
title('\Deltaprecip from Reference Period 2050-2069 60N-60S');
colormap(ax5, color_map);

Z=abs((ea_temporal_mean_PRECT_30N-ea_temporal_mean_PRECT_ref)./(2*[P_se; P_se(1,:)]));
[lon_null,lat_null]=find(Z<1);
lon_null=(lon_null-1)/287*360;
lat_null=(lat_null-1)/191*180-90;

ax6=subplot(3,3,6); %30N
ax6.Position = [0.69 0.33 0.28 0.3];
hold on
worldmap([-90 90],[-180 180])
setm(ax6, 'meridianlabel', 'off', 'parallellabel', 'off')
box on
contourfm(lat,lon2,double(ea_temporal_mean_PRECT_30N')-ea_temporal_mean_PRECT_ref',contours,'LineStyle','none')
scatterm(lat_null,lon_null,.2,'k','filled');
load coastlines
geoshow(coastlat,coastlon,'Color','k')
caxis([Cmin Cmax])
hold off
title('\Deltaprecip from Reference Period 2050-2069 30N');
colormap(ax6, color_map);

Z=abs((ea_temporal_mean_PRECT_60N-ea_temporal_mean_PRECT_ref)./(2*[P_se; P_se(1,:)]));
[lon_null,lat_null]=find(Z<1);
lon_null=(lon_null-1)/287*360;
lat_null=(lat_null-1)/191*180-90;

ax7=subplot(3,3,7); %60N
ax7.Position = [0.03 0 0.32 0.3];
hold on
worldmap([-90 90],[-180 180])
setm(ax7, 'meridianlabel', 'off', 'parallellabel', 'off')
box on
contourfm(lat,lon2,double(ea_temporal_mean_PRECT_60N')-ea_temporal_mean_PRECT_ref',contours,'LineStyle','none')
scatterm(lat_null,lon_null,.2,'k','filled');
load coastlines
geoshow(coastlat,coastlon,'Color','k')
caxis([Cmin Cmax])
cb=colorbar;
ylabel(cb,'\Deltamm/day')
hold off
title('\Deltaprecip from Reference Period 2050-2069 60N');
colormap(ax7, color_map);

%% calculate RMS error for the different strategies
area_weight=AREA/sum(AREA,'all'); % 288 x 192

% for the 2-member-pair ensemble averages, the order is 1+2, 2+3, 1+3

% T
for k = 2:7
    eval(['rms_T_' scenarios{k} ' = squeeze(sqrt(sum(area_weight.*((mu_TREFHT_' scenarios{k} '(1:288,:,:)+alpha_TREFHT(1:288,:,:)).^2),[1 2])));'])
    eval(['ea_rms_T_' scenarios{k} ' = squeeze(sqrt(sum(area_weight.*((mean(mu_TREFHT_' scenarios{k} '(1:288,:,:),3)+mean(alpha_TREFHT(1:288,:,:),3)).^2),[1 2])));'])
    eval(['rms_T_' scenarios{k} '_2(1) = squeeze(sqrt(sum(area_weight.*((mean(mu_TREFHT_' scenarios{k} '(1:288,:,[1 2]),3)+mean(alpha_TREFHT(1:288,:,[1 2]),3)).^2),[1 2])));'])
    eval(['rms_T_' scenarios{k} '_2(2) = squeeze(sqrt(sum(area_weight.*((mean(mu_TREFHT_' scenarios{k} '(1:288,:,[2 3]),3)+mean(alpha_TREFHT(1:288,:,[2 3]),3)).^2),[1 2])));'])
    eval(['rms_T_' scenarios{k} '_2(3) = squeeze(sqrt(sum(area_weight.*((mean(mu_TREFHT_' scenarios{k} '(1:288,:,[1 3]),3)+mean(alpha_TREFHT(1:288,:,[1 3]),3)).^2),[1 2])));'])
end
rms_T_base = squeeze(sqrt(sum(area_weight.*((alpha_TREFHT(1:288,:,:)).^2),[1 2])));
ea_rms_T_base = squeeze(sqrt(sum(area_weight.*((mean(alpha_TREFHT(1:288,:,:),3)).^2),[1 2])));
ea_rms_T_30N = squeeze(sqrt(sum(area_weight.*((mean(mu_TREFHT_30N(1:288,:,:),3)+mean(alpha_TREFHT(1:288,:,:),3)).^2),[1 2])));
ea_rms_T_60N = squeeze(sqrt(sum(area_weight.*((mean(mu_TREFHT_60N(1:288,:,:),3)+mean(alpha_TREFHT(1:288,:,:),3)).^2),[1 2])));
% P
for k = 2:7
    eval(['rms_P_' scenarios{k} ' = squeeze(sqrt(sum(area_weight.*((mu_PRECT_' scenarios{k} '(1:288,:,:)+alpha_PRECT(1:288,:,:)).^2),[1 2])));'])
    eval(['ea_rms_P_' scenarios{k} ' = squeeze(sqrt(sum(area_weight.*((mean(mu_PRECT_' scenarios{k} '(1:288,:,:),3)+mean(alpha_PRECT(1:288,:,:),3)).^2),[1 2])));'])
    eval(['rms_P_' scenarios{k} '_2(1) = squeeze(sqrt(sum(area_weight.*((mean(mu_PRECT_' scenarios{k} '(1:288,:,[1 2]),3)+mean(alpha_PRECT(1:288,:,[1 2]),3)).^2),[1 2])));'])
    eval(['rms_P_' scenarios{k} '_2(2) = squeeze(sqrt(sum(area_weight.*((mean(mu_PRECT_' scenarios{k} '(1:288,:,[2 3]),3)+mean(alpha_PRECT(1:288,:,[2 3]),3)).^2),[1 2])));'])
    eval(['rms_P_' scenarios{k} '_2(3) = squeeze(sqrt(sum(area_weight.*((mean(mu_PRECT_' scenarios{k} '(1:288,:,[1 3]),3)+mean(alpha_PRECT(1:288,:,[1 3]),3)).^2),[1 2])));'])
end
rms_P_base = squeeze(sqrt(sum(area_weight.*((alpha_PRECT(1:288,:,:)).^2),[1 2])));
ea_rms_P_base = squeeze(sqrt(sum(area_weight.*((mean(alpha_PRECT(1:288,:,:),3)).^2),[1 2])));
ea_rms_P_30N = squeeze(sqrt(sum(area_weight.*((mean(mu_PRECT_30N(1:288,:,:),3)+mean(alpha_PRECT(1:288,:,:),3)).^2),[1 2])));
ea_rms_P_60N = squeeze(sqrt(sum(area_weight.*((mean(mu_PRECT_60N(1:288,:,:),3)+mean(alpha_PRECT(1:288,:,:),3)).^2),[1 2])));
addpath('C:\Users\ezrab\Documents\Cornell\Research\Code\SAI_optimization_paper-main\SAI_optimization_paper-main')
fname_landmask='sftlf_CESM-CAM5.1-FV.nc';
ncid=netcdf.open(fname_landmask,'NC_NOWRITE');
landmask=double(netcdf.getVar(ncid,netcdf.inqVarID(ncid,'sftlf')))/100;
landmask_vec=reshape(landmask,288*192,[]);
land_area=AREA.*landmask;
land_weight=land_area/sum(land_area,'all');

% land T
for k = 2:7
    eval(['rms_land_T_' scenarios{k} ' = squeeze(sqrt(sum(land_weight.*((mu_TREFHT_' scenarios{k} '(1:288,:,:)+alpha_TREFHT(1:288,:,:)).^2),[1 2])));'])
    eval(['ea_rms_land_T_' scenarios{k} ' = squeeze(sqrt(sum(land_weight.*((mean(mu_TREFHT_' scenarios{k} '(1:288,:,:),3)+mean(alpha_TREFHT(1:288,:,:),3)).^2),[1 2])));'])
end

% land P-E
for k = 2:7
    eval(['rms_land_PE_' scenarios{k} ' = squeeze(sqrt(sum(land_weight.*((mu_PRECT_' scenarios{k} '(1:288,:,:)-mu_QFLX_' scenarios{k} '(1:288,:,:)+alpha_PRECT(1:288,:,:)-alpha_QFLX(1:288,:,:)).^2),[1 2])));'])
    eval(['ea_rms_land_PE_' scenarios{k} ' = squeeze(sqrt(sum(land_weight.*((mean(mu_PRECT_' scenarios{k} '(1:288,:,:)-mu_QFLX_' scenarios{k} '(1:288,:,:),3)+mean(alpha_PRECT(1:288,:,:)-alpha_QFLX(1:288,:,:),3)).^2),[1 2])));'])
end

load('mu_predicted.mat')
for k = 1:length(varnames)
    eval(['mu_' varnames{k} '_predicted = mu_' varnames{k} '_predicted(1:288,:);'])
end

rms_T_predicted = squeeze(sqrt(sum(area_weight.*((mu_TREFHT_predicted+mean(alpha_TREFHT(1:288,:,:),3)).^2),[1 2])));
rms_P_predicted = squeeze(sqrt(sum(area_weight.*((mu_PRECT_predicted+mean(alpha_PRECT(1:288,:,:),3)).^2),[1 2])));
rms_land_T_predicted = squeeze(sqrt(sum(land_weight.*((mu_TREFHT_predicted+mean(alpha_TREFHT(1:288,:,:),3)).^2),[1 2])));
rms_land_PE_predicted = squeeze(sqrt(sum(land_weight.*((mu_PRECT_predicted-mu_QFLX_predicted+mean(alpha_PRECT(1:288,:,:)-alpha_QFLX(1:288,:,:),3)).^2),[1 2])));

% below here is normalized per standard error
% T std error
for k = 2:7
    eval(['rms_T_' scenarios{k} '_norm = squeeze(sqrt(sum(area_weight.*(((mu_TREFHT_' scenarios{k} '(1:288,:,:)+alpha_TREFHT(1:288,:,:))./T_std).^2),[1 2])));'])
    eval(['ea_rms_T_' scenarios{k} '_norm = squeeze(sqrt(sum(area_weight.*(((mean(mu_TREFHT_' scenarios{k} '(1:288,:,:),3)+mean(alpha_TREFHT(1:288,:,:),3))./T_std).^2),[1 2])));'])
end
rms_T_base_norm = squeeze(sqrt(sum(area_weight.*(((alpha_TREFHT(1:288,:,:))./T_std).^2),[1 2])));
ea_rms_T_base_norm = squeeze(sqrt(sum(area_weight.*(((mean(alpha_TREFHT(1:288,:,:),3))./T_std).^2),[1 2])));
ea_rms_T_30N_norm = squeeze(sqrt(sum(area_weight.*(((mean(mu_TREFHT_30N(1:288,:,:),3)+mean(alpha_TREFHT(1:288,:,:),3))./T_std).^2),[1 2])));
ea_rms_T_60N_norm = squeeze(sqrt(sum(area_weight.*(((mean(mu_TREFHT_60N(1:288,:,:),3)+mean(alpha_TREFHT(1:288,:,:),3))./T_std).^2),[1 2])));

% P std error
for k = 2:7
    eval(['rms_P_' scenarios{k} '_norm = squeeze(sqrt(sum(area_weight.*(((mu_PRECT_' scenarios{k} '(1:288,:,:)+alpha_PRECT(1:288,:,:))./P_std).^2),[1 2])));'])
    eval(['ea_rms_P_' scenarios{k} '_norm = squeeze(sqrt(sum(area_weight.*(((mean(mu_PRECT_' scenarios{k} '(1:288,:,:),3)+mean(alpha_PRECT(1:288,:,:),3))./P_std).^2),[1 2])));'])
end
rms_P_base_norm = squeeze(sqrt(sum(area_weight.*(((alpha_PRECT(1:288,:,:))./P_std).^2),[1 2])));
ea_rms_P_base_norm = squeeze(sqrt(sum(area_weight.*(((mean(alpha_PRECT(1:288,:,:),3))./P_std).^2),[1 2])));
ea_rms_P_30N_norm = squeeze(sqrt(sum(area_weight.*(((mean(mu_PRECT_30N(1:288,:,:),3)+mean(alpha_PRECT(1:288,:,:),3))./P_std).^2),[1 2])));
ea_rms_P_60N_norm = squeeze(sqrt(sum(area_weight.*(((mean(mu_PRECT_60N(1:288,:,:),3)+mean(alpha_PRECT(1:288,:,:),3))./P_std).^2),[1 2])));

% land T std error
for k = 2:7
    eval(['rms_land_T_' scenarios{k} '_norm = squeeze(sqrt(sum(land_weight.*(((mu_TREFHT_' scenarios{k} '(1:288,:,:)+alpha_TREFHT(1:288,:,:))./T_std).^2),[1 2])));'])
    eval(['ea_rms_land_T_' scenarios{k} '_norm = squeeze(sqrt(sum(land_weight.*(((mean(mu_TREFHT_' scenarios{k} '(1:288,:,:),3)+mean(alpha_TREFHT(1:288,:,:),3))./T_std).^2),[1 2])));'])
end
rms_land_T_base_norm = squeeze(sqrt(sum(land_weight.*(((alpha_TREFHT(1:288,:,:))./T_std).^2),[1 2])));
ea_rms_land_T_base_norm = squeeze(sqrt(sum(land_weight.*(((mean(alpha_TREFHT(1:288,:,:),3))./T_std).^2),[1 2])));
ea_rms_land_T_30N_norm = squeeze(sqrt(sum(land_weight.*(((mean(mu_TREFHT_30N(1:288,:,:),3)+mean(alpha_TREFHT(1:288,:,:),3))./T_std).^2),[1 2])));
ea_rms_land_T_60N_norm = squeeze(sqrt(sum(land_weight.*(((mean(mu_TREFHT_60N(1:288,:,:),3)+mean(alpha_TREFHT(1:288,:,:),3))./T_std).^2),[1 2])));

% land P-E std error
for k = 2:7
    eval(['rms_land_PE_' scenarios{k} '_norm = squeeze(sqrt(sum(land_weight.*(((mu_PRECT_' scenarios{k} '(1:288,:,:)-mu_QFLX_' scenarios{k} '(1:288,:,:)+alpha_PRECT(1:288,:,:)-alpha_QFLX(1:288,:,:))./PE_std).^2),[1 2])));'])
    eval(['ea_rms_land_PE_' scenarios{k} '_norm = squeeze(sqrt(sum(land_weight.*(((mean(mu_PRECT_' scenarios{k} '(1:288,:,:)-mu_QFLX_' scenarios{k} '(1:288,:,:),3)+mean(alpha_PRECT(1:288,:,:)-alpha_QFLX(1:288,:,:),3))./PE_std).^2),[1 2])));'])
end

rms_T_norm_predicted = squeeze(sqrt(sum(area_weight.*(((mu_TREFHT_predicted+mean(alpha_TREFHT(1:288,:,:),3))./T_std).^2),[1 2])));
rms_P_norm_predicted = squeeze(sqrt(sum(area_weight.*(((mu_PRECT_predicted+mean(alpha_PRECT(1:288,:,:),3))./P_std).^2),[1 2])));
rms_land_T_norm_predicted = squeeze(sqrt(sum(land_weight.*(((mu_TREFHT_predicted+mean(alpha_TREFHT(1:288,:,:),3))./T_std).^2),[1 2])));
rms_land_PE_norm_predicted = squeeze(sqrt(sum(land_weight.*(((mu_PRECT_predicted-mu_QFLX_predicted+mean(alpha_PRECT(1:288,:,:)-alpha_QFLX(1:288,:,:),3))./PE_std).^2),[1 2])));

%% Setup optimization (3 ens)
load('Delta_T0_cesm2.mat');
% order of strategies: EQ, 15N_15S, 30N_30S, 60N_60S, 60N, multi-objective, 30N
for k = 1:length(varnames)
    eval(['alpha_' varnames{k} '_vec = reshape(mean(alpha_' varnames{k} '(1:288,:,:),3),[55296 1]);'])
    eval(['mu_' varnames{k} '_arr = [reshape(mean(mu_' varnames{k} '_EQ(1:288,:,:),3),[55296 1])' ...
        ' reshape(mean(mu_' varnames{k} '_15N_15S(1:288,:,:),3),[55296 1])' ...
        ' reshape(mean(mu_' varnames{k} '_30N_30S(1:288,:,:),3),[55296 1])' ...
        ' reshape(mean(mu_' varnames{k} '_60N_60S(1:288,:,:),3),[55296 1])' ...
        ' reshape(mean(mu_' varnames{k} '_60N(1:288,:,:),3),[55296 1])' ...
        ' reshape(mean(mu_' varnames{k} '_def(1:288,:,:),3),[55296 1])' ...
        ' reshape(mean(mu_' varnames{k} '_30N(1:288,:,:),3),[55296 1])];'])
end

Delta_arr = [ea_inj_EQ_end, ea_inj_15N_15S_end, ea_inj_30N_30S_end, ea_inj_60N_60S_end, 12, sum(ea_inj_def_end), 12]...
        ./[mean(cooling_EQ) mean(cooling_15N_15S) mean(cooling_30N_30S) mean(cooling_60N_60S) mean(cooling_60N) mean(cooling_def) mean(cooling_30N)];
delta_Global_one_frac = ea_inj_def_end./sum(ea_inj_def_end);
save('delta.mat',"delta_Global_one_frac", "Delta_arr")
area_weight_vec=double(reshape(area_weight,288*192,[]));
area_weight_vec=[area_weight_vec;area_weight_vec];

fname_landmask='sftlf_CESM-CAM5.1-FV.nc';
ncid=netcdf.open(fname_landmask,'NC_NOWRITE');
landmask=double(netcdf.getVar(ncid,netcdf.inqVarID(ncid,'sftlf')))/100;
landmask_vec=reshape(landmask,288*192,[]);

land_area=double(AREA.*landmask);
land_weight=land_area/sum(land_area,'all');
land_weight_vec=reshape(land_weight,288*192,[]);
land_weight_vec=[land_weight_vec;land_weight_vec];

C=[0 0 0 0.5*Delta_arr(4) Delta_arr(5) 0 0 ;
       0 0 0.5*Delta_arr(3) 0 0 delta_Global_one_frac(4)*Delta_arr(6) Delta_arr(7) ;
       0 0.5*Delta_arr(2) 0 0 0 delta_Global_one_frac(3)*Delta_arr(6) 0;
       Delta_arr(1) 0 0 0 0 0 0;
       0 0.5*Delta_arr(2) 0 0 0 delta_Global_one_frac(2)*Delta_arr(6) 0 ;
       0 0 0.5*Delta_arr(3) 0 0 delta_Global_one_frac(1)*Delta_arr(6) 0;
       0 0 0 0.5*Delta_arr(4) 0 0 0];

%% perform optimization T, T1, T2, ITCZ, SSI
inj_each_lat_arr_5DOF=zeros(7,1);
% residuals
T1_residual_arr=zeros(1,3);
T2_residual_arr=zeros(1,3);
ITCZ_residual_arr=zeros(1,3);
SSI_residual_arr=zeros(1,3);

% beta values
beta_T1=mu_T1_arr+ea_alpha_T1;%1x7
beta_T2=mu_T2_arr+ea_alpha_T2;%1x7
beta_ITCZ=mu_ITCZ_arr+ea_alpha_ITCZ; %1X7
beta_SSI=mu_SSI_arr+alpha_SSI;%1x7

% beta matrix
%beta_matrix=[beta_T1;beta_T2;beta_ITCZ;beta_SSI];

% beta_matrix_w=[beta_P0/std(beta_P0);
%     beta_T1/std(beta_T1);
%     beta_T2/std(beta_T2);
%     beta_AMOC/std(beta_AMOC)];
beta_matrix_w=[beta_T1/T1_std_err;
    beta_T2/T2_std_err;
    beta_ITCZ/ITCZ_std_err;
    beta_SSI/SSI_std_err];

for i=1
    Delta_T0=i;
    % force injection at each latitude to be nonnegative
%     % force injection at equator to be 0
%     C1=[-Delta_arr(1) 0 0 0 0 0 0];
    % f=0.05*Delta_arr';
    f=zeros(7,1);
    % epsilon=0.001*ones(7,1);
    
    % inj_each_lat_arr=zeros(7,9);
    H=double((beta_matrix_w)'*(beta_matrix_w));
    
%     fun=@(x) sqrt((beta_matrix*x)'*(beta_matrix*x));

%     % force injection at equator to be 0
%     A=-[C;C1];
%     b=zeros(8,1);

    A=-C;
    b=zeros(7,1);
    % enforce global mean temperature rise is offset
    Aeq=ones(1,7);
    beq=Delta_T0;
    % b=epsilon;
%     x0=[0.1 0.1 0.3 0 0.3 0.1 0.1]';
    lb=[];
    ub=[];
%     [x,fval]=fmincon(fun,x0,A,b,Aeq,beq,lb,ub)
    options = optimoptions(@quadprog,'Algorithm','active-set');
    [x,fval,exitflag,output]=quadprog(H,f,A,b,Aeq,beq,[],[],[0;0;0;0;0;1;0],options);
    % x is the fraction of injection for each SAI strategy
    % for 1 degC cooling, save the solution
    if i==1
        save('SAI_strategy_choices_T1_T2_ITCZ_SSI.mat','x');
    end
    % injection rate for 60N,30N,15N,0N,15S,30S,60S
    inj_each_lat=C*x;
    inj_each_lat_arr_5DOF(:,i)=inj_each_lat;
    % climate impact
%     climate_impact=sqrt(x'*beta_matrix_w'*beta_matrix_w*x);
    T1_residual_arr(i)=beta_T1*x;
    T2_residual_arr(i)=beta_T2*x;
    ITCZ_residual_arr(i)=beta_ITCZ*x;
    SSI_residual_arr(i)=beta_SSI*x;
    
end

alpha_vec = [alpha_TREFHT_vec; alpha_PRECT_vec];
beta_matrix=[mu_TREFHT_arr;mu_PRECT_arr]+alpha_vec;
beta_matrix_opt=beta_matrix.*sqrt(area_weight_vec);
temp_residual_rms_5DOF=sqrt(x'*beta_matrix_opt(1:55296,:)'*beta_matrix_opt(1:55296,:)*x);
prec_residual_rms_5DOF=sqrt(x'*beta_matrix_opt(55297:110592,:)'*beta_matrix_opt(55297:110592,:)*x);

% % save injection rates
% writetable(table(inj_each_lat_arr),'injection_rates_different_cooling_T1_T2_ITCZ_SSI.csv');
% % save residuals
% writetable(table([T1_residual_arr;T2_residual_arr;ITCZ_residual_arr;SSI_residual_arr]),'residual_different_cooling_T1_T2_ITCZ_SSI.csv')
T1_T2_ITCZ_SSI_opt=[T1_residual_arr(1),T2_residual_arr(1),ITCZ_residual_arr(1),SSI_residual_arr(1)];
save('residual_T1_T2_ITCZ_SSI.mat','T1_T2_ITCZ_SSI_opt');

% strategies
strategies={'Optimized','Multi-Objective','EQ','15N+15S','30N+30S','60N+60S','60N','30N'};
X = categorical(strategies);
X = reordercats(X,strategies);

T1_data=[T1_T2_ITCZ_SSI_opt(1),beta_T1(6), beta_T1(1:5), beta_T1(7)];
T2_data=[T1_T2_ITCZ_SSI_opt(2),beta_T2(6), beta_T2(1:5),beta_T2(7)];
ITCZ_data=[T1_T2_ITCZ_SSI_opt(3),beta_ITCZ(6), beta_ITCZ(1:5),beta_ITCZ(7)];
SSI_data=[T1_T2_ITCZ_SSI_opt(4),beta_SSI(6), beta_SSI(1:5),beta_SSI(7)];

% colors: grey, orange, yellow, green, blue, cyan, purple, lightgreen
colors=[0.5,0.5,0.5;
    0.7176 0.2745 1.0000;
    0.8500 0.3250 0.0980;
    0.9290, 0.6940, 0.1250;
    0.4660 0.6740 0.1880;
    0 0.4470 0.7410;
    0 1 1;
    0 1 0];

figure(1)
clf
set(gcf,'PaperPosition',[0.5 1 13 6]);
set(gcf,'Units','inches');
set(gcf,'Position',get(gcf,'PaperPosition'))
set(gca,'FontSize',11)
idx=1;

subplot(2,3,1)
hold on
b1=bar(X,T1_data,0.6,'FaceColor','flat');
b1.CData=colors(1:8,:);
% plot the optimal solution with a star shape
scatter(X(1),T1_data(1),40,'k',"pentagram")
plot(xlim,2*[T1_std_err T1_std_err],'k--','LineWidth',0.5)
plot(xlim,2*[-T1_std_err -T1_std_err],'k--','LineWidth',0.5)
ylabel('\Delta T1 [K]')
ylim([-0.6 0.1])
grid on
hold off 
add_abc_tofig(idx,-0.01,-0.01,12);
idx=idx+1;

subplot(2,3,2)
hold on
b2=bar(X,T2_data,0.6,'FaceColor','flat');
b2.CData=colors(1:8,:);
% plot the optimal solution with a star shape
scatter(X(1),T2_data(1),40,'k',"pentagram")
plot(xlim,2*[T2_std_err T2_std_err],'k--','LineWidth',0.5)
plot(xlim,2*[-T2_std_err -T2_std_err],'k--','LineWidth',0.5)
ylabel('\Delta T2 [K]')
grid on
hold off 
add_abc_tofig(idx,-0.01,-0.01,12);
idx=idx+1;

subplot(2,3,3)
hold on
b3=bar(X,ITCZ_data,0.6,'FaceColor','flat');
b3.CData=colors(1:8,:);
% plot the optimal solution with a star shape
scatter(X(1),ITCZ_data(1),40,'k',"pentagram")
plot(xlim,2*[ITCZ_std_err ITCZ_std_err],'k--','LineWidth',0.5)
plot(xlim,2*[-ITCZ_std_err -ITCZ_std_err],'k--','LineWidth',0.5)
ylabel('\Delta ITCZ [degree]')
grid on
hold off 
add_abc_tofig(idx,-0.01,-0.01,12);
idx=idx+1;

subplot(2,3,4)
hold on
b4=bar(X,SSI_data/1000000,0.6,'FaceColor','flat');
b4.CData=colors(1:8,:);
% plot the optimal solution with a star shape
scatter(X(1),SSI_data(1)/1000000,40,'k',"pentagram")
plot(xlim,2*[SSI_std_err SSI_std_err]/1000000,'k--','LineWidth',0.5)
plot(xlim,2*[-SSI_std_err -SSI_std_err]/1000000,'k--','LineWidth',0.5)
ylabel('\Delta SSI [km^2]')
grid on
hold off 
add_abc_tofig(idx,-0.01,-0.01,12);
idx=idx+1;

inj_matrix = [flip(inj_each_lat_arr_5DOF)';
              0 ea_inj_def_end(1) ea_inj_def_end(2) 0 ea_inj_def_end(3) ea_inj_def_end(4) 0];
lats = categorical({'60S'; '30S'; '15S'; 'EQ'; '15N'; '30N'; '60N'});
lats = reordercats(lats,{'60S'; '30S'; '15S'; 'EQ'; '15N'; '30N'; '60N'});
ax=subplot(2,3,5);
ax.Position=[0.43 0.16 0.5 0.29];
hold on
b5=bar(lats,inj_matrix,'FaceColor','flat');
b5(1).CData = colors(1,:);
b5(2).CData = colors(2,:);
%xticklabels([-60 -30 -15 0 15 30 60])
ylabel('injection rate [Tg/yr]')
grid on
hold off 
add_abc_tofig(idx,-0.01,-0.01,12);
idx=idx+1;
legend('Optimized','Multi-Objective','Location','eastoutside')

saveas(gcf, 'C:\Users\ezrab\Documents\Cornell\Research\Figures\SAI_optimization\Fig_T1_T2_ITCZ_SSI_opt_vs_other_strategies.png')


%% Perform optimization rms TREFHT vs rms PRECT(3 ens)

alpha_vec = [alpha_TREFHT_vec; alpha_PRECT_vec];

Delta_T0=1;
inj_each_lat_arrTvsP_=zeros(7,11);
TvsP_x = zeros(7,11);
TvsP_temp_residual_rms_arr=zeros(1,11);
TvsP_prec_residual_rms_arr=zeros(1,11);
% temp_residual_rms_arr_land=zeros(1,11);
% prec_residual_rms_arr_land=zeros(1,11);

for i=0:0.1:1
    % weighting for temperature and precipitation
    W_t=i*ones(55296,1);
    W_p=(1-i)*ones(55296,1);

    beta_matrix=[mu_TREFHT_arr;mu_PRECT_arr]+alpha_vec;
    beta_matrix_TvsP=beta_matrix.*sqrt(area_weight_vec.*[W_t;W_p]);
    % f=0.05*Delta_arr';
    f=zeros(7,1);
    % epsilon=0.001*ones(7,1);
    
    H=(beta_matrix_TvsP)'*(beta_matrix_TvsP);
    
    % fun=@(x) sqrt((beta_matrix*x)'*(beta_matrix*x));
    
    A=-C;
    b=zeros(7,1);
    Aeq=ones(1,7);
    beq=Delta_T0;
    % b=epsilon;
    % x0=[0.1 0.1 0.3 0 0.3 0.1 0.1];
    lb=[];
    ub=[];
    % [x,fval]=fmincon(fun,x0,A,b,lb,ub)
    [x,fval,exitflag,output]=quadprog(H,f,A,b,Aeq,beq,[],[]);
    % save fraction of injection for each SAI strategy
    if i==0
        opt_frac_TvsP_prec = x;
    elseif i==1
        opt_frac_TvsP_temp = x;
    end
    TvsP_x(:,floor(10*i+1)) = x;
    % injection rate for 60N,30N,15N,0N,15S,30S,60S
    inj_each_lat=C*x;
    TvsP_inj_each_lat_arr(:,floor(10*i+1))=inj_each_lat;
    % climate impact
    beta_matrix_opt=beta_matrix.*sqrt(area_weight_vec);
    temp_residual_rms=sqrt(x'*beta_matrix_opt(1:55296,:)'*beta_matrix_opt(1:55296,:)*x);
    TvsP_temp_residual_rms_arr(floor(10*i+1))=temp_residual_rms;
    prec_residual_rms=sqrt(x'*beta_matrix_opt(55297:110592,:)'*beta_matrix_opt(55297:110592,:)*x);
    TvsP_prec_residual_rms_arr(floor(10*i+1))=prec_residual_rms;
end

figure()
hold on
scatter(TvsP_temp_residual_rms_arr,TvsP_prec_residual_rms_arr,'k')
plot(TvsP_temp_residual_rms_arr,TvsP_prec_residual_rms_arr,'k')

grid on
box on
xlabel('r.m.s. temperature change [°C]')
ylabel('r.m.s. precipitation change [mm/day]')
saveas(gcf,'C:\Users\ezrab\Documents\Cornell\Research\Figures\SAI_optimization\Fig_P_rms_vs_T_rms_different_weighting_3ens.png')

%% Perform optimization rms TREFHT vs rms PRECT (individual ensemble members)

for k = 1:length(varnames)
    eval(['alpha_' varnames{k} '_vecs = zeros(55296,n_base);'])
    eval(['mu_' varnames{k} '_arrs = zeros(55296,7,n_base);'])
    for i = 1:n_base %assuming all strategies have the same # of ensemble members
        eval(['alpha_' varnames{k} '_vecs(:,i) = reshape(alpha_' varnames{k} '(1:288,:,i),[55296 1]);'])
        eval(['mu_' varnames{k} '_arrs(:,:,i) = [reshape(mu_' varnames{k} '_EQ(1:288,:,i),[55296 1])' ...
            ' reshape(mu_' varnames{k} '_15N_15S(1:288,:,i),[55296 1])' ...
            ' reshape(mu_' varnames{k} '_30N_30S(1:288,:,i),[55296 1])' ...
            ' reshape(mu_' varnames{k} '_60N_60S(1:288,:,i),[55296 1])' ...
            ' reshape(mu_' varnames{k} '_60N(1:288,:,i),[55296 1])' ...
            ' reshape(mu_' varnames{k} '_def(1:288,:,i),[55296 1])' ...
            ' reshape(mu_' varnames{k} '_30N(1:288,:,i),[55296 1])];'])
    end
end

alpha_vecs = [alpha_TREFHT_vecs; alpha_PRECT_vecs];

Delta_arrs = [inj_EQ_end, inj_15N_15S_end, inj_30N_30S_end, inj_60N_60S_end, [12;12;12], sum(inj_def_end,2), [12;12;12]]...
        ./[cooling_EQ cooling_15N_15S cooling_30N_30S cooling_60N_60S cooling_60N cooling_def cooling_30N];
delta_Global_one_frac_per_ens = inj_def_end./sum(inj_def_end,2);

C_ens = zeros(7,7,3);
C_ens(:,:,1)=[0 0 0 0.5*Delta_arrs(1,4) Delta_arrs(1,5) 0 0 ;
       0 0 0.5*Delta_arrs(1,3) 0 0 delta_Global_one_frac_per_ens(1,4)*Delta_arrs(1,6) Delta_arrs(1,7) ;
       0 0.5*Delta_arrs(1,2) 0 0 0 delta_Global_one_frac_per_ens(1,3)*Delta_arrs(1,6) 0;
       Delta_arrs(1,1) 0 0 0 0 0 0;
       0 0.5*Delta_arrs(1,2) 0 0 0 delta_Global_one_frac_per_ens(1,2)*Delta_arrs(1,6) 0 ;
       0 0 0.5*Delta_arrs(1,3) 0 0 delta_Global_one_frac_per_ens(1,1)*Delta_arrs(1,6) 0;
       0 0 0 0.5*Delta_arrs(1,4) 0 0 0];
C_ens(:,:,2)=[0 0 0 0.5*Delta_arrs(2,4) Delta_arrs(2,5) 0 0 ;
       0 0 0.5*Delta_arrs(2,3) 0 0 delta_Global_one_frac_per_ens(2,4)*Delta_arrs(2,6) Delta_arrs(2,7) ;
       0 0.5*Delta_arrs(2,2) 0 0 0 delta_Global_one_frac_per_ens(2,3)*Delta_arrs(2,6) 0;
       Delta_arrs(2,1) 0 0 0 0 0 0;
       0 0.5*Delta_arrs(2,2) 0 0 0 delta_Global_one_frac_per_ens(2,2)*Delta_arrs(2,6) 0 ;
       0 0 0.5*Delta_arrs(2,3) 0 0 delta_Global_one_frac_per_ens(2,1)*Delta_arrs(2,6) 0;
       0 0 0 0.5*Delta_arrs(2,4) 0 0 0];
C_ens(:,:,3)=[0 0 0 0.5*Delta_arrs(3,4) Delta_arrs(3,5) 0 0 ;
       0 0 0.5*Delta_arrs(3,3) 0 0 delta_Global_one_frac_per_ens(3,4)*Delta_arrs(3,6) Delta_arrs(3,7) ;
       0 0.5*Delta_arrs(3,2) 0 0 0 delta_Global_one_frac_per_ens(3,3)*Delta_arrs(3,6) 0;
       Delta_arrs(3,1) 0 0 0 0 0 0;
       0 0.5*Delta_arrs(3,2) 0 0 0 delta_Global_one_frac_per_ens(3,2)*Delta_arrs(3,6) 0 ;
       0 0 0.5*Delta_arrs(3,3) 0 0 delta_Global_one_frac_per_ens(3,1)*Delta_arrs(3,6) 0;
       0 0 0 0.5*Delta_arrs(3,4) 0 0 0];

Delta_T0=1;
inj_each_lat_arrTvsP_001=zeros(7,11);
TvsP_x_001 = zeros(7,11);
inj_each_lat_arrTvsP_002=zeros(7,11);
TvsP_x_002 = zeros(7,11);
inj_each_lat_arrTvsP_003=zeros(7,11);
TvsP_x_003 = zeros(7,11);

TvsP_temp_residual_rms_arrs=zeros(3,11);
TvsP_prec_residual_rms_arrs=zeros(3,11);
TvsP_inj_each_lat_arrs = zeros(7,11,3);

for i=0:0.1:1
    % weighting for temperature and precipitation
    W_t=i*ones(55296,1);
    W_p=(1-i)*ones(55296,1);

    beta_matrix_001=[mu_TREFHT_arrs(:,:,1);mu_PRECT_arrs(:,:,1)]+alpha_vecs(:,1);
    beta_matrix_TvsP_001=beta_matrix_001.*sqrt(area_weight_vec.*[W_t;W_p]);
    beta_matrix_002=[mu_TREFHT_arrs(:,:,2);mu_PRECT_arrs(:,:,2)]+alpha_vecs(:,2);
    beta_matrix_TvsP_002=beta_matrix_002.*sqrt(area_weight_vec.*[W_t;W_p]);
    beta_matrix_003=[mu_TREFHT_arrs(:,:,3);mu_PRECT_arrs(:,:,3)]+alpha_vecs(:,3);
    beta_matrix_TvsP_003=beta_matrix_003.*sqrt(area_weight_vec.*[W_t;W_p]);
    % f=0.05*Delta_arr';
    f=zeros(7,1);
    % epsilon=0.001*ones(7,1);
    
    H_001=(beta_matrix_TvsP_001)'*(beta_matrix_TvsP_001);
    H_002=(beta_matrix_TvsP_002)'*(beta_matrix_TvsP_002);
    H_003=(beta_matrix_TvsP_003)'*(beta_matrix_TvsP_003);
    
    % fun=@(x) sqrt((beta_matrix*x)'*(beta_matrix*x));
    
    A_001=-C_ens(:,:,1);
    A_002=-C_ens(:,:,2);
    A_003=-C_ens(:,:,3);
    b=zeros(7,1);
    Aeq=ones(1,7);
    beq=Delta_T0;
    % b=epsilon;
    % x0=[0.1 0.1 0.3 0 0.3 0.1 0.1];
    lb=[];
    ub=[];
    % [x,fval]=fmincon(fun,x0,A,b,lb,ub)
    [x_001,fval,exitflag,output]=quadprog(H_001,f,A_001,b,Aeq,beq,[],[]);
    [x_002,fval,exitflag,output]=quadprog(H_002,f,A_002,b,Aeq,beq,[],[]);
    [x_003,fval,exitflag,output]=quadprog(H_003,f,A_003,b,Aeq,beq,[],[]);
    % save fraction of injection for each SAI strategy
    TvsP_x_001(:,floor(10*i+1)) = x_001;
    TvsP_x_002(:,floor(10*i+1)) = x_002;
    TvsP_x_003(:,floor(10*i+1)) = x_003;
    % injection rate for 60N,30N,15N,0N,15S,30S,60S
    inj_each_lat_001=C_ens(:,:,1)*x_001;
    TvsP_inj_each_lat_arrs(:,floor(10*i+1),1)=inj_each_lat_001;
    inj_each_lat_002=C_ens(:,:,2)*x_002;
    TvsP_inj_each_lat_arrs(:,floor(10*i+1),2)=inj_each_lat_002;
    inj_each_lat_003=C_ens(:,:,3)*x_003;
    TvsP_inj_each_lat_arrs(:,floor(10*i+1),3)=inj_each_lat_003;
    % climate impact
    beta_matrix_opt_001=beta_matrix_001.*sqrt(area_weight_vec);
    temp_residual_rms=sqrt(x_001'*beta_matrix_opt_001(1:55296,:)'*beta_matrix_opt_001(1:55296,:)*x_001);
    TvsP_temp_residual_rms_arrs(1,floor(10*i+1))=temp_residual_rms;
    prec_residual_rms=sqrt(x_001'*beta_matrix_opt_001(55297:110592,:)'*beta_matrix_opt_001(55297:110592,:)*x_001);
    TvsP_prec_residual_rms_arrs(1,floor(10*i+1))=prec_residual_rms;

    beta_matrix_opt_002=beta_matrix_002.*sqrt(area_weight_vec);
    temp_residual_rms=sqrt(x_002'*beta_matrix_opt_002(1:55296,:)'*beta_matrix_opt_002(1:55296,:)*x_002);
    TvsP_temp_residual_rms_arrs(2,floor(10*i+1))=temp_residual_rms;
    prec_residual_rms=sqrt(x_002'*beta_matrix_opt_002(55297:110592,:)'*beta_matrix_opt_002(55297:110592,:)*x_002);
    TvsP_prec_residual_rms_arrs(2,floor(10*i+1))=prec_residual_rms;

    beta_matrix_opt_003=beta_matrix_003.*sqrt(area_weight_vec);
    temp_residual_rms=sqrt(x_003'*beta_matrix_opt_003(1:55296,:)'*beta_matrix_opt_003(1:55296,:)*x_003);
    TvsP_temp_residual_rms_arrs(3,floor(10*i+1))=temp_residual_rms;
    prec_residual_rms=sqrt(x_003'*beta_matrix_opt_003(55297:110592,:)'*beta_matrix_opt_003(55297:110592,:)*x_003);
    TvsP_prec_residual_rms_arrs(3,floor(10*i+1))=prec_residual_rms;
end

%% maps of optimal temp and prec diff from ref period
%create optimal mus
% T optimal
mu_TREFHT_T_opt = opt_frac_TvsP_temp(1)*ea_mu_TREFHT_EQ + opt_frac_TvsP_temp(2)*ea_mu_TREFHT_15N_15S ...
+ opt_frac_TvsP_temp(3)*ea_mu_TREFHT_30N_30S + opt_frac_TvsP_temp(4)*ea_mu_TREFHT_60N_60S ...
+ opt_frac_TvsP_temp(5)*ea_mu_TREFHT_60N + opt_frac_TvsP_temp(6)*ea_mu_TREFHT_def + opt_frac_TvsP_temp(7)*ea_mu_TREFHT_30N;

mu_PRECT_T_opt = opt_frac_TvsP_temp(1)*ea_mu_PRECT_EQ + opt_frac_TvsP_temp(2)*ea_mu_PRECT_15N_15S ...
+ opt_frac_TvsP_temp(3)*ea_mu_PRECT_30N_30S + opt_frac_TvsP_temp(4)*ea_mu_PRECT_60N_60S ...
+ opt_frac_TvsP_temp(5)*ea_mu_PRECT_60N + opt_frac_TvsP_temp(6)*ea_mu_PRECT_def + opt_frac_TvsP_temp(7)*ea_mu_PRECT_30N;

% balanced
mu_TREFHT_bal_opt = TvsP_x(1,6)*ea_mu_TREFHT_EQ + TvsP_x(2,6)*ea_mu_TREFHT_15N_15S ...
+ TvsP_x(3,6)*ea_mu_TREFHT_30N_30S + TvsP_x(4,6)*ea_mu_TREFHT_60N_60S ...
+ TvsP_x(5,6)*ea_mu_TREFHT_60N + TvsP_x(6,6)*ea_mu_TREFHT_def + TvsP_x(7,6)*ea_mu_TREFHT_30N;

mu_PRECT_bal_opt = TvsP_x(1,6)*ea_mu_PRECT_EQ + TvsP_x(2,6)*ea_mu_PRECT_15N_15S ...
+ TvsP_x(3,6)*ea_mu_PRECT_30N_30S + TvsP_x(4,6)*ea_mu_PRECT_60N_60S ...
+ TvsP_x(5,6)*ea_mu_PRECT_60N + TvsP_x(6,6)*ea_mu_PRECT_def + TvsP_x(7,6)*ea_mu_PRECT_30N;

% P optimal
mu_TREFHT_P_opt = opt_frac_TvsP_prec(1)*ea_mu_TREFHT_EQ + opt_frac_TvsP_prec(2)*ea_mu_TREFHT_15N_15S ...
+ opt_frac_TvsP_prec(3)*ea_mu_TREFHT_30N_30S + opt_frac_TvsP_prec(4)*ea_mu_TREFHT_60N_60S ...
+ opt_frac_TvsP_prec(5)*ea_mu_TREFHT_60N + opt_frac_TvsP_prec(6)*ea_mu_TREFHT_def + opt_frac_TvsP_prec(7)*ea_mu_TREFHT_30N;

mu_PRECT_P_opt = opt_frac_TvsP_prec(1)*ea_mu_PRECT_EQ + opt_frac_TvsP_prec(2)*ea_mu_PRECT_15N_15S ...
+ opt_frac_TvsP_prec(3)*ea_mu_PRECT_30N_30S + opt_frac_TvsP_prec(4)*ea_mu_PRECT_60N_60S ...
+ opt_frac_TvsP_prec(5)*ea_mu_PRECT_60N + opt_frac_TvsP_prec(6)*ea_mu_PRECT_def + opt_frac_TvsP_prec(7)*ea_mu_PRECT_30N;

Tmin = -3;
Tmax = 3;
Pmin = -2;
Pmax = 2;
[T_contours, T_color_map] = custom_color_map(Tmin,Tmax,25,'*RdBu');
[P_contours, P_color_map] = custom_color_map(Pmin,Pmax,11,'BrBG');

figure('Renderer', 'painters', 'Position', [20 60 800 800])
clf

% T opt TREFHT
Z=abs((ea_alpha_TREFHT+mu_TREFHT_T_opt)./(2*[T_se; T_se(1,:)]));
[lon_null,lat_null]=find(Z<1);
lon_null=(lon_null-1)/287*360;
lat_null=(lat_null-1)/191*180-90;

ax1=subplot(3,2,1);
ax1.Position = [.07 .7 .42 .25];
hold on
set(gca,'FontSize',11)
worldmap([-90 90],[-180 180])
setm(ax1, 'meridianlabel', 'off', 'parallellabel', 'off')
box on
contourfm(lat,lon2,ea_alpha_TREFHT'+mu_TREFHT_T_opt',T_contours,'LineStyle','none')
scatterm(lat_null,lon_null,.2,'k','filled');
load coastlines
geoshow(coastlat,coastlon,'Color','k')
caxis([Tmin Tmax])
hold off
title('\DeltaT , Optimal T strategy');
colormap(ax1, T_color_map);
add_abc_tofig(1,0.04,-0.04,12);

% T opt PRECT
Z=abs((ea_alpha_PRECT+mu_PRECT_T_opt)./(2*[P_se; P_se(1,:)]));
[lon_null,lat_null]=find(Z<1);
lon_null=(lon_null-1)/287*360;
lat_null=(lat_null-1)/191*180-90;

ax2=subplot(3,2,2);
ax2.Position = [.53 .7 .42 .25];
hold on
set(gca,'FontSize',11)
worldmap([-90 90],[-180 180])
setm(ax2, 'meridianlabel', 'off', 'parallellabel', 'off')
box on
contourfm(lat,lon2,ea_alpha_PRECT'+mu_PRECT_T_opt',P_contours,'LineStyle','none')
scatterm(lat_null,lon_null,.2,'k','filled');
load coastlines
geoshow(coastlat,coastlon,'Color','k')
caxis([Pmin Pmax])
hold off
title('\DeltaP , Optimal T strategy');
colormap(ax2, P_color_map);
add_abc_tofig(2,0.04,-0.04,12);

% balanced TREFHT
Z=abs((ea_alpha_TREFHT+mu_TREFHT_bal_opt)./(2*[T_se; T_se(1,:)]));
[lon_null,lat_null]=find(Z<1);
lon_null=(lon_null-1)/287*360;
lat_null=(lat_null-1)/191*180-90;

ax3=subplot(3,2,3);
ax3.Position = [.07 .41 .42 .25];
hold on
set(gca,'FontSize',11)
worldmap([-90 90],[-180 180])
setm(ax3, 'meridianlabel', 'off', 'parallellabel', 'off')
box on
contourfm(lat,lon2,ea_alpha_TREFHT'+mu_TREFHT_bal_opt',T_contours,'LineStyle','none')
scatterm(lat_null,lon_null,.2,'k','filled');
load coastlines
geoshow(coastlat,coastlon,'Color','k')
caxis([Tmin Tmax])
hold off
title('\DeltaT , Balanced strategy');
colormap(ax3, T_color_map);
add_abc_tofig(3,0.04,-0.04,12);

% balanced PRECT
Z=abs((ea_alpha_PRECT+mu_PRECT_bal_opt)./(2*[P_se; P_se(1,:)]));
[lon_null,lat_null]=find(Z<1);
lon_null=(lon_null-1)/287*360;
lat_null=(lat_null-1)/191*180-90;

ax4=subplot(3,2,4);
ax4.Position = [.53 .41 .42 .25];
hold on
set(gca,'FontSize',11)
worldmap([-90 90],[-180 180])
setm(ax4, 'meridianlabel', 'off', 'parallellabel', 'off')
box on
contourfm(lat,lon2,ea_alpha_PRECT'+mu_PRECT_bal_opt',P_contours,'LineStyle','none')
scatterm(lat_null,lon_null,.2,'k','filled');
load coastlines
geoshow(coastlat,coastlon,'Color','k')
caxis([Pmin Pmax])
hold off
title('\DeltaP , Balanced strategy');
colormap(ax4, P_color_map);
add_abc_tofig(4,0.04,-0.04,12);

% P opt TREFHT
Z=abs((ea_alpha_TREFHT+mu_TREFHT_P_opt)./(2*[T_se; T_se(1,:)]));
[lon_null,lat_null]=find(Z<1);
lon_null=(lon_null-1)/287*360;
lat_null=(lat_null-1)/191*180-90;

ax5=subplot(3,2,5);
ax5.Position = [.07 .02 .42 .35];
hold on
set(gca,'FontSize',11)
worldmap([-90 90],[-180 180])
setm(ax5, 'meridianlabel', 'off', 'parallellabel', 'off')
box on
contourfm(lat,lon2,ea_alpha_TREFHT'+mu_TREFHT_P_opt',T_contours,'LineStyle','none')
scatterm(lat_null,lon_null,.2,'k','filled');
load coastlines
geoshow(coastlat,coastlon,'Color','k')
caxis([Tmin Tmax])
cb=colorbar;
ylabel(cb,'K','FontSize',10,'Rotation',0,'FontWeight','bold')
cb.Location = 'southoutside';
hold off
title('\DeltaT , Optimal P strategy');
colormap(ax5, T_color_map);
add_abc_tofig(5,0.04,-0.04,12);

% P opt PRECT
Z=abs((ea_alpha_PRECT+mu_PRECT_P_opt)./(2*[P_se; P_se(1,:)]));
[lon_null,lat_null]=find(Z<1);
lon_null=(lon_null-1)/287*360;
lat_null=(lat_null-1)/191*180-90;

ax6=subplot(3,2,6);
ax6.Position = [.53 .02 .42 .35];
hold on
set(gca,'FontSize',11)
worldmap([-90 90],[-180 180])
setm(ax6, 'meridianlabel', 'off', 'parallellabel', 'off')
box on
contourfm(lat,lon2,ea_alpha_PRECT'+mu_PRECT_P_opt',P_contours,'LineStyle','none')
scatterm(lat_null,lon_null,.2,'k','filled');
load coastlines
geoshow(coastlat,coastlon,'Color','k')
caxis([Pmin Pmax])
cb=colorbar;
ylabel(cb,'mm/day','FontSize',12,'Rotation',0,'FontWeight','bold')
cb.Location = 'southoutside';
hold off
title('\DeltaP , Optimal P strategy');
colormap(ax6, P_color_map);
add_abc_tofig(6,0.04,-0.04,12);

%% Perform optimization rms TREFHT vs rms PRECT normalized per standard deviation (3 ens)

alpha_vec = [alpha_TREFHT_vec; alpha_PRECT_vec];

Delta_T0=1;
inj_each_lat_arrTvsP_norm=zeros(7,11);
TvsP_norm_x = zeros(7,11);
TvsP_norm_temp_residual_rms_arr=zeros(1,11);
TvsP_norm_prec_residual_rms_arr=zeros(1,11);
% temp_residual_rms_arr_land=zeros(1,11);
% prec_residual_rms_arr_land=zeros(1,11);

for i=0:0.1:1
    % weighting for temperature and precipitation
    W_t=i*ones(55296,1);
    W_p=(1-i)*ones(55296,1);

    beta_matrix=([mu_TREFHT_arr;mu_PRECT_arr]+alpha_vec)./([reshape(T_std,[55296 1]); reshape(P_std,[55296 1])]);
    beta_matrix_TvsP_norm=beta_matrix.*sqrt(area_weight_vec.*[W_t;W_p]);
    % f=0.05*Delta_arr';
    f=zeros(7,1);
    % epsilon=0.001*ones(7,1);
    
    H=(beta_matrix_TvsP_norm)'*(beta_matrix_TvsP_norm);
    
    % fun=@(x) sqrt((beta_matrix*x)'*(beta_matrix*x));
    
    A=-C;
    b=zeros(7,1);
    Aeq=ones(1,7);
    beq=Delta_T0;
    % b=epsilon;
    % x0=[0.1 0.1 0.3 0 0.3 0.1 0.1];
    lb=[];
    ub=[];
    % [x,fval]=fmincon(fun,x0,A,b,lb,ub)
    [x,fval,exitflag,output]=quadprog(H,f,A,b,Aeq,beq,[],[]);
    TvsP_norm_x(:,floor(10*i+1)) = x;
    % injection rate for 60N,30N,15N,0N,15S,30S,60S
    inj_each_lat=C*x;
    TvsP_norm_inj_each_lat_arr(:,floor(10*i+1))=inj_each_lat;
    % climate impact
    beta_matrix_opt=beta_matrix.*sqrt(area_weight_vec);
    temp_residual_rms=sqrt(x'*beta_matrix_opt(1:55296,:)'*beta_matrix_opt(1:55296,:)*x);
    TvsP_norm_temp_residual_rms_arr(floor(10*i+1))=temp_residual_rms;
    prec_residual_rms=sqrt(x'*beta_matrix_opt(55297:110592,:)'*beta_matrix_opt(55297:110592,:)*x);
    TvsP_norm_prec_residual_rms_arr(floor(10*i+1))=prec_residual_rms;

    % land only
%     beta_matrix_land_opt=beta_matrix.*sqrt(land_weight_vec);
%     temp_residual_rms_land=sqrt(x'*beta_matrix_land_opt(1:55296,:)'*beta_matrix_land_opt(1:55296,:)*x);
%     temp_residual_rms_arr_land(floor(10*i+1))=temp_residual_rms_land;
%     prec_residual_rms_land=sqrt(x'*beta_matrix_land_opt(55297:110592,:)'*beta_matrix_land_opt(55297:110592,:)*x);
%     prec_residual_rms_arr_land(floor(10*i+1))=prec_residual_rms_land;
end

% save injection rates and r.m.s.
writetable(table(TvsP_norm_inj_each_lat_arr),'injection_rates_different_weighting.csv');
writetable(table([TvsP_norm_temp_residual_rms_arr;TvsP_norm_prec_residual_rms_arr]),'residual_different_weighting.csv')

figure()
hold on
scatter(TvsP_norm_temp_residual_rms_arr,TvsP_norm_prec_residual_rms_arr,'k')
plot(TvsP_norm_temp_residual_rms_arr,TvsP_norm_prec_residual_rms_arr,'k')

grid on
box on
xlabel('r.m.s. temperature change [\sigma]')
ylabel('r.m.s. precipitation change [\sigma]')
saveas(gcf,'C:\Users\ezrab\Documents\Cornell\Research\Figures\SAI_optimization\Fig_P_rms_vs_T_rms_norm_different_weighting_3ens.png')

%% maps of optimal temp and prec diff from ref period, normalized per std deviation
%create optimal mus
% T optimal
mu_TREFHT_T_norm_opt = TvsP_norm_x(1,11)*ea_mu_TREFHT_EQ + TvsP_norm_x(2,11)*ea_mu_TREFHT_15N_15S ...
+ TvsP_norm_x(3,11)*ea_mu_TREFHT_30N_30S + TvsP_norm_x(4,11)*ea_mu_TREFHT_60N_60S ...
+ TvsP_norm_x(5,11)*ea_mu_TREFHT_60N + TvsP_norm_x(6,11)*ea_mu_TREFHT_def + TvsP_norm_x(7,11)*ea_mu_TREFHT_30N;

mu_PRECT_T_norm_opt = TvsP_norm_x(1,11)*ea_mu_PRECT_EQ + TvsP_norm_x(2,11)*ea_mu_PRECT_15N_15S ...
+ TvsP_norm_x(3,11)*ea_mu_PRECT_30N_30S + TvsP_norm_x(4,11)*ea_mu_PRECT_60N_60S ...
+ TvsP_norm_x(5,11)*ea_mu_PRECT_60N + TvsP_norm_x(6,11)*ea_mu_PRECT_def + TvsP_norm_x(7,11)*ea_mu_PRECT_30N;

% balanced
mu_TREFHT_bal_norm_opt = TvsP_norm_x(1,6)*ea_mu_TREFHT_EQ + TvsP_norm_x(2,6)*ea_mu_TREFHT_15N_15S ...
+ TvsP_norm_x(3,6)*ea_mu_TREFHT_30N_30S + TvsP_norm_x(4,6)*ea_mu_TREFHT_60N_60S ...
+ TvsP_norm_x(5,6)*ea_mu_TREFHT_60N + TvsP_norm_x(6,6)*ea_mu_TREFHT_def + TvsP_norm_x(7,6)*ea_mu_TREFHT_30N;

mu_PRECT_bal_norm_opt = TvsP_norm_x(1,6)*ea_mu_PRECT_EQ + TvsP_norm_x(2,6)*ea_mu_PRECT_15N_15S ...
+ TvsP_norm_x(3,6)*ea_mu_PRECT_30N_30S + TvsP_norm_x(4,6)*ea_mu_PRECT_60N_60S ...
+ TvsP_norm_x(5,6)*ea_mu_PRECT_60N + TvsP_norm_x(6,6)*ea_mu_PRECT_def + TvsP_norm_x(7,6)*ea_mu_PRECT_30N;

% P optimal
mu_TREFHT_P_norm_opt = TvsP_norm_x(1,1)*ea_mu_TREFHT_EQ + TvsP_norm_x(2,1)*ea_mu_TREFHT_15N_15S ...
+ TvsP_norm_x(3,1)*ea_mu_TREFHT_30N_30S + TvsP_norm_x(4,1)*ea_mu_TREFHT_60N_60S ...
+ TvsP_norm_x(5,1)*ea_mu_TREFHT_60N + TvsP_norm_x(6,1)*ea_mu_TREFHT_def + TvsP_norm_x(7,1)*ea_mu_TREFHT_30N;

mu_PRECT_P_norm_opt = TvsP_norm_x(1,1)*ea_mu_PRECT_EQ + TvsP_norm_x(2,1)*ea_mu_PRECT_15N_15S ...
+ TvsP_norm_x(3,1)*ea_mu_PRECT_30N_30S + TvsP_norm_x(4,1)*ea_mu_PRECT_60N_60S ...
+ TvsP_norm_x(5,1)*ea_mu_PRECT_60N + TvsP_norm_x(6,1)*ea_mu_PRECT_def + TvsP_norm_x(7,1)*ea_mu_PRECT_30N;

Tmin = -4;
Tmax = 4;
Pmin = -3;
Pmax = 3;
[T_contours, T_color_map] = custom_color_map(Tmin,Tmax,25,'*RdBu');
[P_contours, P_color_map] = custom_color_map(Pmin,Pmax,11,'BrBG');

figure('Renderer', 'painters', 'Position', [20 60 800 800])
clf

% T opt TREFHT
Z=abs((ea_alpha_TREFHT+mu_TREFHT_T_norm_opt)./(2*[T_se; T_se(1,:)]));
[lon_null,lat_null]=find(Z<1);
lon_null=(lon_null-1)/287*360;
lat_null=(lat_null-1)/191*180-90;

ax1=subplot(3,2,1);
ax1.Position = [.07 .7 .42 .25];
hold on
set(gca,'FontSize',11)
worldmap([-90 90],[-180 180])
setm(ax1, 'meridianlabel', 'off', 'parallellabel', 'off')
box on
contourfm(lat,lon2,((ea_alpha_TREFHT+mu_TREFHT_T_norm_opt)./[T_std; T_std(1,:)])', T_contours,'LineStyle','none');
scatterm(lat_null,lon_null,.2,'k','filled');
load coastlines
geoshow(coastlat,coastlon,'Color','k')
caxis([Tmin Tmax])
hold off
title('\DeltaT , Optimal T strategy');
colormap(ax1, T_color_map);
add_abc_tofig(1,0.04,-0.04,12);

% T opt PRECT
Z=abs((ea_alpha_PRECT+mu_PRECT_T_norm_opt)./(2*[P_se; P_se(1,:)]));
[lon_null,lat_null]=find(Z<1);
lon_null=(lon_null-1)/287*360;
lat_null=(lat_null-1)/191*180-90;

ax2=subplot(3,2,2);
ax2.Position = [.53 .7 .42 .25];
hold on
set(gca,'FontSize',11)
worldmap([-90 90],[-180 180])
setm(ax2, 'meridianlabel', 'off', 'parallellabel', 'off')
box on
contourfm(lat,lon2,((ea_alpha_PRECT+mu_PRECT_T_norm_opt)./[P_std; P_std(1,:)])',P_contours,'LineStyle','none')
scatterm(lat_null,lon_null,.2,'k','filled');
load coastlines
geoshow(coastlat,coastlon,'Color','k')
caxis([Pmin Pmax])
hold off
title('\DeltaP , Optimal T strategy');
colormap(ax2, P_color_map);
add_abc_tofig(2,0.04,-0.04,12);

% balanced TREFHT
Z=abs((ea_alpha_TREFHT+mu_TREFHT_bal_norm_opt)./(2*[T_se; T_se(1,:)]));
[lon_null,lat_null]=find(Z<1);
lon_null=(lon_null-1)/287*360;
lat_null=(lat_null-1)/191*180-90;

ax3=subplot(3,2,3);
ax3.Position = [.07 .41 .42 .25];
hold on
set(gca,'FontSize',11)
worldmap([-90 90],[-180 180])
setm(ax3, 'meridianlabel', 'off', 'parallellabel', 'off')
box on
contourfm(lat,lon2,((ea_alpha_TREFHT+mu_TREFHT_bal_norm_opt)./[T_std; T_std(1,:)])',T_contours,'LineStyle','none')
scatterm(lat_null,lon_null,.2,'k','filled');
load coastlines
geoshow(coastlat,coastlon,'Color','k')
caxis([Tmin Tmax])
hold off
title('\DeltaT , Balanced strategy');
colormap(ax3, T_color_map);
add_abc_tofig(3,0.04,-0.04,12);

% balanced PRECT
Z=abs((ea_alpha_PRECT+mu_PRECT_bal_norm_opt)./(2*[P_se; P_se(1,:)]));
[lon_null,lat_null]=find(Z<1);
lon_null=(lon_null-1)/287*360;
lat_null=(lat_null-1)/191*180-90;

ax4=subplot(3,2,4);
ax4.Position = [.53 .41 .42 .25];
hold on
set(gca,'FontSize',11)
worldmap([-90 90],[-180 180])
setm(ax4, 'meridianlabel', 'off', 'parallellabel', 'off')
box on
contourfm(lat,lon2,((ea_alpha_PRECT+mu_PRECT_bal_norm_opt)./[P_std; P_std(1,:)])',P_contours,'LineStyle','none')
scatterm(lat_null,lon_null,.2,'k','filled');
load coastlines
geoshow(coastlat,coastlon,'Color','k')
caxis([Pmin Pmax])
hold off
title('\DeltaP , Balanced strategy');
colormap(ax4, P_color_map);
add_abc_tofig(4,0.04,-0.04,12);

% P opt TREFHT
Z=abs((ea_alpha_TREFHT+mu_TREFHT_P_norm_opt)./(2*[T_se; T_se(1,:)]));
[lon_null,lat_null]=find(Z<1);
lon_null=(lon_null-1)/287*360;
lat_null=(lat_null-1)/191*180-90;

ax5=subplot(3,2,5);
ax5.Position = [.07 .02 .42 .35];
hold on
set(gca,'FontSize',11)
worldmap([-90 90],[-180 180])
setm(ax5, 'meridianlabel', 'off', 'parallellabel', 'off')
box on
contourfm(lat,lon2,((ea_alpha_TREFHT+mu_TREFHT_P_norm_opt)./[T_std; T_std(1,:)])',T_contours,'LineStyle','none')
scatterm(lat_null,lon_null,.2,'k','filled');
load coastlines
geoshow(coastlat,coastlon,'Color','k')
caxis([Tmin Tmax])
cb=colorbar;
ylabel(cb,'\sigma','FontSize',10,'Rotation',0,'FontWeight','bold')
cb.Location = 'southoutside';
hold off
title('\DeltaT , Optimal P strategy');
colormap(ax5, T_color_map);
add_abc_tofig(5,0.04,-0.04,12);

% P opt PRECT
Z=abs((ea_alpha_PRECT+mu_PRECT_P_norm_opt)./(2*[P_se; P_se(1,:)]));
[lon_null,lat_null]=find(Z<1);
lon_null=(lon_null-1)/287*360;
lat_null=(lat_null-1)/191*180-90;

ax6=subplot(3,2,6);
ax6.Position = [.53 .02 .42 .35];
hold on
set(gca,'FontSize',11)
worldmap([-90 90],[-180 180])
setm(ax6, 'meridianlabel', 'off', 'parallellabel', 'off')
box on
contourfm(lat,lon2,((ea_alpha_PRECT+mu_PRECT_P_norm_opt)./[P_std; P_std(1,:)])',P_contours,'LineStyle','none')
scatterm(lat_null,lon_null,.2,'k','filled');
load coastlines
geoshow(coastlat,coastlon,'Color','k')
caxis([Pmin Pmax])
cb=colorbar;
ylabel(cb,'\sigma','FontSize',12,'Rotation',0,'FontWeight','bold')
cb.Location = 'southoutside';
hold off
title('\DeltaP , Optimal P strategy');
colormap(ax6, P_color_map);
add_abc_tofig(6,0.04,-0.04,12);

%% maps of optimal temp and prec diff from ref period, normalized per std deviation, shown in absolute units
%create optimal mus
% T optimal
mu_TREFHT_T_norm_opt = TvsP_norm_x(1,11)*ea_mu_TREFHT_EQ + TvsP_norm_x(2,11)*ea_mu_TREFHT_15N_15S ...
+ TvsP_norm_x(3,11)*ea_mu_TREFHT_30N_30S + TvsP_norm_x(4,11)*ea_mu_TREFHT_60N_60S ...
+ TvsP_norm_x(5,11)*ea_mu_TREFHT_60N + TvsP_norm_x(6,11)*ea_mu_TREFHT_def + TvsP_norm_x(7,11)*ea_mu_TREFHT_30N;

mu_PRECT_T_norm_opt = TvsP_norm_x(1,11)*ea_mu_PRECT_EQ + TvsP_norm_x(2,11)*ea_mu_PRECT_15N_15S ...
+ TvsP_norm_x(3,11)*ea_mu_PRECT_30N_30S + TvsP_norm_x(4,11)*ea_mu_PRECT_60N_60S ...
+ TvsP_norm_x(5,11)*ea_mu_PRECT_60N + TvsP_norm_x(6,11)*ea_mu_PRECT_def + TvsP_norm_x(7,11)*ea_mu_PRECT_30N;

% balanced
mu_TREFHT_bal_norm_opt = TvsP_norm_x(1,6)*ea_mu_TREFHT_EQ + TvsP_norm_x(2,6)*ea_mu_TREFHT_15N_15S ...
+ TvsP_norm_x(3,6)*ea_mu_TREFHT_30N_30S + TvsP_norm_x(4,6)*ea_mu_TREFHT_60N_60S ...
+ TvsP_norm_x(5,6)*ea_mu_TREFHT_60N + TvsP_norm_x(6,6)*ea_mu_TREFHT_def + TvsP_norm_x(7,6)*ea_mu_TREFHT_30N;

mu_PRECT_bal_norm_opt = TvsP_norm_x(1,6)*ea_mu_PRECT_EQ + TvsP_norm_x(2,6)*ea_mu_PRECT_15N_15S ...
+ TvsP_norm_x(3,6)*ea_mu_PRECT_30N_30S + TvsP_norm_x(4,6)*ea_mu_PRECT_60N_60S ...
+ TvsP_norm_x(5,6)*ea_mu_PRECT_60N + TvsP_norm_x(6,6)*ea_mu_PRECT_def + TvsP_norm_x(7,6)*ea_mu_PRECT_30N;

% P optimal
mu_TREFHT_P_norm_opt = TvsP_norm_x(1,1)*ea_mu_TREFHT_EQ + TvsP_norm_x(2,1)*ea_mu_TREFHT_15N_15S ...
+ TvsP_norm_x(3,1)*ea_mu_TREFHT_30N_30S + TvsP_norm_x(4,1)*ea_mu_TREFHT_60N_60S ...
+ TvsP_norm_x(5,1)*ea_mu_TREFHT_60N + TvsP_norm_x(6,1)*ea_mu_TREFHT_def + TvsP_norm_x(7,1)*ea_mu_TREFHT_30N;

mu_PRECT_P_norm_opt = TvsP_norm_x(1,1)*ea_mu_PRECT_EQ + TvsP_norm_x(2,1)*ea_mu_PRECT_15N_15S ...
+ TvsP_norm_x(3,1)*ea_mu_PRECT_30N_30S + TvsP_norm_x(4,1)*ea_mu_PRECT_60N_60S ...
+ TvsP_norm_x(5,1)*ea_mu_PRECT_60N + TvsP_norm_x(6,1)*ea_mu_PRECT_def + TvsP_norm_x(7,1)*ea_mu_PRECT_30N;

Tmin = -3;
Tmax = 3;
Pmin = -2;
Pmax = 2;
[T_contours, T_color_map] = custom_color_map(Tmin,Tmax,25,'*RdBu');
[P_contours, P_color_map] = custom_color_map(Pmin,Pmax,11,'BrBG');

figure('Renderer', 'painters', 'Position', [20 60 800 800])
clf

% T opt TREFHT
Z=abs((ea_alpha_TREFHT+mu_TREFHT_T_norm_opt)./(2*[T_se; T_se(1,:)]));
[lon_null,lat_null]=find(Z<1);
lon_null=(lon_null-1)/287*360;
lat_null=(lat_null-1)/191*180-90;

ax1=subplot(3,2,1);
ax1.Position = [.07 .7 .42 .25];
hold on
set(gca,'FontSize',11)
worldmap([-90 90],[-180 180])
setm(ax1, 'meridianlabel', 'off', 'parallellabel', 'off')
box on
contourfm(lat,lon2,(ea_alpha_TREFHT+mu_TREFHT_T_norm_opt)', T_contours,'LineStyle','none');
scatterm(lat_null,lon_null,.2,'k','filled');
load coastlines
geoshow(coastlat,coastlon,'Color','k')
caxis([Tmin Tmax])
hold off
title('\DeltaT , Optimal T strategy');
colormap(ax1, T_color_map);
add_abc_tofig(1,0.04,-0.04,12);

% T opt PRECT
Z=abs((ea_alpha_PRECT+mu_PRECT_T_norm_opt)./(2*[P_se; P_se(1,:)]));
[lon_null,lat_null]=find(Z<1);
lon_null=(lon_null-1)/287*360;
lat_null=(lat_null-1)/191*180-90;

ax2=subplot(3,2,2);
ax2.Position = [.53 .7 .42 .25];
hold on
set(gca,'FontSize',11)
worldmap([-90 90],[-180 180])
setm(ax2, 'meridianlabel', 'off', 'parallellabel', 'off')
box on
contourfm(lat,lon2,(ea_alpha_PRECT+mu_PRECT_T_norm_opt)',P_contours,'LineStyle','none')
scatterm(lat_null,lon_null,.2,'k','filled');
load coastlines
geoshow(coastlat,coastlon,'Color','k')
caxis([Pmin Pmax])
hold off
title('\DeltaP , Optimal T strategy');
colormap(ax2, P_color_map);
add_abc_tofig(2,0.04,-0.04,12);

% balanced TREFHT
Z=abs((ea_alpha_TREFHT+mu_TREFHT_bal_norm_opt)./(2*[T_se; T_se(1,:)]));
[lon_null,lat_null]=find(Z<1);
lon_null=(lon_null-1)/287*360;
lat_null=(lat_null-1)/191*180-90;

ax3=subplot(3,2,3);
ax3.Position = [.07 .41 .42 .25];
hold on
set(gca,'FontSize',11)
worldmap([-90 90],[-180 180])
setm(ax3, 'meridianlabel', 'off', 'parallellabel', 'off')
box on
contourfm(lat,lon2,(ea_alpha_TREFHT+mu_TREFHT_bal_norm_opt)',T_contours,'LineStyle','none')
scatterm(lat_null,lon_null,.2,'k','filled');
load coastlines
geoshow(coastlat,coastlon,'Color','k')
caxis([Tmin Tmax])
hold off
title('\DeltaT , Balanced strategy');
colormap(ax3, T_color_map);
add_abc_tofig(3,0.04,-0.04,12);

% balanced PRECT
Z=abs((ea_alpha_PRECT+mu_PRECT_bal_norm_opt)./(2*[P_se; P_se(1,:)]));
[lon_null,lat_null]=find(Z<1);
lon_null=(lon_null-1)/287*360;
lat_null=(lat_null-1)/191*180-90;

ax4=subplot(3,2,4);
ax4.Position = [.53 .41 .42 .25];
hold on
set(gca,'FontSize',11)
worldmap([-90 90],[-180 180])
setm(ax4, 'meridianlabel', 'off', 'parallellabel', 'off')
box on
contourfm(lat,lon2,(ea_alpha_PRECT+mu_PRECT_bal_norm_opt)',P_contours,'LineStyle','none')
scatterm(lat_null,lon_null,.2,'k','filled');
load coastlines
geoshow(coastlat,coastlon,'Color','k')
caxis([Pmin Pmax])
hold off
title('\DeltaP , Balanced strategy');
colormap(ax4, P_color_map);
add_abc_tofig(4,0.04,-0.04,12);

% P opt TREFHT
Z=abs((ea_alpha_TREFHT+mu_TREFHT_P_norm_opt)./(2*[T_se; T_se(1,:)]));
[lon_null,lat_null]=find(Z<1);
lon_null=(lon_null-1)/287*360;
lat_null=(lat_null-1)/191*180-90;

ax5=subplot(3,2,5);
ax5.Position = [.07 .02 .42 .35];
hold on
set(gca,'FontSize',11)
worldmap([-90 90],[-180 180])
setm(ax5, 'meridianlabel', 'off', 'parallellabel', 'off')
box on
contourfm(lat,lon2,(ea_alpha_TREFHT+mu_TREFHT_P_norm_opt)',T_contours,'LineStyle','none')
scatterm(lat_null,lon_null,.2,'k','filled');
load coastlines
geoshow(coastlat,coastlon,'Color','k')
caxis([Tmin Tmax])
cb=colorbar;
ylabel(cb,'^\circC','FontSize',10,'Rotation',0,'FontWeight','bold')
cb.Location = 'southoutside';
hold off
title('\DeltaT , Optimal P strategy');
colormap(ax5, T_color_map);
add_abc_tofig(5,0.04,-0.04,12);

% P opt PRECT
Z=abs((ea_alpha_PRECT+mu_PRECT_P_norm_opt)./(2*[P_se; P_se(1,:)]));
[lon_null,lat_null]=find(Z<1);
lon_null=(lon_null-1)/287*360;
lat_null=(lat_null-1)/191*180-90;

ax6=subplot(3,2,6);
ax6.Position = [.53 .02 .42 .35];
hold on
set(gca,'FontSize',11)
worldmap([-90 90],[-180 180])
setm(ax6, 'meridianlabel', 'off', 'parallellabel', 'off')
box on
contourfm(lat,lon2,(ea_alpha_PRECT+mu_PRECT_P_norm_opt)',P_contours,'LineStyle','none')
scatterm(lat_null,lon_null,.2,'k','filled');
load coastlines
geoshow(coastlat,coastlon,'Color','k')
caxis([Pmin Pmax])
cb=colorbar;
ylabel(cb,'mm/day','FontSize',12,'Rotation',0,'FontWeight','bold')
cb.Location = 'southoutside';
hold off
title('\DeltaP , Optimal P strategy');
colormap(ax6, P_color_map);
add_abc_tofig(6,0.04,-0.04,12);

%% Perform optimization rms TREFHT vs rms PRECT (3 ens, no EQ)

alpha_vec = [alpha_TREFHT_vec; alpha_PRECT_vec];

Delta_T0=1;
inj_each_lat_arr_noEQ=zeros(6,11);
temp_residual_rms_arr_noEQ=zeros(1,11);
prec_residual_rms_arr_noEQ=zeros(1,11);

C_noEQ=[0 0 0.5*Delta_arr(4) Delta_arr(5) 0 0 ;
       0 0.5*Delta_arr(3) 0 0 delta_Global_one_frac(4)*Delta_arr(6) Delta_arr(7) ;
       0.5*Delta_arr(2) 0 0 0 delta_Global_one_frac(3)*Delta_arr(6) 0;
       0.5*Delta_arr(2) 0 0 0 delta_Global_one_frac(2)*Delta_arr(6) 0 ;
       0 0.5*Delta_arr(3) 0 0 delta_Global_one_frac(1)*Delta_arr(6) 0;
       0 0 0.5*Delta_arr(4) 0 0 0];

for i=0:0.1:1
    % weighting for temperature and precipitation
    W_t=i*ones(55296,1);
    W_p=(1-i)*ones(55296,1);

    beta_matrix=[mu_TREFHT_arr(:,2:7);mu_PRECT_arr(:,2:7)]+alpha_vec;
    beta_matrix_w=beta_matrix.*sqrt(area_weight_vec.*[W_t;W_p]);
    % f=0.05*Delta_arr';
    f=zeros(6,1);
    % epsilon=0.001*ones(7,1);
    
    H=(beta_matrix_w)'*(beta_matrix_w);
    
    % fun=@(x) sqrt((beta_matrix*x)'*(beta_matrix*x));
    
    A=-C_noEQ;
    b=zeros(6,1);
    Aeq=ones(1,6);
    beq=Delta_T0;
    % b=epsilon;
    % x0=[0.1 0.1 0.3 0 0.3 0.1 0.1];
    lb=[];
    ub=[];
    % [x,fval]=fmincon(fun,x0,A,b,lb,ub)
    [x,fval,exitflag,output]=quadprog(H,f,A,b,Aeq,beq,[],[]);
    % injection rate for 60N,30N,15N,0N,15S,30S,60S
    inj_each_lat=C_noEQ*x;
    inj_each_lat_arr_noEQ(:,floor(10*i+1))=inj_each_lat;
    % climate impact
    beta_matrix_opt=beta_matrix.*sqrt(area_weight_vec);
    temp_residual_rms_noEQ=sqrt(x'*beta_matrix_opt(1:55296,:)'*beta_matrix_opt(1:55296,:)*x);
    temp_residual_rms_arr_noEQ(floor(10*i+1))=temp_residual_rms_noEQ;
    prec_residual_rms_noEQ=sqrt(x'*beta_matrix_opt(55297:110592,:)'*beta_matrix_opt(55297:110592,:)*x);
    prec_residual_rms_arr_noEQ(floor(10*i+1))=prec_residual_rms_noEQ;
end

% save injection rates and r.m.s.
writetable(table(inj_each_lat_arr_noEQ),'injection_rates_different_weighting_noEQ.csv');
writetable(table([temp_residual_rms_arr_noEQ;prec_residual_rms_arr_noEQ]),'residual_different_weighting_noEQ.csv')

figure()
hold on
scatter(temp_residual_rms_arr_noEQ,prec_residual_rms_arr_noEQ,'k')
plot(temp_residual_rms_arr_noEQ,prec_residual_rms_arr_noEQ,'k')

grid on
box on
xlabel('r.m.s. temperature change [°C]')
ylabel('r.m.s. precipitation change [mm/day]')
saveas(gcf,'C:\Users\ezrab\Documents\Cornell\Research\Figures\SAI_optimization\Fig_P_rms_vs_T_rms_different_weighting_3ens_noEQ.png')

%% Perform optimization rms land TREFHT vs rms land PRECT-QFLX (3 ens)

alpha_vec = [alpha_TREFHT_vec; alpha_PRECT_vec-alpha_QFLX_vec];

Delta_T0=1;
landT_landPE_inj_each_lat_arr=zeros(7,11);
landTvslandPE_x = zeros(7,11);
landT_landPE_temp_residual_rms_arr=zeros(1,11);
landT_landPE_PE_residual_rms_arr=zeros(1,11);
% temp_residual_rms_arr_land=zeros(1,11);
% prec_residual_rms_arr_land=zeros(1,11);

for i=0:0.1:1
    % weighting for temperature and precipitation
    W_t=i*ones(55296,1);
    W_p=(1-i)*ones(55296,1);

    beta_matrix=[mu_TREFHT_arr;mu_PRECT_arr-mu_QFLX_arr]+alpha_vec;
%     beta_matrix=beta_matrix*Delta_T0;
    beta_matrix_landT_landPE=beta_matrix.*sqrt(land_weight_vec.*[W_t;W_p]);
    % f=0.05*Delta_arr';
    f=zeros(7,1);
    % epsilon=0.001*ones(7,1);
    
    H=(beta_matrix_landT_landPE)'*(beta_matrix_landT_landPE);
    
    % fun=@(x) sqrt((beta_matrix*x)'*(beta_matrix*x));
    
    A=-C;
    b=zeros(7,1);
    Aeq=ones(1,7);
    beq=Delta_T0;
    % b=epsilon;
    % x0=[0.1 0.1 0.3 0 0.3 0.1 0.1];
    lb=[];
    ub=[];
    % [x,fval]=fmincon(fun,x0,A,b,lb,ub)
    [x,fval,exitflag,output]=quadprog(H,f,A,b,Aeq,beq,[],[]);
    landTvslandPE_x(:,floor(10*i+1)) = x;
   % injection rate for 60N,30N,15N,0N,15S,30S,60S
    landT_landPE_inj_each_lat=C*x;
    landT_landPE_inj_each_lat_arr(:,floor(10*i+1))=landT_landPE_inj_each_lat;

    % land only
    beta_matrix_land_opt=beta_matrix.*sqrt(land_weight_vec);
    landT_landPE_temp_residual_rms=sqrt(x'*beta_matrix_land_opt(1:55296,:)'*beta_matrix_land_opt(1:55296,:)*x);
    landT_landPE_temp_residual_rms_arr(floor(10*i+1))=landT_landPE_temp_residual_rms;
    landT_landPE_PE_residual_rms=sqrt(x'*beta_matrix_land_opt(55297:110592,:)'*beta_matrix_land_opt(55297:110592,:)*x);
    landT_landPE_PE_residual_rms_arr(floor(10*i+1))=landT_landPE_PE_residual_rms;
end

% save injection rates and r.m.s.
writetable(table(landT_landPE_inj_each_lat_arr),'landT_landPE_injection_rates_different_weighting.csv');
writetable(table([landT_landPE_temp_residual_rms_arr;landT_landPE_PE_residual_rms_arr]),'landT_landPE_residual_different_weighting.csv')

figure()
hold on
scatter(landT_landPE_temp_residual_rms_arr,landT_landPE_PE_residual_rms_arr,'k')
plot(landT_landPE_temp_residual_rms_arr,landT_landPE_PE_residual_rms_arr,'k')

grid on
box on
xlabel('r.m.s. temperature change [°C]')
ylabel('r.m.s. P-E change [mm/day]')
saveas(gcf,'C:\Users\ezrab\Documents\Cornell\Research\Figures\SAI_optimization\Fig_landPE_rms_vs_landT_rms_different_weighting_3ens.png')

%% maps of optimal landT and landP-E diff from ref period
%create optimal mus
% land T optimal
mu_TREFHT_landT_opt = landTvslandPE_x(1,11)*ea_mu_TREFHT_EQ + landTvslandPE_x(2,11)*ea_mu_TREFHT_15N_15S ...
+ landTvslandPE_x(3,11)*ea_mu_TREFHT_30N_30S + landTvslandPE_x(4,11)*ea_mu_TREFHT_60N_60S ...
+ landTvslandPE_x(5,11)*ea_mu_TREFHT_60N + landTvslandPE_x(6,11)*ea_mu_TREFHT_def + landTvslandPE_x(7,11)*ea_mu_TREFHT_30N;

mu_PE_landT_opt = landTvslandPE_x(1,11)*(ea_mu_PRECT_EQ-ea_mu_QFLX_EQ) + landTvslandPE_x(2,11)*(ea_mu_PRECT_15N_15S-ea_mu_QFLX_15N_15S) ...
+ landTvslandPE_x(3,11)*(ea_mu_PRECT_30N_30S-ea_mu_QFLX_30N_30S) + landTvslandPE_x(4,11)*(ea_mu_PRECT_60N_60S-ea_mu_QFLX_60N_60S) ...
+ landTvslandPE_x(5,11)*(ea_mu_PRECT_60N-ea_mu_QFLX_60N) + landTvslandPE_x(6,11)*(ea_mu_PRECT_def-ea_mu_QFLX_def) + landTvslandPE_x(7,11)*(ea_mu_PRECT_30N-ea_mu_QFLX_30N);

% balanced
mu_TREFHT_bal_opt = landTvslandPE_x(1,6)*ea_mu_TREFHT_EQ + landTvslandPE_x(2,6)*ea_mu_TREFHT_15N_15S ...
+ landTvslandPE_x(3,6)*ea_mu_TREFHT_30N_30S + landTvslandPE_x(4,6)*ea_mu_TREFHT_60N_60S ...
+ landTvslandPE_x(5,6)*ea_mu_TREFHT_60N + landTvslandPE_x(6,6)*ea_mu_TREFHT_def + landTvslandPE_x(7,6)*ea_mu_TREFHT_30N;

mu_PE_bal_opt = landTvslandPE_x(1,6)*(ea_mu_PRECT_EQ-ea_mu_QFLX_EQ) + landTvslandPE_x(2,6)*(ea_mu_PRECT_15N_15S-ea_mu_QFLX_15N_15S) ...
+ landTvslandPE_x(3,6)*(ea_mu_PRECT_30N_30S-ea_mu_QFLX_30N_30S) + landTvslandPE_x(4,6)*(ea_mu_PRECT_60N_60S-ea_mu_QFLX_60N_60S) ...
+ landTvslandPE_x(5,6)*(ea_mu_PRECT_60N-ea_mu_QFLX_60N) + landTvslandPE_x(6,6)*(ea_mu_PRECT_def-ea_mu_QFLX_def) + landTvslandPE_x(7,6)*(ea_mu_PRECT_30N-ea_mu_QFLX_30N);

% land P-E optimal
mu_TREFHT_landPE_opt = landTvslandPE_x(1,1)*ea_mu_TREFHT_EQ + landTvslandPE_x(2,1)*ea_mu_TREFHT_15N_15S ...
+ landTvslandPE_x(3,1)*ea_mu_TREFHT_30N_30S + landTvslandPE_x(4,1)*ea_mu_TREFHT_60N_60S ...
+ landTvslandPE_x(5,1)*ea_mu_TREFHT_60N + landTvslandPE_x(6,1)*ea_mu_TREFHT_def + landTvslandPE_x(7,1)*ea_mu_TREFHT_30N;

mu_PE_landPE_opt = landTvslandPE_x(1,1)*(ea_mu_PRECT_EQ-ea_mu_QFLX_EQ) + landTvslandPE_x(2,1)*(ea_mu_PRECT_15N_15S-ea_mu_QFLX_15N_15S) ...
+ landTvslandPE_x(3,1)*(ea_mu_PRECT_30N_30S-ea_mu_QFLX_30N_30S) + landTvslandPE_x(4,1)*(ea_mu_PRECT_60N_60S-ea_mu_QFLX_60N_60S) ...
+ landTvslandPE_x(5,1)*(ea_mu_PRECT_60N-ea_mu_QFLX_60N) + landTvslandPE_x(6,1)*(ea_mu_PRECT_def-ea_mu_QFLX_def) + landTvslandPE_x(7,1)*(ea_mu_PRECT_30N-ea_mu_QFLX_30N);

Tmin = -3;
Tmax = 3;
Pmin = -2;
Pmax = 2;
[T_contours, T_color_map] = custom_color_map(Tmin,Tmax,25,'*RdBu');
[P_contours, P_color_map] = custom_color_map(Pmin,Pmax,11,'BrBG');

figure('Renderer', 'painters', 'Position', [20 60 800 800])
clf

% landT opt TREFHT
Z=abs((ea_alpha_TREFHT+mu_TREFHT_landT_opt)./(2*[T_se; T_se(1,:)]));
[lon_null,lat_null]=find(Z<1);
lon_null=(lon_null-1)/287*360;
lat_null=(lat_null-1)/191*180-90;

ax1=subplot(3,2,1);
ax1.Position = [.07 .7 .42 .25];
hold on
set(gca,'FontSize',11)
worldmap([-90 90],[-180 180])
setm(ax1, 'meridianlabel', 'off', 'parallellabel', 'off')
box on
contourfm(lat,lon2,ea_alpha_TREFHT'+mu_TREFHT_landT_opt',T_contours,'LineStyle','none')
scatterm(lat_null,lon_null,.2,'k','filled');
load coastlines
geoshow(coastlat,coastlon,'Color','k')
caxis([Tmin Tmax])
hold off
title('\DeltaT , Optimal land T strategy');
colormap(ax1, T_color_map);
add_abc_tofig(1,0.04,-0.04,12);

% T opt P-E
Z=abs((ea_alpha_PRECT-ea_alpha_QFLX+mu_PE_landT_opt)./(2*[PE_se; PE_se(1,:)]));
[lon_null,lat_null]=find(Z<1);
lon_null=(lon_null-1)/287*360;
lat_null=(lat_null-1)/191*180-90;

ax2=subplot(3,2,2);
ax2.Position = [.53 .7 .42 .25];
hold on
set(gca,'FontSize',11)
worldmap([-90 90],[-180 180])
setm(ax2, 'meridianlabel', 'off', 'parallellabel', 'off')
box on
contourfm(lat,lon2,ea_alpha_PRECT'-ea_alpha_QFLX'+mu_PE_landT_opt',P_contours,'LineStyle','none')
scatterm(lat_null,lon_null,.2,'k','filled');
load coastlines
geoshow(coastlat,coastlon,'Color','k')
caxis([Pmin Pmax])
hold off
title('\DeltaP-E , Optimal land T strategy');
colormap(ax2, P_color_map);
add_abc_tofig(2,0.04,-0.04,12);

% balanced TREFHT
Z=abs((ea_alpha_TREFHT+mu_TREFHT_bal_opt)./(2*[T_se; T_se(1,:)]));
[lon_null,lat_null]=find(Z<1);
lon_null=(lon_null-1)/287*360;
lat_null=(lat_null-1)/191*180-90;

ax3=subplot(3,2,3);
ax3.Position = [.07 .41 .42 .25];
hold on
set(gca,'FontSize',11)
worldmap([-90 90],[-180 180])
setm(ax3, 'meridianlabel', 'off', 'parallellabel', 'off')
box on
contourfm(lat,lon2,ea_alpha_TREFHT'+mu_TREFHT_bal_opt',T_contours,'LineStyle','none')
scatterm(lat_null,lon_null,.2,'k','filled');
load coastlines
geoshow(coastlat,coastlon,'Color','k')
caxis([Tmin Tmax])
hold off
title('\DeltaT , Balanced strategy');
colormap(ax3, T_color_map);
add_abc_tofig(3,0.04,-0.04,12);

% balanced P-E
Z=abs((ea_alpha_PRECT-ea_alpha_QFLX+mu_PE_bal_opt)./(2*[PE_se; PE_se(1,:)]));
[lon_null,lat_null]=find(Z<1);
lon_null=(lon_null-1)/287*360;
lat_null=(lat_null-1)/191*180-90;

ax4=subplot(3,2,4);
ax4.Position = [.53 .41 .42 .25];
hold on
set(gca,'FontSize',11)
worldmap([-90 90],[-180 180])
setm(ax4, 'meridianlabel', 'off', 'parallellabel', 'off')
box on
contourfm(lat,lon2,ea_alpha_PRECT'-ea_alpha_QFLX'+mu_PE_bal_opt',P_contours,'LineStyle','none')
scatterm(lat_null,lon_null,.2,'k','filled');
load coastlines
geoshow(coastlat,coastlon,'Color','k')
caxis([Pmin Pmax])
hold off
title('\DeltaP-E , Balanced strategy');
colormap(ax4, P_color_map);
add_abc_tofig(4,0.04,-0.04,12);

% land P-E opt TREFHT
Z=abs((ea_alpha_TREFHT+mu_TREFHT_landPE_opt)./(2*[T_se; T_se(1,:)]));
[lon_null,lat_null]=find(Z<1);
lon_null=(lon_null-1)/287*360;
lat_null=(lat_null-1)/191*180-90;

ax5=subplot(3,2,5);
ax5.Position = [.07 .02 .42 .35];
hold on
set(gca,'FontSize',11)
worldmap([-90 90],[-180 180])
setm(ax5, 'meridianlabel', 'off', 'parallellabel', 'off')
box on
contourfm(lat,lon2,ea_alpha_TREFHT'+mu_TREFHT_landPE_opt',T_contours,'LineStyle','none')
scatterm(lat_null,lon_null,.2,'k','filled');
load coastlines
geoshow(coastlat,coastlon,'Color','k')
caxis([Tmin Tmax])
cb=colorbar;
ylabel(cb,'K','FontSize',10,'Rotation',0,'FontWeight','bold')
cb.Location = 'southoutside';
hold off
title('\DeltaT , Optimal land P-E strategy');
colormap(ax5, T_color_map);
add_abc_tofig(5,0.04,-0.04,12);

% P opt PRECT
Z=abs((ea_alpha_PRECT-ea_alpha_QFLX+mu_PE_landPE_opt)./(2*[PE_se; PE_se(1,:)]));
[lon_null,lat_null]=find(Z<1);
lon_null=(lon_null-1)/287*360;
lat_null=(lat_null-1)/191*180-90;

ax6=subplot(3,2,6);
ax6.Position = [.53 .02 .42 .35];
hold on
set(gca,'FontSize',11)
worldmap([-90 90],[-180 180])
setm(ax6, 'meridianlabel', 'off', 'parallellabel', 'off')
box on
contourfm(lat,lon2,ea_alpha_PRECT'-ea_alpha_QFLX'+mu_PE_landPE_opt',P_contours,'LineStyle','none')
scatterm(lat_null,lon_null,.2,'k','filled');
load coastlines
geoshow(coastlat,coastlon,'Color','k')
caxis([Pmin Pmax])
cb=colorbar;
ylabel(cb,'mm/day','FontSize',12,'Rotation',0,'FontWeight','bold')
cb.Location = 'southoutside';
hold off
title('\DeltaP-E , Optimal land P-E strategy');
colormap(ax6, P_color_map);
add_abc_tofig(6,0.04,-0.04,12);

%% Perform optimization rms land TREFHT vs rms land PRECT-QFLX (3 ens) normalized per standard deviation

alpha_vec = [alpha_TREFHT_vec; alpha_PRECT_vec-alpha_QFLX_vec];

Delta_T0=1;
landT_landPE_norm_inj_each_lat_arr=zeros(7,11);
landTvslandPE_norm_x = zeros(7,11);
landT_landPE_norm_temp_residual_rms_arr=zeros(1,11);
landT_landPE_norm_PE_residual_rms_arr=zeros(1,11);
% temp_residual_rms_arr_land=zeros(1,11);
% prec_residual_rms_arr_land=zeros(1,11);

for i=0:0.1:1
    % weighting for temperature and precipitation
    W_t=i*ones(55296,1);
    W_p=(1-i)*ones(55296,1);

    beta_matrix=([mu_TREFHT_arr;mu_PRECT_arr-mu_QFLX_arr]+alpha_vec)./([reshape(T_std,[55296 1]); reshape(PE_std,[55296 1])]);
%     beta_matrix=beta_matrix*Delta_T0;
    beta_matrix_landT_landPE=beta_matrix.*sqrt(land_weight_vec.*[W_t;W_p]);
    % f=0.05*Delta_arr';
    f=zeros(7,1);
    % epsilon=0.001*ones(7,1);
    
    H=(beta_matrix_landT_landPE)'*(beta_matrix_landT_landPE);
    
    % fun=@(x) sqrt((beta_matrix*x)'*(beta_matrix*x));
    
    A=-C;
    b=zeros(7,1);
    Aeq=ones(1,7);
    beq=Delta_T0;
    % b=epsilon;
    % x0=[0.1 0.1 0.3 0 0.3 0.1 0.1];
    lb=[];
    ub=[];
    % [x,fval]=fmincon(fun,x0,A,b,lb,ub)
    [x,fval,exitflag,output]=quadprog(H,f,A,b,Aeq,beq,[],[]);
    landTvslandPE_norm_x(:,floor(10*i+1)) = x;
    % injection rate for 60N,30N,15N,0N,15S,30S,60S
    landT_landPE_norm_inj_each_lat=C*x;
    landT_landPE_norm_inj_each_lat_arr(:,floor(10*i+1))=landT_landPE_norm_inj_each_lat;

    % land only
    beta_matrix_land_opt=beta_matrix.*sqrt(land_weight_vec);
    landT_landPE_norm_temp_residual_rms=sqrt(x'*beta_matrix_land_opt(1:55296,:)'*beta_matrix_land_opt(1:55296,:)*x);
    landT_landPE_norm_temp_residual_rms_arr(floor(10*i+1))=landT_landPE_norm_temp_residual_rms;
    landT_landPE_norm_PE_residual_rms=sqrt(x'*beta_matrix_land_opt(55297:110592,:)'*beta_matrix_land_opt(55297:110592,:)*x);
    landT_landPE_norm_PE_residual_rms_arr(floor(10*i+1))=landT_landPE_norm_PE_residual_rms;
end

% save injection rates and r.m.s.
writetable(table(landT_landPE_norm_inj_each_lat_arr),'landT_landPE_injection_rates_different_weighting.csv');
writetable(table([landT_landPE_norm_temp_residual_rms_arr;landT_landPE_norm_PE_residual_rms_arr]),'landT_landPE_residual_different_weighting.csv')

figure()
hold on
scatter(landT_landPE_norm_temp_residual_rms_arr,landT_landPE_norm_PE_residual_rms_arr,'k')
plot(landT_landPE_norm_temp_residual_rms_arr,landT_landPE_norm_PE_residual_rms_arr,'k')

grid on
box on
xlabel('r.m.s. land temperature change [\sigma]')
ylabel('r.m.s. land P-E change [\sigma]')
saveas(gcf,'C:\Users\ezrab\Documents\Cornell\Research\Figures\SAI_optimization\Fig_landPE_rms_vs_landT_rms_norm_different_weighting_3ens.png')

%% maps of optimal landT and landP-E diff from ref period normalized per standard deviation
%create optimal mus
% land T optimal
mu_TREFHT_landT_norm_opt = landTvslandPE_norm_x(1,11)*ea_mu_TREFHT_EQ + landTvslandPE_norm_x(2,11)*ea_mu_TREFHT_15N_15S ...
+ landTvslandPE_norm_x(3,11)*ea_mu_TREFHT_30N_30S + landTvslandPE_norm_x(4,11)*ea_mu_TREFHT_60N_60S ...
+ landTvslandPE_norm_x(5,11)*ea_mu_TREFHT_60N + landTvslandPE_norm_x(6,11)*ea_mu_TREFHT_def + landTvslandPE_norm_x(7,11)*ea_mu_TREFHT_30N;

mu_PE_landT_norm_opt = landTvslandPE_norm_x(1,11)*(ea_mu_PRECT_EQ-ea_mu_QFLX_EQ) + landTvslandPE_norm_x(2,11)*(ea_mu_PRECT_15N_15S-ea_mu_QFLX_15N_15S) ...
+ landTvslandPE_norm_x(3,11)*(ea_mu_PRECT_30N_30S-ea_mu_QFLX_30N_30S) + landTvslandPE_norm_x(4,11)*(ea_mu_PRECT_60N_60S-ea_mu_QFLX_60N_60S) ...
+ landTvslandPE_norm_x(5,11)*(ea_mu_PRECT_60N-ea_mu_QFLX_60N) + landTvslandPE_norm_x(6,11)*(ea_mu_PRECT_def-ea_mu_QFLX_def) + landTvslandPE_norm_x(7,11)*(ea_mu_PRECT_30N-ea_mu_QFLX_30N);

% balanced
mu_TREFHT_bal_norm_opt = landTvslandPE_norm_x(1,6)*ea_mu_TREFHT_EQ + landTvslandPE_norm_x(2,6)*ea_mu_TREFHT_15N_15S ...
+ landTvslandPE_norm_x(3,6)*ea_mu_TREFHT_30N_30S + landTvslandPE_norm_x(4,6)*ea_mu_TREFHT_60N_60S ...
+ landTvslandPE_norm_x(5,6)*ea_mu_TREFHT_60N + landTvslandPE_norm_x(6,6)*ea_mu_TREFHT_def + landTvslandPE_norm_x(7,6)*ea_mu_TREFHT_30N;

mu_PE_bal_norm_opt = landTvslandPE_norm_x(1,6)*(ea_mu_PRECT_EQ-ea_mu_QFLX_EQ) + landTvslandPE_norm_x(2,6)*(ea_mu_PRECT_15N_15S-ea_mu_QFLX_15N_15S) ...
+ landTvslandPE_norm_x(3,6)*(ea_mu_PRECT_30N_30S-ea_mu_QFLX_30N_30S) + landTvslandPE_norm_x(4,6)*(ea_mu_PRECT_60N_60S-ea_mu_QFLX_60N_60S) ...
+ landTvslandPE_norm_x(5,6)*(ea_mu_PRECT_60N-ea_mu_QFLX_60N) + landTvslandPE_norm_x(6,6)*(ea_mu_PRECT_def-ea_mu_QFLX_def) + landTvslandPE_norm_x(7,6)*(ea_mu_PRECT_30N-ea_mu_QFLX_30N);

% land P-E optimal
mu_TREFHT_landPE_norm_opt = landTvslandPE_norm_x(1,1)*ea_mu_TREFHT_EQ + landTvslandPE_norm_x(2,1)*ea_mu_TREFHT_15N_15S ...
+ landTvslandPE_norm_x(3,1)*ea_mu_TREFHT_30N_30S + landTvslandPE_norm_x(4,1)*ea_mu_TREFHT_60N_60S ...
+ landTvslandPE_norm_x(5,1)*ea_mu_TREFHT_60N + landTvslandPE_norm_x(6,1)*ea_mu_TREFHT_def + landTvslandPE_norm_x(7,1)*ea_mu_TREFHT_30N;

mu_PE_landPE_norm_opt = landTvslandPE_norm_x(1,1)*(ea_mu_PRECT_EQ-ea_mu_QFLX_EQ) + landTvslandPE_norm_x(2,1)*(ea_mu_PRECT_15N_15S-ea_mu_QFLX_15N_15S) ...
+ landTvslandPE_norm_x(3,1)*(ea_mu_PRECT_30N_30S-ea_mu_QFLX_30N_30S) + landTvslandPE_norm_x(4,1)*(ea_mu_PRECT_60N_60S-ea_mu_QFLX_60N_60S) ...
+ landTvslandPE_norm_x(5,1)*(ea_mu_PRECT_60N-ea_mu_QFLX_60N) + landTvslandPE_norm_x(6,1)*(ea_mu_PRECT_def-ea_mu_QFLX_def) + landTvslandPE_norm_x(7,1)*(ea_mu_PRECT_30N-ea_mu_QFLX_30N);

Tmin = -4;
Tmax = 4;
Pmin = -3;
Pmax = 3;
[T_contours, T_color_map] = custom_color_map(Tmin,Tmax,25,'*RdBu');
[P_contours, P_color_map] = custom_color_map(Pmin,Pmax,11,'BrBG');

figure('Renderer', 'painters', 'Position', [20 60 800 800])
clf

% landT opt TREFHT
Z=abs((ea_alpha_TREFHT+mu_TREFHT_landT_norm_opt)./(2*[T_se; T_se(1,:)]));
[lon_null,lat_null]=find(Z<1);
lon_null=(lon_null-1)/287*360;
lat_null=(lat_null-1)/191*180-90;

ax1=subplot(3,2,1);
ax1.Position = [.07 .7 .42 .25];
hold on
set(gca,'FontSize',11)
worldmap([-90 90],[-180 180])
setm(ax1, 'meridianlabel', 'off', 'parallellabel', 'off')
box on
contourfm(lat,lon2,((ea_alpha_TREFHT+mu_TREFHT_landT_norm_opt)./[T_std; T_std(1,:)])',T_contours,'LineStyle','none')
scatterm(lat_null,lon_null,.2,'k','filled');
load coastlines
geoshow(coastlat,coastlon,'Color','k')
caxis([Tmin Tmax])
hold off
title('\DeltaT , Optimal land T strategy');
colormap(ax1, T_color_map);
add_abc_tofig(1,0.04,-0.04,12);

% T opt P-E
Z=abs((ea_alpha_PRECT-ea_alpha_QFLX+mu_PE_landT_norm_opt)./(2*[PE_se; PE_se(1,:)]));
[lon_null,lat_null]=find(Z<1);
lon_null=(lon_null-1)/287*360;
lat_null=(lat_null-1)/191*180-90;

ax2=subplot(3,2,2);
ax2.Position = [.53 .7 .42 .25];
hold on
set(gca,'FontSize',11)
worldmap([-90 90],[-180 180])
setm(ax2, 'meridianlabel', 'off', 'parallellabel', 'off')
box on
contourfm(lat,lon2,((ea_alpha_PRECT-ea_alpha_QFLX+mu_PE_landT_norm_opt)./[PE_std; PE_std(1,:)])',P_contours,'LineStyle','none')
scatterm(lat_null,lon_null,.2,'k','filled');
load coastlines
geoshow(coastlat,coastlon,'Color','k')
caxis([Pmin Pmax])
hold off
title('\DeltaP-E , Optimal land T strategy');
colormap(ax2, P_color_map);
add_abc_tofig(2,0.04,-0.04,12);

% balanced TREFHT
Z=abs((ea_alpha_TREFHT+mu_TREFHT_bal_norm_opt)./(2*[T_se; T_se(1,:)]));
[lon_null,lat_null]=find(Z<1);
lon_null=(lon_null-1)/287*360;
lat_null=(lat_null-1)/191*180-90;

ax3=subplot(3,2,3);
ax3.Position = [.07 .41 .42 .25];
hold on
set(gca,'FontSize',11)
worldmap([-90 90],[-180 180])
setm(ax3, 'meridianlabel', 'off', 'parallellabel', 'off')
box on
contourfm(lat,lon2,((ea_alpha_TREFHT+mu_TREFHT_bal_norm_opt)./[T_std; T_std(1,:)])',T_contours,'LineStyle','none')
scatterm(lat_null,lon_null,.2,'k','filled');
load coastlines
geoshow(coastlat,coastlon,'Color','k')
caxis([Tmin Tmax])
hold off
title('\DeltaT , Balanced strategy');
colormap(ax3, T_color_map);
add_abc_tofig(3,0.04,-0.04,12);

% balanced P-E
Z=abs((ea_alpha_PRECT-ea_alpha_QFLX+mu_PE_bal_norm_opt)./(2*[PE_se; PE_se(1,:)]));
[lon_null,lat_null]=find(Z<1);
lon_null=(lon_null-1)/287*360;
lat_null=(lat_null-1)/191*180-90;

ax4=subplot(3,2,4);
ax4.Position = [.53 .41 .42 .25];
hold on
set(gca,'FontSize',11)
worldmap([-90 90],[-180 180])
setm(ax4, 'meridianlabel', 'off', 'parallellabel', 'off')
box on
contourfm(lat,lon2,((ea_alpha_PRECT-ea_alpha_QFLX+mu_PE_bal_norm_opt)./[PE_std; PE_std(1,:)])',P_contours,'LineStyle','none')
scatterm(lat_null,lon_null,.2,'k','filled');
load coastlines
geoshow(coastlat,coastlon,'Color','k')
caxis([Pmin Pmax])
hold off
title('\DeltaP-E , Balanced strategy');
colormap(ax4, P_color_map);
add_abc_tofig(4,0.04,-0.04,12);

% land P-E opt TREFHT
Z=abs((ea_alpha_TREFHT+mu_TREFHT_landPE_norm_opt)./(2*[T_se; T_se(1,:)]));
[lon_null,lat_null]=find(Z<1);
lon_null=(lon_null-1)/287*360;
lat_null=(lat_null-1)/191*180-90;

ax5=subplot(3,2,5);
ax5.Position = [.07 .02 .42 .35];
hold on
set(gca,'FontSize',11)
worldmap([-90 90],[-180 180])
setm(ax5, 'meridianlabel', 'off', 'parallellabel', 'off')
box on
contourfm(lat,lon2,((ea_alpha_TREFHT+mu_TREFHT_landPE_norm_opt)./[T_std; T_std(1,:)])',T_contours,'LineStyle','none')
scatterm(lat_null,lon_null,.2,'k','filled');
load coastlines
geoshow(coastlat,coastlon,'Color','k')
caxis([Tmin Tmax])
cb=colorbar;
ylabel(cb,'\sigma','FontSize',10,'Rotation',0,'FontWeight','bold')
cb.Location = 'southoutside';
hold off
title('\DeltaT , Optimal land P-E strategy');
colormap(ax5, T_color_map);
add_abc_tofig(5,0.04,-0.04,12);

% P opt PRECT
Z=abs((ea_alpha_PRECT-ea_alpha_QFLX+mu_PE_landPE_norm_opt)./(2*[PE_se; PE_se(1,:)]));
[lon_null,lat_null]=find(Z<1);
lon_null=(lon_null-1)/287*360;
lat_null=(lat_null-1)/191*180-90;

ax6=subplot(3,2,6);
ax6.Position = [.53 .02 .42 .35];
hold on
set(gca,'FontSize',11)
worldmap([-90 90],[-180 180])
setm(ax6, 'meridianlabel', 'off', 'parallellabel', 'off')
box on
contourfm(lat,lon2,((ea_alpha_PRECT-ea_alpha_QFLX+mu_PE_landPE_norm_opt)./[PE_std; PE_std(1,:)])',P_contours,'LineStyle','none')
scatterm(lat_null,lon_null,.2,'k','filled');
load coastlines
geoshow(coastlat,coastlon,'Color','k')
caxis([Pmin Pmax])
cb=colorbar;
ylabel(cb,'\sigma','FontSize',12,'Rotation',0,'FontWeight','bold')
cb.Location = 'southoutside';
hold off
title('\DeltaP-E , Optimal land P-E strategy');
colormap(ax6, P_color_map);
add_abc_tofig(6,0.04,-0.04,12);

%% maps of optimal landT and landP-E diff from ref period normalized per standard deviation, shown in absolute units
%create optimal mus
% land T optimal
mu_TREFHT_landT_norm_opt = landTvslandPE_norm_x(1,11)*ea_mu_TREFHT_EQ + landTvslandPE_norm_x(2,11)*ea_mu_TREFHT_15N_15S ...
+ landTvslandPE_norm_x(3,11)*ea_mu_TREFHT_30N_30S + landTvslandPE_norm_x(4,11)*ea_mu_TREFHT_60N_60S ...
+ landTvslandPE_norm_x(5,11)*ea_mu_TREFHT_60N + landTvslandPE_norm_x(6,11)*ea_mu_TREFHT_def + landTvslandPE_norm_x(7,11)*ea_mu_TREFHT_30N;

mu_PE_landT_norm_opt = landTvslandPE_norm_x(1,11)*(ea_mu_PRECT_EQ-ea_mu_QFLX_EQ) + landTvslandPE_norm_x(2,11)*(ea_mu_PRECT_15N_15S-ea_mu_QFLX_15N_15S) ...
+ landTvslandPE_norm_x(3,11)*(ea_mu_PRECT_30N_30S-ea_mu_QFLX_30N_30S) + landTvslandPE_norm_x(4,11)*(ea_mu_PRECT_60N_60S-ea_mu_QFLX_60N_60S) ...
+ landTvslandPE_norm_x(5,11)*(ea_mu_PRECT_60N-ea_mu_QFLX_60N) + landTvslandPE_norm_x(6,11)*(ea_mu_PRECT_def-ea_mu_QFLX_def) + landTvslandPE_norm_x(7,11)*(ea_mu_PRECT_30N-ea_mu_QFLX_30N);

% balanced
mu_TREFHT_bal_norm_opt = landTvslandPE_norm_x(1,6)*ea_mu_TREFHT_EQ + landTvslandPE_norm_x(2,6)*ea_mu_TREFHT_15N_15S ...
+ landTvslandPE_norm_x(3,6)*ea_mu_TREFHT_30N_30S + landTvslandPE_norm_x(4,6)*ea_mu_TREFHT_60N_60S ...
+ landTvslandPE_norm_x(5,6)*ea_mu_TREFHT_60N + landTvslandPE_norm_x(6,6)*ea_mu_TREFHT_def + landTvslandPE_norm_x(7,6)*ea_mu_TREFHT_30N;

mu_PE_bal_norm_opt = landTvslandPE_norm_x(1,6)*(ea_mu_PRECT_EQ-ea_mu_QFLX_EQ) + landTvslandPE_norm_x(2,6)*(ea_mu_PRECT_15N_15S-ea_mu_QFLX_15N_15S) ...
+ landTvslandPE_norm_x(3,6)*(ea_mu_PRECT_30N_30S-ea_mu_QFLX_30N_30S) + landTvslandPE_norm_x(4,6)*(ea_mu_PRECT_60N_60S-ea_mu_QFLX_60N_60S) ...
+ landTvslandPE_norm_x(5,6)*(ea_mu_PRECT_60N-ea_mu_QFLX_60N) + landTvslandPE_norm_x(6,6)*(ea_mu_PRECT_def-ea_mu_QFLX_def) + landTvslandPE_norm_x(7,6)*(ea_mu_PRECT_30N-ea_mu_QFLX_30N);

% land P-E optimal
mu_TREFHT_landPE_norm_opt = landTvslandPE_norm_x(1,1)*ea_mu_TREFHT_EQ + landTvslandPE_norm_x(2,1)*ea_mu_TREFHT_15N_15S ...
+ landTvslandPE_norm_x(3,1)*ea_mu_TREFHT_30N_30S + landTvslandPE_norm_x(4,1)*ea_mu_TREFHT_60N_60S ...
+ landTvslandPE_norm_x(5,1)*ea_mu_TREFHT_60N + landTvslandPE_norm_x(6,1)*ea_mu_TREFHT_def + landTvslandPE_norm_x(7,1)*ea_mu_TREFHT_30N;

mu_PE_landPE_norm_opt = landTvslandPE_norm_x(1,1)*(ea_mu_PRECT_EQ-ea_mu_QFLX_EQ) + landTvslandPE_norm_x(2,1)*(ea_mu_PRECT_15N_15S-ea_mu_QFLX_15N_15S) ...
+ landTvslandPE_norm_x(3,1)*(ea_mu_PRECT_30N_30S-ea_mu_QFLX_30N_30S) + landTvslandPE_norm_x(4,1)*(ea_mu_PRECT_60N_60S-ea_mu_QFLX_60N_60S) ...
+ landTvslandPE_norm_x(5,1)*(ea_mu_PRECT_60N-ea_mu_QFLX_60N) + landTvslandPE_norm_x(6,1)*(ea_mu_PRECT_def-ea_mu_QFLX_def) + landTvslandPE_norm_x(7,1)*(ea_mu_PRECT_30N-ea_mu_QFLX_30N);

Tmin = -3;
Tmax = 3;
Pmin = -2;
Pmax = 2;
[T_contours, T_color_map] = custom_color_map(Tmin,Tmax,25,'*RdBu');
[P_contours, P_color_map] = custom_color_map(Pmin,Pmax,11,'BrBG');

figure('Renderer', 'painters', 'Position', [20 60 800 800])
clf

% landT opt TREFHT
Z=abs((ea_alpha_TREFHT+mu_TREFHT_landT_norm_opt)./(2*[T_se; T_se(1,:)]));
[lon_null,lat_null]=find(Z<1);
lon_null=(lon_null-1)/287*360;
lat_null=(lat_null-1)/191*180-90;

ax1=subplot(3,2,1);
ax1.Position = [.07 .7 .42 .25];
hold on
set(gca,'FontSize',11)
worldmap([-90 90],[-180 180])
setm(ax1, 'meridianlabel', 'off', 'parallellabel', 'off')
box on
contourfm(lat,lon2,(ea_alpha_TREFHT+mu_TREFHT_landT_norm_opt)',T_contours,'LineStyle','none')
scatterm(lat_null,lon_null,.2,'k','filled');
load coastlines
geoshow(coastlat,coastlon,'Color','k')
caxis([Tmin Tmax])
hold off
title('\DeltaT , Optimal land T strategy');
colormap(ax1, T_color_map);
add_abc_tofig(1,0.04,-0.04,12);

% T opt P-E
Z=abs((ea_alpha_PRECT-ea_alpha_QFLX+mu_PE_landT_norm_opt)./(2*[PE_se; PE_se(1,:)]));
[lon_null,lat_null]=find(Z<1);
lon_null=(lon_null-1)/287*360;
lat_null=(lat_null-1)/191*180-90;

ax2=subplot(3,2,2);
ax2.Position = [.53 .7 .42 .25];
hold on
set(gca,'FontSize',11)
worldmap([-90 90],[-180 180])
setm(ax2, 'meridianlabel', 'off', 'parallellabel', 'off')
box on
contourfm(lat,lon2,(ea_alpha_PRECT-ea_alpha_QFLX+mu_PE_landT_norm_opt)',P_contours,'LineStyle','none')
scatterm(lat_null,lon_null,.2,'k','filled');
load coastlines
geoshow(coastlat,coastlon,'Color','k')
caxis([Pmin Pmax])
hold off
title('\DeltaP-E , Optimal land T strategy');
colormap(ax2, P_color_map);
add_abc_tofig(2,0.04,-0.04,12);

% balanced TREFHT
Z=abs((ea_alpha_TREFHT+mu_TREFHT_bal_norm_opt)./(2*[T_se; T_se(1,:)]));
[lon_null,lat_null]=find(Z<1);
lon_null=(lon_null-1)/287*360;
lat_null=(lat_null-1)/191*180-90;

ax3=subplot(3,2,3);
ax3.Position = [.07 .41 .42 .25];
hold on
set(gca,'FontSize',11)
worldmap([-90 90],[-180 180])
setm(ax3, 'meridianlabel', 'off', 'parallellabel', 'off')
box on
contourfm(lat,lon2,(ea_alpha_TREFHT+mu_TREFHT_bal_norm_opt)',T_contours,'LineStyle','none')
scatterm(lat_null,lon_null,.2,'k','filled');
load coastlines
geoshow(coastlat,coastlon,'Color','k')
caxis([Tmin Tmax])
hold off
title('\DeltaT , Balanced strategy');
colormap(ax3, T_color_map);
add_abc_tofig(3,0.04,-0.04,12);

% balanced P-E
Z=abs((ea_alpha_PRECT-ea_alpha_QFLX+mu_PE_bal_norm_opt)./(2*[PE_se; PE_se(1,:)]));
[lon_null,lat_null]=find(Z<1);
lon_null=(lon_null-1)/287*360;
lat_null=(lat_null-1)/191*180-90;

ax4=subplot(3,2,4);
ax4.Position = [.53 .41 .42 .25];
hold on
set(gca,'FontSize',11)
worldmap([-90 90],[-180 180])
setm(ax4, 'meridianlabel', 'off', 'parallellabel', 'off')
box on
contourfm(lat,lon2,(ea_alpha_PRECT-ea_alpha_QFLX+mu_PE_bal_norm_opt)',P_contours,'LineStyle','none')
scatterm(lat_null,lon_null,.2,'k','filled');
load coastlines
geoshow(coastlat,coastlon,'Color','k')
caxis([Pmin Pmax])
hold off
title('\DeltaP-E , Balanced strategy');
colormap(ax4, P_color_map);
add_abc_tofig(4,0.04,-0.04,12);

% land P-E opt TREFHT
Z=abs((ea_alpha_TREFHT+mu_TREFHT_landPE_norm_opt)./(2*[T_se; T_se(1,:)]));
[lon_null,lat_null]=find(Z<1);
lon_null=(lon_null-1)/287*360;
lat_null=(lat_null-1)/191*180-90;

ax5=subplot(3,2,5);
ax5.Position = [.07 .02 .42 .35];
hold on
set(gca,'FontSize',11)
worldmap([-90 90],[-180 180])
setm(ax5, 'meridianlabel', 'off', 'parallellabel', 'off')
box on
contourfm(lat,lon2,(ea_alpha_TREFHT+mu_TREFHT_landPE_norm_opt)',T_contours,'LineStyle','none')
scatterm(lat_null,lon_null,.2,'k','filled');
load coastlines
geoshow(coastlat,coastlon,'Color','k')
caxis([Tmin Tmax])
cb=colorbar;
ylabel(cb,'^\circC','FontSize',10,'Rotation',0,'FontWeight','bold')
cb.Location = 'southoutside';
hold off
title('\DeltaT , Optimal land P-E strategy');
colormap(ax5, T_color_map);
add_abc_tofig(5,0.04,-0.04,12);

% P opt PRECT
Z=abs((ea_alpha_PRECT-ea_alpha_QFLX+mu_PE_landPE_norm_opt)./(2*[PE_se; PE_se(1,:)]));
[lon_null,lat_null]=find(Z<1);
lon_null=(lon_null-1)/287*360;
lat_null=(lat_null-1)/191*180-90;

ax6=subplot(3,2,6);
ax6.Position = [.53 .02 .42 .35];
hold on
set(gca,'FontSize',11)
worldmap([-90 90],[-180 180])
setm(ax6, 'meridianlabel', 'off', 'parallellabel', 'off')
box on
contourfm(lat,lon2,(ea_alpha_PRECT-ea_alpha_QFLX+mu_PE_landPE_norm_opt)',P_contours,'LineStyle','none')
scatterm(lat_null,lon_null,.2,'k','filled');
load coastlines
geoshow(coastlat,coastlon,'Color','k')
caxis([Pmin Pmax])
cb=colorbar;
ylabel(cb,'mm/day','FontSize',12,'Rotation',0,'FontWeight','bold')
cb.Location = 'southoutside';
hold off
title('\DeltaP-E , Optimal land P-E strategy');
colormap(ax6, P_color_map);
add_abc_tofig(6,0.04,-0.04,12);

%% COLORS FOR RMS PLOTS

colors=[0.5,0.5,0.5;
    0.7176 0.2745 1.0000;
    0.8500 0.3250 0.0980;
    0.9290, 0.6940, 0.1250;
    0.4660 0.6740 0.1880;
    0 0.4470 0.7410;
    0 1 0;
    0 1 1];
%% plot RMS T vs P (all members) with pareto
figure('Renderer', 'painters', 'Position', [100 200 420 460])
clf
ax = axes;
hold on
plot(TvsP_temp_residual_rms_arr, TvsP_prec_residual_rms_arr,'k','HandleVisibility','off')
scatter(TvsP_temp_residual_rms_arr, TvsP_prec_residual_rms_arr,'k','filled','marker','square','HandleVisibility','off')
scatter(ea_rms_T_def, ea_rms_P_def,"filled","MarkerFaceColor",colors(2,:))
scatter(ea_rms_T_EQ, ea_rms_P_EQ,"filled","MarkerFaceColor",colors(3,:))
scatter(ea_rms_T_15N_15S, ea_rms_P_15N_15S,"filled","MarkerFaceColor",colors(4,:))
scatter(ea_rms_T_30N_30S, ea_rms_P_30N_30S,"filled","MarkerFaceColor",colors(5,:))
scatter(ea_rms_T_60N_60S, ea_rms_P_60N_60S,"filled","MarkerFaceColor",colors(6,:))
%scatter(rms_T_predicted, rms_P_predicted,"filled","MarkerFaceColor",colors(2,:),'Marker','square')
hold off
grid on
xlabel('r.m.s. temperature change [^\circ C]')
ylabel('r.m.s. precipitation change [mm/day]')
legend('Multi-Objective','EQ','15N+15S', '30N+30S', '60N+60S',['Multi-Objective' newline 'predicted'],'Location','southoutside')

%% plot RMS T vs P (all members) with predicted
figure('Renderer', 'painters', 'Position', [100 200 420 500])
clf
ax = axes;
hold on
plot(TvsP_temp_residual_rms_arr, TvsP_prec_residual_rms_arr,'k','HandleVisibility','off')
scatter(TvsP_temp_residual_rms_arr, TvsP_prec_residual_rms_arr,'k','filled','marker','square','HandleVisibility','off')
scatter(ea_rms_T_def, ea_rms_P_def,"filled","MarkerFaceColor",colors(2,:))
scatter(ea_rms_T_EQ, ea_rms_P_EQ,"filled","MarkerFaceColor",colors(3,:))
scatter(ea_rms_T_15N_15S, ea_rms_P_15N_15S,"filled","MarkerFaceColor",colors(4,:))
scatter(ea_rms_T_30N_30S, ea_rms_P_30N_30S,"filled","MarkerFaceColor",colors(5,:))
scatter(ea_rms_T_60N_60S, ea_rms_P_60N_60S,"filled","MarkerFaceColor",colors(6,:))
scatter(rms_T_predicted, rms_P_predicted,"filled","MarkerFaceColor",colors(2,:),'Marker','square')
hold off
grid on
xlabel('r.m.s. temperature change [^\circ C]')
ylabel('r.m.s. precipitation change [mm/day]')
legend('Multi-Objective','EQ','15N+15S', '30N+30S', '60N+60S',['Multi-Objective' newline 'predicted'],'Location','southoutside')

%% plot RMS T vs P normalized per std error (all members) pareto predicted
figure('Renderer', 'painters', 'Position', [100 200 420 500])
clf
ax = axes;
hold on
plot(TvsP_norm_temp_residual_rms_arr, TvsP_norm_prec_residual_rms_arr,'k','HandleVisibility','off')
scatter(TvsP_norm_temp_residual_rms_arr, TvsP_norm_prec_residual_rms_arr,'k','filled','marker','square','HandleVisibility','off')
scatter(ea_rms_T_def_norm, ea_rms_P_def_norm,"filled","SizeData",200,"MarkerFaceColor",colors(2,:),'MarkerEdgeColor','k')
scatter(ea_rms_T_EQ_norm, ea_rms_P_EQ_norm,"filled","SizeData",200,"MarkerFaceColor",colors(3,:),'MarkerEdgeColor','k')
scatter(ea_rms_T_15N_15S_norm, ea_rms_P_15N_15S_norm,"filled","SizeData",200,"SizeData",200,"MarkerFaceColor",colors(4,:),'MarkerEdgeColor','k')
scatter(ea_rms_T_30N_30S_norm, ea_rms_P_30N_30S_norm,"filled","SizeData",200,"MarkerFaceColor",colors(5,:),'MarkerEdgeColor','k')
scatter(ea_rms_T_60N_60S_norm, ea_rms_P_60N_60S_norm,"filled","SizeData",200,"MarkerFaceColor",colors(6,:),'MarkerEdgeColor','k')
scatter(rms_T_norm_predicted, rms_P_norm_predicted,"filled","MarkerFaceColor",colors(2,:),'Marker','square')
for i = 1:3
    scatter(rms_T_def_norm(i), rms_P_def_norm(i),"filled","MarkerFaceColor",colors(2,:),'MarkerEdgeColor','k')
    scatter(rms_T_EQ_norm(i), rms_P_EQ_norm(i),"filled","MarkerFaceColor",colors(3,:),'MarkerEdgeColor','k')
    scatter(rms_T_15N_15S_norm(i), rms_P_15N_15S_norm(i),"filled","MarkerFaceColor",colors(4,:),'MarkerEdgeColor','k')
    scatter(rms_T_30N_30S_norm(i), rms_P_30N_30S_norm(i),"filled","MarkerFaceColor",colors(5,:),'MarkerEdgeColor','k')
    scatter(rms_T_60N_60S_norm(i), rms_P_60N_60S_norm(i),"filled","MarkerFaceColor",colors(6,:),'MarkerEdgeColor','k')
end
hold off
grid on
xlabel('r.m.s. temperature change [\sigma]')
ylabel('r.m.s. precipitation change [\sigma]')
legend('Multi-Objective','EQ','15N+15S', '30N+30S', '60N+60S',['Multi-Objective' newline 'predicted'],'Location','southoutside')

%% plot RMS T vs P normalized per standard deviation (all members) pareto
figure('Renderer', 'painters', 'Position', [100 200 420 460])
clf
ax = axes;
hold on
plot(TvsP_norm_temp_residual_rms_arr, TvsP_norm_prec_residual_rms_arr,'k','HandleVisibility','off')
scatter(TvsP_norm_temp_residual_rms_arr, TvsP_norm_prec_residual_rms_arr,'k','filled','marker','square','HandleVisibility','off')
scatter(ea_rms_T_def_norm, ea_rms_P_def_norm,"filled","MarkerFaceColor",colors(2,:))
scatter(ea_rms_T_EQ_norm, ea_rms_P_EQ_norm,"filled","MarkerFaceColor",colors(3,:))
scatter(ea_rms_T_15N_15S_norm, ea_rms_P_15N_15S_norm,"filled","MarkerFaceColor",colors(4,:))
scatter(ea_rms_T_30N_30S_norm, ea_rms_P_30N_30S_norm,"filled","MarkerFaceColor",colors(5,:))
scatter(ea_rms_T_60N_60S_norm, ea_rms_P_60N_60S_norm,"filled","MarkerFaceColor",colors(6,:))
hold off
grid on
xlabel('r.m.s. temperature change [\sigma]')
ylabel('r.m.s. precipitation change [\sigma]')
legend('Multi-Objective','EQ','15N+15S', '30N+30S', '60N+60S',['Multi-Objective' newline 'predicted'],'Location','southoutside')

%% plot RMS T vs P (all members) with predicted, all strategies + baseline
figure('Renderer', 'painters', 'Position', [100 200 560 320])
clf
ax = axes;
hold on
plot(TvsP_temp_residual_rms_arr, TvsP_prec_residual_rms_arr,'k','HandleVisibility','off')
scatter(TvsP_temp_residual_rms_arr, TvsP_prec_residual_rms_arr,'k','filled','marker','square','HandleVisibility','off')
scatter(ea_rms_T_base, ea_rms_P_base,"filled","MarkerFaceColor",colors(1,:))
scatter(ea_rms_T_def, ea_rms_P_def,"filled","MarkerFaceColor",colors(2,:))
scatter(ea_rms_T_EQ, ea_rms_P_EQ,"filled","MarkerFaceColor",colors(3,:))
scatter(ea_rms_T_15N_15S, ea_rms_P_15N_15S,"filled","MarkerFaceColor",colors(4,:))
scatter(ea_rms_T_30N_30S, ea_rms_P_30N_30S,"filled","MarkerFaceColor",colors(5,:))
scatter(ea_rms_T_60N_60S, ea_rms_P_60N_60S,"filled","MarkerFaceColor",colors(6,:))
scatter(ea_rms_T_30N, ea_rms_P_30N,"filled","MarkerFaceColor",colors(7,:))
scatter(ea_rms_T_60N, ea_rms_P_60N,"filled","MarkerFaceColor",colors(8,:))
scatter(rms_T_predicted, rms_P_predicted,"filled","MarkerFaceColor",colors(2,:),'Marker','square')
hold off
grid on
xlabel('r.m.s. temperature change [^\circ C]')
ylabel('r.m.s. precipitation change [mm/day]')
legend('SSP2-4.5','Multi-Objective','EQ','15N+15S', '30N+30S', '60N+60S','30N','60N',['Multi-Objective' newline 'predicted'],'Location','eastoutside')

%% plot RMS T vs P for the strategies, individual members w pareto front
figure('Renderer', 'painters', 'Position', [100 200 450 550])
clf
ax = axes;
hold on
plot(TvsP_temp_residual_rms_arr, TvsP_prec_residual_rms_arr,'k','HandleVisibility','off')
scatter(TvsP_temp_residual_rms_arr, TvsP_prec_residual_rms_arr,'k','filled','marker','square','HandleVisibility','off')
scatter(ea_rms_T_def, ea_rms_P_def,"filled","MarkerFaceColor",colors(2,:),"SizeData",200,'MarkerEdgeColor','k')
scatter(ea_rms_T_EQ, ea_rms_P_EQ,"filled","MarkerFaceColor",colors(3,:),"SizeData",200,'MarkerEdgeColor','k')
scatter(ea_rms_T_15N_15S, ea_rms_P_15N_15S,"filled","MarkerFaceColor",colors(4,:),"SizeData",200,'MarkerEdgeColor','k')
scatter(ea_rms_T_30N_30S, ea_rms_P_30N_30S,"filled","MarkerFaceColor",colors(5,:),"SizeData",200,'MarkerEdgeColor','k')
scatter(ea_rms_T_60N_60S, ea_rms_P_60N_60S,"filled","MarkerFaceColor",colors(6,:),"SizeData",200,'MarkerEdgeColor','k')
scatter(rms_T_predicted, rms_P_predicted,"filled","MarkerFaceColor",colors(2,:),'MarkerEdgeColor','k','Marker','square')
for i = 1:3
    scatter(rms_T_def(i), rms_P_def(i),"filled","MarkerFaceColor",colors(2,:),'MarkerEdgeColor','k')
    scatter(rms_T_EQ(i), rms_P_EQ(i),"filled","MarkerFaceColor",colors(3,:),'MarkerEdgeColor','k')
    scatter(rms_T_15N_15S(i), rms_P_15N_15S(i),"filled","MarkerFaceColor",colors(4,:),'MarkerEdgeColor','k')
    scatter(rms_T_30N_30S(i), rms_P_30N_30S(i),"filled","MarkerFaceColor",colors(5,:),'MarkerEdgeColor','k')
    scatter(rms_T_60N_60S(i), rms_P_60N_60S(i),"filled","MarkerFaceColor",colors(6,:),'MarkerEdgeColor','k')
end
hold off
grid on
xlabel('r.m.s. temperature change [^\circ C]')
ylabel('r.m.s. precipitation change [mm/day]')
legend('Multi-Objective','EQ','15N+15S', '30N+30S', '60N+60S',['Multi-Objective' newline 'predicted'],'Location','southoutside')

%% plot RMS T vs P (all members) with pareto /AND/ RMS T vs P for the strategies, individual members w pareto front
figure('Renderer', 'painters', 'Position', [100 200 1120 400])
clf
ax1 = subplot(1,2,1);
ax1.Position = [0.07 0.15 0.3 0.8];
hold on
plot(TvsP_temp_residual_rms_arr, TvsP_prec_residual_rms_arr,'k','HandleVisibility','off')
scatter(TvsP_temp_residual_rms_arr, TvsP_prec_residual_rms_arr,'k','filled','marker','square','HandleVisibility','off')
scatter(ea_rms_T_base, ea_rms_P_base,'k',"filled")
scatter(ea_rms_T_def, ea_rms_P_def,"filled","MarkerFaceColor",colors(2,:))
scatter(ea_rms_T_EQ, ea_rms_P_EQ,"filled","MarkerFaceColor",colors(3,:))
scatter(ea_rms_T_15N_15S, ea_rms_P_15N_15S,"filled","MarkerFaceColor",colors(4,:))
scatter(ea_rms_T_30N_30S, ea_rms_P_30N_30S,"filled","MarkerFaceColor",colors(5,:))
scatter(ea_rms_T_60N_60S, ea_rms_P_60N_60S,"filled","MarkerFaceColor",colors(6,:))
scatter(ea_rms_T_30N, ea_rms_P_30N,"filled","MarkerFaceColor",colors(7,:))
scatter(ea_rms_T_60N, ea_rms_P_60N,"filled","MarkerFaceColor",colors(8,:))
%scatter(rms_T_predicted, rms_P_predicted,"filled","MarkerFaceColor",colors(2,:),'Marker','square')
hold off
grid on
axis([0 1.2 0 0.7])
xlabel('r.m.s. temperature change [^\circ C]')
ylabel('r.m.s. precipitation change [mm/day]')

ax2 = subplot(1,2,2);
ax2.Position = [0.45 0.15 0.5 0.8];
hold on
scatter(ea_rms_T_base, ea_rms_P_base,'k',"filled")
scatter(ea_rms_T_def, ea_rms_P_def,"filled","MarkerFaceColor",colors(2,:),"SizeData",200,'MarkerEdgeColor','k')
scatter(ea_rms_T_EQ, ea_rms_P_EQ,"filled","MarkerFaceColor",colors(3,:),"SizeData",200,'MarkerEdgeColor','k')
scatter(ea_rms_T_15N_15S, ea_rms_P_15N_15S,"filled","MarkerFaceColor",colors(4,:),"SizeData",200,'MarkerEdgeColor','k')
scatter(ea_rms_T_30N_30S, ea_rms_P_30N_30S,"filled","MarkerFaceColor",colors(5,:),"SizeData",200,'MarkerEdgeColor','k')
scatter(ea_rms_T_60N_60S, ea_rms_P_60N_60S,"filled","MarkerFaceColor",colors(6,:),"SizeData",200,'MarkerEdgeColor','k')
scatter(ea_rms_T_30N, ea_rms_P_30N,"filled","MarkerFaceColor",colors(7,:),"SizeData",200,'MarkerEdgeColor','k')
scatter(ea_rms_T_60N, ea_rms_P_60N,"filled","MarkerFaceColor",colors(8,:),"SizeData",200,'MarkerEdgeColor','k')
scatter(rms_T_predicted, rms_P_predicted,"filled","MarkerFaceColor",colors(2,:),'MarkerEdgeColor','k','Marker','square')
scatter(temp_residual_rms_5DOF, prec_residual_rms_5DOF,"filled","MarkerFaceColor",[0.5 0.5 0.5],'MarkerEdgeColor','k','Marker','square')
plot(TvsP_temp_residual_rms_arrs(1,:), TvsP_prec_residual_rms_arrs(1,:),'k','HandleVisibility','off')
scatter(TvsP_temp_residual_rms_arrs(1,:), TvsP_prec_residual_rms_arrs(1,:),'k','marker','^','HandleVisibility','off')
plot(TvsP_temp_residual_rms_arrs(2,:), TvsP_prec_residual_rms_arrs(2,:),'k','HandleVisibility','off')
scatter(TvsP_temp_residual_rms_arrs(2,:), TvsP_prec_residual_rms_arrs(2,:),'k','marker','^','HandleVisibility','off')
plot(TvsP_temp_residual_rms_arrs(3,:), TvsP_prec_residual_rms_arrs(3,:),'k','HandleVisibility','off')
scatter(TvsP_temp_residual_rms_arrs(3,:), TvsP_prec_residual_rms_arrs(3,:),'k','marker','^','HandleVisibility','off')
% plot(TvsP_temp_residual_rms_arr_2(1,:), TvsP_prec_residual_rms_arr_2(1,:),'k','HandleVisibility','off')
% scatter(TvsP_temp_residual_rms_arr_2(1,:), TvsP_prec_residual_rms_arr_2(1,:),'k','marker','^','HandleVisibility','off')
% plot(TvsP_temp_residual_rms_arr_2(2,:), TvsP_prec_residual_rms_arr_2(2,:),'k','HandleVisibility','off')
% scatter(TvsP_temp_residual_rms_arr_2(2,:), TvsP_prec_residual_rms_arr_2(2,:),'k','marker','^','HandleVisibility','off')
% plot(TvsP_temp_residual_rms_arr_2(3,:), TvsP_prec_residual_rms_arr_2(3,:),'k','HandleVisibility','off')
% scatter(TvsP_temp_residual_rms_arr_2(3,:), TvsP_prec_residual_rms_arr_2(3,:),'k','marker','^','HandleVisibility','off')
plot(TvsP_temp_residual_rms_arr, TvsP_prec_residual_rms_arr,'k','HandleVisibility','off')
scatter(TvsP_temp_residual_rms_arr, TvsP_prec_residual_rms_arr,'k','filled','marker','square')
for i = 1:3
    scatter(rms_T_def(i), rms_P_def(i),"filled","MarkerFaceColor",colors(2,:),'MarkerEdgeColor','k')
    scatter(rms_T_EQ(i), rms_P_EQ(i),"filled","MarkerFaceColor",colors(3,:),'MarkerEdgeColor','k')
    scatter(rms_T_15N_15S(i), rms_P_15N_15S(i),"filled","MarkerFaceColor",colors(4,:),'MarkerEdgeColor','k')
    scatter(rms_T_30N_30S(i), rms_P_30N_30S(i),"filled","MarkerFaceColor",colors(5,:),'MarkerEdgeColor','k')
    scatter(rms_T_60N_60S(i), rms_P_60N_60S(i),"filled","MarkerFaceColor",colors(6,:),'MarkerEdgeColor','k')
end
hold off
grid on
axis([0.2 0.6 0.2 0.45])
xlabel('r.m.s. temperature change [^\circ C]')
%ylabel('r.m.s. precipitation change [mm/day]')
legend('SSP2-4.5','Multi-Objective','EQ','15N+15S', '30N+30S', '60N+60S','30N','60N',['Multi-Objective' newline 'predicted'],'5DOF optimized','Pareto front','Location','eastoutside')

%% plot RMS T vs P (all members) with pareto /AND/ RMS T vs P for the strategies, each member w its own subplot
figure('Renderer', 'painters', 'Position', [100 200 1320 350])
clf
ax1 = subplot(1,3,1);
ax1.Position = [0.07 0.15 0.2 0.8];
hold on
plot(TvsP_temp_residual_rms_arrs(1,:), TvsP_prec_residual_rms_arrs(1,:),'k','HandleVisibility','off')
scatter(TvsP_temp_residual_rms_arrs(1,:), TvsP_prec_residual_rms_arrs(1,:),'k','marker','^','HandleVisibility','off')
scatter(rms_T_def(1), rms_P_def(1),"filled","MarkerFaceColor",colors(2,:),'MarkerEdgeColor','k')
scatter(rms_T_EQ(1), rms_P_EQ(1),"filled","MarkerFaceColor",colors(3,:),'MarkerEdgeColor','k')
scatter(rms_T_15N_15S(1), rms_P_15N_15S(1),"filled","MarkerFaceColor",colors(4,:),'MarkerEdgeColor','k')
scatter(rms_T_30N_30S(1), rms_P_30N_30S(1),"filled","MarkerFaceColor",colors(5,:),'MarkerEdgeColor','k')
scatter(rms_T_60N_60S(1), rms_P_60N_60S(1),"filled","MarkerFaceColor",colors(6,:),'MarkerEdgeColor','k')
hold off
grid on
axis([0.2 0.6 0.2 0.45])
xlabel('r.m.s. temperature change [^\circ C]')
ylabel('r.m.s. precipitation change [mm/day]')

ax2 = subplot(1,3,2);
ax2.Position = [0.37 0.15 0.2 0.8];
hold on
plot(TvsP_temp_residual_rms_arrs(2,:), TvsP_prec_residual_rms_arrs(2,:),'k','HandleVisibility','off')
scatter(TvsP_temp_residual_rms_arrs(2,:), TvsP_prec_residual_rms_arrs(2,:),'k','marker','^','HandleVisibility','off')
scatter(rms_T_def(2), rms_P_def(2),"filled","MarkerFaceColor",colors(2,:),'MarkerEdgeColor','k')
scatter(rms_T_EQ(2), rms_P_EQ(2),"filled","MarkerFaceColor",colors(3,:),'MarkerEdgeColor','k')
scatter(rms_T_15N_15S(2), rms_P_15N_15S(2),"filled","MarkerFaceColor",colors(4,:),'MarkerEdgeColor','k')
scatter(rms_T_30N_30S(2), rms_P_30N_30S(2),"filled","MarkerFaceColor",colors(5,:),'MarkerEdgeColor','k')
scatter(rms_T_60N_60S(2), rms_P_60N_60S(2),"filled","MarkerFaceColor",colors(6,:),'MarkerEdgeColor','k')
hold off
grid on
axis([0.2 0.6 0.2 0.45])
xlabel('r.m.s. temperature change [^\circ C]')
ylabel('r.m.s. precipitation change [mm/day]')

ax3 = subplot(1,3,3);
ax3.Position = [0.67 0.15 0.32 0.8];
hold on
plot(TvsP_temp_residual_rms_arrs(3,:), TvsP_prec_residual_rms_arrs(3,:),'k','HandleVisibility','off')
scatter(TvsP_temp_residual_rms_arrs(3,:), TvsP_prec_residual_rms_arrs(3,:),'k','marker','^')
scatter(rms_T_def(3), rms_P_def(3),"filled","MarkerFaceColor",colors(2,:),'MarkerEdgeColor','k')
scatter(rms_T_EQ(3), rms_P_EQ(3),"filled","MarkerFaceColor",colors(3,:),'MarkerEdgeColor','k')
scatter(rms_T_15N_15S(3), rms_P_15N_15S(3),"filled","MarkerFaceColor",colors(4,:),'MarkerEdgeColor','k')
scatter(rms_T_30N_30S(3), rms_P_30N_30S(3),"filled","MarkerFaceColor",colors(5,:),'MarkerEdgeColor','k')
scatter(rms_T_60N_60S(3), rms_P_60N_60S(3),"filled","MarkerFaceColor",colors(6,:),'MarkerEdgeColor','k')
hold off
grid on
axis([0.2 0.6 0.2 0.45])
xlabel('r.m.s. temperature change [^\circ C]')
ylabel('r.m.s. precipitation change [mm/day]')
legend('Pareto front','Multi-Objective','EQ','15N+15S', '30N+30S', '60N+60S','Location','eastoutside')

%% plot RMS T vs P for the strategies, individual members w pareto front, normalized per standard error
figure('Renderer', 'painters', 'Position', [100 200 650 400])
clf
ax = axes;
hold on
plot(TvsP_norm_temp_residual_rms_arr, TvsP_norm_prec_residual_rms_arr,'k','HandleVisibility','off')
scatter(TvsP_norm_temp_residual_rms_arr, TvsP_norm_prec_residual_rms_arr,'k','filled','marker','square','HandleVisibility','off')
scatter(ea_rms_T_def_norm, ea_rms_P_def_norm,"filled","MarkerFaceColor",colors(2,:),"SizeData",200,'MarkerEdgeColor','k')
scatter(ea_rms_T_EQ_norm, ea_rms_P_EQ_norm,"filled","MarkerFaceColor",colors(3,:),"SizeData",200,'MarkerEdgeColor','k')
scatter(ea_rms_T_15N_15S_norm, ea_rms_P_15N_15S_norm,"filled","MarkerFaceColor",colors(4,:),"SizeData",200,'MarkerEdgeColor','k')
scatter(ea_rms_T_30N_30S_norm, ea_rms_P_30N_30S_norm,"filled","MarkerFaceColor",colors(5,:),"SizeData",200,'MarkerEdgeColor','k')
scatter(ea_rms_T_60N_60S_norm, ea_rms_P_60N_60S_norm,"filled","MarkerFaceColor",colors(6,:),"SizeData",200,'MarkerEdgeColor','k')
scatter(rms_T_norm_predicted, rms_P_norm_predicted,"filled","MarkerFaceColor",colors(2,:),'MarkerEdgeColor','k','Marker','square')
for i = 1:3
    scatter(rms_T_def_norm(i), rms_P_def_norm(i),"filled","MarkerFaceColor",colors(2,:),'MarkerEdgeColor','k')
    scatter(rms_T_EQ_norm(i), rms_P_EQ_norm(i),"filled","MarkerFaceColor",colors(3,:),'MarkerEdgeColor','k')
    scatter(rms_T_15N_15S_norm(i), rms_P_15N_15S_norm(i),"filled","MarkerFaceColor",colors(4,:),'MarkerEdgeColor','k')
    scatter(rms_T_30N_30S_norm(i), rms_P_30N_30S_norm(i),"filled","MarkerFaceColor",colors(5,:),'MarkerEdgeColor','k')
    scatter(rms_T_60N_60S_norm(i), rms_P_60N_60S_norm(i),"filled","MarkerFaceColor",colors(6,:),'MarkerEdgeColor','k')
end
hold off
grid on
xlabel('r.m.s. temperature change [^\circ C]')
ylabel('r.m.s. precipitation change [mm/day]')
legend('Multi-Objective','EQ','15N+15S', '30N+30S', '60N+60S',['Multi-Objective' newline 'predicted'],'Location','eastoutside')

%% plot RMS T vs P (all members) with predicted, individual members, and no-EQ pareto front
figure('Renderer', 'painters', 'Position', [100 200 700 400])
clf
ax = axes;
hold on
plot(TvsP_temp_residual_rms_arr, TvsP_prec_residual_rms_arr,'k','HandleVisibility','off')
scatter(TvsP_temp_residual_rms_arr, TvsP_prec_residual_rms_arr,'k','filled','marker','square','HandleVisibility','off')
plot(temp_residual_rms_arr_noEQ, prec_residual_rms_arr_noEQ,'k','HandleVisibility','off')
scatter(temp_residual_rms_arr_noEQ, prec_residual_rms_arr_noEQ,'k','filled','marker','pentagram','HandleVisibility','off')
scatter(ea_rms_T_def, ea_rms_P_def,"filled","MarkerFaceColor",colors(2,:),"SizeData",200,'MarkerEdgeColor','k')
scatter(ea_rms_T_EQ, ea_rms_P_EQ,"filled","MarkerFaceColor",colors(3,:),"SizeData",200,'MarkerEdgeColor','k')
scatter(ea_rms_T_15N_15S, ea_rms_P_15N_15S,"filled","MarkerFaceColor",colors(4,:),"SizeData",200,'MarkerEdgeColor','k')
scatter(ea_rms_T_30N_30S, ea_rms_P_30N_30S,"filled","MarkerFaceColor",colors(5,:),"SizeData",200,'MarkerEdgeColor','k')
scatter(ea_rms_T_60N_60S, ea_rms_P_60N_60S,"filled","MarkerFaceColor",colors(6,:),"SizeData",200,'MarkerEdgeColor','k')
scatter(rms_T_predicted, rms_P_predicted,"filled","MarkerFaceColor",colors(2,:),'Marker','square')
for i = 1:3
    scatter(rms_T_def(i), rms_P_def(i),"filled","MarkerFaceColor",colors(2,:),'MarkerEdgeColor','k')
    scatter(rms_T_EQ(i), rms_P_EQ(i),"filled","MarkerFaceColor",colors(3,:),'MarkerEdgeColor','k')
    scatter(rms_T_15N_15S(i), rms_P_15N_15S(i),"filled","MarkerFaceColor",colors(4,:),'MarkerEdgeColor','k')
    scatter(rms_T_30N_30S(i), rms_P_30N_30S(i),"filled","MarkerFaceColor",colors(5,:),'MarkerEdgeColor','k')
    scatter(rms_T_60N_60S(i), rms_P_60N_60S(i),"filled","MarkerFaceColor",colors(6,:),'MarkerEdgeColor','k')
end
hold off
grid on
xlabel('r.m.s. temperature change [^\circ C]')
ylabel('r.m.s. precipitation change [mm/day]')
legend('Multi-Objective','EQ','15N+15S', '30N+30S', '60N+60S','Multi-Objective predicted', 'Location', 'eastoutside')

%% plot RMS T vs P (all members) with predicted and no-EQ pareto front
figure('Renderer', 'painters', 'Position', [100 200 560 320])
clf
ax = axes;
hold on
plot(TvsP_temp_residual_rms_arr, TvsP_prec_residual_rms_arr,'k','HandleVisibility','off')
scatter(TvsP_temp_residual_rms_arr, TvsP_prec_residual_rms_arr,'k','filled','marker','square','HandleVisibility','off')
plot(temp_residual_rms_arr_noEQ, prec_residual_rms_arr_noEQ,'k','HandleVisibility','off')
scatter(temp_residual_rms_arr_noEQ, prec_residual_rms_arr_noEQ,'k','filled','marker','pentagram','HandleVisibility','off')
scatter(ea_rms_T_def, ea_rms_P_def,"filled","MarkerFaceColor",colors(2,:))
scatter(ea_rms_T_EQ, ea_rms_P_EQ,"filled","MarkerFaceColor",colors(3,:))
scatter(ea_rms_T_15N_15S, ea_rms_P_15N_15S,"filled","MarkerFaceColor",colors(4,:))
scatter(ea_rms_T_30N_30S, ea_rms_P_30N_30S,"filled","MarkerFaceColor",colors(5,:))
scatter(ea_rms_T_60N_60S, ea_rms_P_60N_60S,"filled","MarkerFaceColor",colors(6,:))
scatter(rms_T_predicted, rms_P_predicted,"filled","MarkerFaceColor",colors(2,:),'Marker','square')
hold off
grid on
xlabel('r.m.s. temperature change [^\circ C]')
ylabel('r.m.s. precipitation change [mm/day]')
legend('Multi-Objective','EQ','15N+15S', '30N+30S', '60N+60S','Multi-Objective predicted','Location','eastoutside')

%% plot RMS land T vs land P-E (all members) pareto
figure('Renderer', 'painters', 'Position', [100 200 420 460])
clf
ax = axes;
hold on
plot(landT_landPE_temp_residual_rms_arr, landT_landPE_PE_residual_rms_arr,'k','HandleVisibility','off')
scatter(landT_landPE_temp_residual_rms_arr, landT_landPE_PE_residual_rms_arr,'k','filled','marker','square','HandleVisibility','off')
scatter(ea_rms_land_T_def, ea_rms_land_PE_def,"filled","MarkerFaceColor",colors(2,:))
scatter(ea_rms_land_T_EQ, ea_rms_land_PE_EQ,"filled","MarkerFaceColor",colors(3,:))
scatter(ea_rms_land_T_15N_15S, ea_rms_land_PE_15N_15S,"filled","MarkerFaceColor",colors(4,:))
scatter(ea_rms_land_T_30N_30S, ea_rms_land_PE_30N_30S,"filled","MarkerFaceColor",colors(5,:))
scatter(ea_rms_land_T_60N_60S, ea_rms_land_PE_60N_60S,"filled","MarkerFaceColor",colors(6,:))
%scatter(rms_T_predicted, rms_P_predicted,"filled","MarkerFaceColor",colors(2,:),'Marker','square')
hold off
grid on
xlabel('r.m.s. land T change [^\circ C]')
ylabel('r.m.s. land P-E change [mm/day]')
legend('Multi-Objective','EQ','15N+15S', '30N+30S', '60N+60S','Location','southoutside')

%% plot RMS land T vs land P-E normalized per std (all members) with predicted
figure('Renderer', 'painters', 'Position', [100 200 420 500])
clf
ax = axes;
hold on
plot(landT_landPE_norm_temp_residual_rms_arr, landT_landPE_norm_PE_residual_rms_arr,'k','HandleVisibility','off')
scatter(landT_landPE_norm_temp_residual_rms_arr, landT_landPE_norm_PE_residual_rms_arr,'k','filled','marker','square','HandleVisibility','off')
scatter(ea_rms_land_T_def_norm, ea_rms_land_PE_def_norm,"filled","SizeData",200,"MarkerFaceColor",colors(2,:),'MarkerEdgeColor','k')
scatter(ea_rms_land_T_EQ_norm, ea_rms_land_PE_EQ_norm,"filled","SizeData",200,"MarkerFaceColor",colors(3,:),'MarkerEdgeColor','k')
scatter(ea_rms_land_T_15N_15S_norm, ea_rms_land_PE_15N_15S_norm,"filled","SizeData",200,"MarkerFaceColor",colors(4,:),'MarkerEdgeColor','k')
scatter(ea_rms_land_T_30N_30S_norm, ea_rms_land_PE_30N_30S_norm,"filled","SizeData",200,"MarkerFaceColor",colors(5,:),'MarkerEdgeColor','k')
scatter(ea_rms_land_T_60N_60S_norm, ea_rms_land_PE_60N_60S_norm,"filled","SizeData",200,"MarkerFaceColor",colors(6,:),'MarkerEdgeColor','k')
scatter(rms_land_T_norm_predicted, rms_land_PE_norm_predicted,"filled","MarkerFaceColor",colors(2,:),'Marker','square')
for i = 1:3
    scatter(rms_land_T_def_norm(i), rms_land_PE_def_norm(i),"filled","MarkerFaceColor",colors(2,:),'MarkerEdgeColor','k')
    scatter(rms_land_T_EQ_norm(i), rms_land_PE_EQ_norm(i),"filled","MarkerFaceColor",colors(3,:),'MarkerEdgeColor','k')
    scatter(rms_land_T_15N_15S_norm(i), rms_land_PE_15N_15S_norm(i),"filled","MarkerFaceColor",colors(4,:),'MarkerEdgeColor','k')
    scatter(rms_land_T_30N_30S_norm(i), rms_land_PE_30N_30S_norm(i),"filled","MarkerFaceColor",colors(5,:),'MarkerEdgeColor','k')
    scatter(rms_land_T_60N_60S_norm(i), rms_land_PE_60N_60S_norm(i),"filled","MarkerFaceColor",colors(6,:),'MarkerEdgeColor','k')
end
hold off
grid on
%axis([.25 .85 .2 .32])
xlabel('r.m.s. land T change [\sigma]')
ylabel('r.m.s. land P-E change [\sigma]')
legend('Multi-Objective','EQ','15N+15S', '30N+30S', '60N+60S',['Multi-Objective' newline 'predicted'],'Location','southoutside')

%% plot RMS land T vs land P-E (all members) with predicted
figure('Renderer', 'painters', 'Position', [100 200 420 500])
clf
ax = axes;
hold on
plot(landT_landPE_temp_residual_rms_arr, landT_landPE_PE_residual_rms_arr,'k','HandleVisibility','off')
scatter(landT_landPE_temp_residual_rms_arr, landT_landPE_PE_residual_rms_arr,'k','filled','marker','square','HandleVisibility','off')
scatter(ea_rms_land_T_def, ea_rms_land_PE_def,"filled","MarkerFaceColor",colors(2,:))
scatter(ea_rms_land_T_EQ, ea_rms_land_PE_EQ,"filled","MarkerFaceColor",colors(3,:))
scatter(ea_rms_land_T_15N_15S, ea_rms_land_PE_15N_15S,"filled","MarkerFaceColor",colors(4,:))
scatter(ea_rms_land_T_30N_30S, ea_rms_land_PE_30N_30S,"filled","MarkerFaceColor",colors(5,:))
scatter(ea_rms_land_T_60N_60S, ea_rms_land_PE_60N_60S,"filled","MarkerFaceColor",colors(6,:))
scatter(rms_land_T_predicted, rms_land_PE_predicted,"filled","MarkerFaceColor",colors(2,:),'Marker','square')
hold off
grid on
xlabel('r.m.s. land T change [^\circ C]')
ylabel('r.m.s. land P-E change [mm/day]')
legend('Multi-Objective','EQ','15N+15S', '30N+30S', '60N+60S',['Multi-Objective' newline 'predicted'],'Location','southoutside')

%% plot RMS land T vs land P-E for the strategies, individual members w pareto front
figure('Renderer', 'painters', 'Position', [100 200 450 550])
clf
ax = axes;
hold on
plot(landT_landPE_temp_residual_rms_arr, landT_landPE_PE_residual_rms_arr,'k','HandleVisibility','off')
scatter(landT_landPE_temp_residual_rms_arr, landT_landPE_PE_residual_rms_arr,'k','filled','marker','square','HandleVisibility','off')
scatter(ea_rms_land_T_def, ea_rms_land_PE_def,"filled","MarkerFaceColor",colors(2,:),"SizeData",200,'MarkerEdgeColor','k')
scatter(ea_rms_land_T_EQ, ea_rms_land_PE_EQ,"filled","MarkerFaceColor",colors(3,:),"SizeData",200,'MarkerEdgeColor','k')
scatter(ea_rms_land_T_15N_15S, ea_rms_land_PE_15N_15S,"filled","MarkerFaceColor",colors(4,:),"SizeData",200,'MarkerEdgeColor','k')
scatter(ea_rms_land_T_30N_30S, ea_rms_land_PE_30N_30S,"filled","MarkerFaceColor",colors(5,:),"SizeData",200,'MarkerEdgeColor','k')
scatter(ea_rms_land_T_60N_60S, ea_rms_land_PE_60N_60S,"filled","MarkerFaceColor",colors(6,:),"SizeData",200,'MarkerEdgeColor','k')
scatter(rms_land_T_predicted, rms_land_PE_predicted,"filled","MarkerFaceColor",colors(2,:),'Marker','square','MarkerEdgeColor','k')
for i = 1:3
    scatter(rms_land_T_def(i), rms_land_PE_def(i),"filled","MarkerFaceColor",colors(2,:),'MarkerEdgeColor','k')
    scatter(rms_land_T_EQ(i), rms_land_PE_EQ(i),"filled","MarkerFaceColor",colors(3,:),'MarkerEdgeColor','k')
    scatter(rms_land_T_15N_15S(i), rms_land_PE_15N_15S(i),"filled","MarkerFaceColor",colors(4,:),'MarkerEdgeColor','k')
    scatter(rms_land_T_30N_30S(i), rms_land_PE_30N_30S(i),"filled","MarkerFaceColor",colors(5,:),'MarkerEdgeColor','k')
    scatter(rms_land_T_60N_60S(i), rms_land_PE_60N_60S(i),"filled","MarkerFaceColor",colors(6,:),'MarkerEdgeColor','k')
end
hold off
grid on
xlabel('r.m.s. land T change [^\circ C]')
ylabel('r.m.s. land P-E change [mm/day]')
legend('Multi-Objective','EQ','15N+15S', '30N+30S', '60N+60S',['Multi-Objective' newline 'predicted'],'Location','southoutside')
%% plot injection rates throughout the pareto front rms T vs P AND land T vs land P-E
figure('Renderer', 'painters', 'Position', [100 200 1160 420])
clf
ax1 = subplot(1,2,1);
ax1.Position = [.07 .1 .38 .75];
hold on
plot(1:11,TvsP_inj_each_lat_arr(1,:),"Color",colors(6,:),"Marker","square","MarkerFaceColor","auto",'LineWidth',2) %60N
plot(1:11,TvsP_inj_each_lat_arr(2,:),"Color",colors(5,:),"Marker","square","MarkerFaceColor","auto",'LineWidth',2) %30N
plot(1:11,TvsP_inj_each_lat_arr(3,:),"Color",colors(4,:),"Marker","square","MarkerFaceColor","auto",'LineWidth',2) %15N
plot(1:11,TvsP_inj_each_lat_arr(4,:),"Color",colors(3,:),"Marker","square","MarkerFaceColor","auto",'LineWidth',2) %EQ
plot(1:11,TvsP_inj_each_lat_arr(5,:),'--',"Color",colors(4,:),"Marker","square","MarkerFaceColor","auto",'LineWidth',2) %15S
plot(1:11,TvsP_inj_each_lat_arr(6,:),'--',"Color",colors(5,:),"Marker","square","MarkerFaceColor","auto",'LineWidth',2) %30S
plot(1:11,TvsP_inj_each_lat_arr(7,:),'--',"Color",colors(6,:),"Marker","square","MarkerFaceColor","auto",'LineWidth',2) %60S
hold off
grid on
axis([0 12 0 11])
set(gca,'FontSize',12)
ylabel('Injection rate [Tg/yr]')
xticklabels([])
xlabel('<- r.m.s. P                       r.m.s T ->')
title('r.m.s. T vs r.m.s P')

ax2 = subplot(1,2,2);
ax2.Position = [.52 .1 .46 .75];
hold on
plot(1:11,landT_landPE_inj_each_lat_arr(1,:),"Color",colors(6,:),"Marker","square","MarkerFaceColor","auto",'LineWidth',2) %60N
plot(1:11,landT_landPE_inj_each_lat_arr(2,:),"Color",colors(5,:),"Marker","square","MarkerFaceColor","auto",'LineWidth',2) %30N
plot(1:11,landT_landPE_inj_each_lat_arr(3,:),"Color",colors(4,:),"Marker","square","MarkerFaceColor","auto",'LineWidth',2) %15N
plot(1:11,landT_landPE_inj_each_lat_arr(4,:),"Color",colors(3,:),"Marker","square","MarkerFaceColor","auto",'LineWidth',2) %EQ
plot(1:11,landT_landPE_inj_each_lat_arr(5,:),'--',"Color",colors(4,:),"Marker","square","MarkerFaceColor","auto",'LineWidth',2) %15S
plot(1:11,landT_landPE_inj_each_lat_arr(6,:),'--',"Color",colors(5,:),"Marker","square","MarkerFaceColor","auto",'LineWidth',2) %30S
plot(1:11,landT_landPE_inj_each_lat_arr(7,:),'--',"Color",colors(6,:),"Marker","square","MarkerFaceColor","auto",'LineWidth',2) %60S
hold off
grid on
axis([0 12 0 11])
set(gca,'FontSize',12)
xticklabels([])
xlabel('<- r.m.s. land P-E                       r.m.s land T ->')
legend('60N','30N','15N','EQ','15S','30S','60S','Location','eastoutside')
title('r.m.s. land T vs r.m.s land P-E')

sgtitle('Injection Rates at Different Latitudes')

%% plot injection rates throughout the pareto front rms T vs P
%figure('Renderer', 'painters', 'Position', [100 200 760 420])
figure('Renderer', 'painters', 'Position', [100 200 760 560])
clf
ax = axes;
hold on
plot(1:11,TvsP_inj_each_lat_arr(1,:),"Color",colors(6,:),"Marker","square","MarkerFaceColor","auto",'LineWidth',2) %60N
plot(1:11,TvsP_inj_each_lat_arr(2,:),"Color",colors(5,:),"Marker","square","MarkerFaceColor","auto",'LineWidth',2) %30N
plot(1:11,TvsP_inj_each_lat_arr(3,:),"Color",colors(4,:),"Marker","square","MarkerFaceColor","auto",'LineWidth',2) %15N
plot(1:11,TvsP_inj_each_lat_arr(4,:),"Color",colors(3,:),"Marker","square","MarkerFaceColor","auto",'LineWidth',2) %EQ
plot(1:11,TvsP_inj_each_lat_arr(5,:),'--',"Color",colors(4,:),"Marker","^","MarkerFaceColor","auto",'LineWidth',2) %15S
plot(1:11,TvsP_inj_each_lat_arr(6,:),'--',"Color",colors(5,:),"Marker","^","MarkerFaceColor","auto",'LineWidth',2) %30S
plot(1:11,TvsP_inj_each_lat_arr(7,:),'--',"Color",colors(6,:),"Marker","^","MarkerFaceColor","auto",'LineWidth',2) %60S
plot(1:11,sum(TvsP_inj_each_lat_arr),'k',"Marker","o","MarkerFaceColor","auto",'LineWidth',3) %total
hold off
grid on
set(gca,'FontSize',12)
ylabel('Injection rate [Tg/yr]')
xticklabels([])
xlabel('<- r.m.s. P                       r.m.s T ->')
legend('60N','30N','15N','EQ','15S','30S','60S','total','Location','eastoutside')
title({'Injection Rates at Different Latitudes', 'r.m.s. T vs r.m.s P'})

%% plot injection rates throughout the pareto front rms T vs P normalized per standard deviation
figure('Renderer', 'painters', 'Position', [100 200 760 420])
clf
ax = axes;
hold on
plot(1:11,TvsP_norm_inj_each_lat_arr(1,:),"Color",colors(6,:),"Marker","square","MarkerFaceColor","auto",'LineWidth',2) %60N
plot(1:11,TvsP_norm_inj_each_lat_arr(2,:),"Color",colors(5,:),"Marker","square","MarkerFaceColor","auto",'LineWidth',2) %30N
plot(1:11,TvsP_norm_inj_each_lat_arr(3,:),"Color",colors(4,:),"Marker","square","MarkerFaceColor","auto",'LineWidth',2) %15N
plot(1:11,TvsP_norm_inj_each_lat_arr(4,:),"Color",colors(3,:),"Marker","square","MarkerFaceColor","auto",'LineWidth',2) %EQ
plot(1:11,TvsP_norm_inj_each_lat_arr(5,:),'--',"Color",colors(4,:),"Marker","^","MarkerFaceColor","auto",'LineWidth',2) %15S
plot(1:11,TvsP_norm_inj_each_lat_arr(6,:),'--',"Color",colors(5,:),"Marker","^","MarkerFaceColor","auto",'LineWidth',2) %30S
plot(1:11,TvsP_norm_inj_each_lat_arr(7,:),'--',"Color",colors(6,:),"Marker","^","MarkerFaceColor","auto",'LineWidth',2) %60S
plot(1:11,sum(TvsP_norm_inj_each_lat_arr),'k',"Marker","o","MarkerFaceColor","auto",'LineWidth',3) %total
hold off
grid on
set(gca,'FontSize',12)
ylabel('Injection rate [Tg/yr]')
xticklabels([])
xlabel('<- r.m.s. P                       r.m.s T ->')
legend('60N','30N','15N','EQ','15S','30S','60S','total','Location','eastoutside')
title({'Injection Rates at Different Latitudes', 'r.m.s. T vs r.m.s P normalized per std'})

%% plot injection rates throughout the pareto front rms landT vs landPE
%figure('Renderer', 'painters', 'Position', [100 200 760 420])
figure('Renderer', 'painters', 'Position', [100 200 760 560])
clf
ax = axes;
hold on
plot(1:11,landT_landPE_inj_each_lat_arr(1,:),"Color",colors(6,:),"Marker","square","MarkerFaceColor","auto",'LineWidth',2) %60N
plot(1:11,landT_landPE_inj_each_lat_arr(2,:),"Color",colors(5,:),"Marker","square","MarkerFaceColor","auto",'LineWidth',2) %30N
plot(1:11,landT_landPE_inj_each_lat_arr(3,:),"Color",colors(4,:),"Marker","square","MarkerFaceColor","auto",'LineWidth',2) %15N
plot(1:11,landT_landPE_inj_each_lat_arr(4,:),"Color",colors(3,:),"Marker","square","MarkerFaceColor","auto",'LineWidth',2) %EQ
plot(1:11,landT_landPE_inj_each_lat_arr(5,:),'--',"Color",colors(4,:),"Marker","^","MarkerFaceColor","auto",'LineWidth',2) %15S
plot(1:11,landT_landPE_inj_each_lat_arr(6,:),'--',"Color",colors(5,:),"Marker","^","MarkerFaceColor","auto",'LineWidth',2) %30S
plot(1:11,landT_landPE_inj_each_lat_arr(7,:),'--',"Color",colors(6,:),"Marker","^","MarkerFaceColor","auto",'LineWidth',2) %60S
plot(1:11,sum(landT_landPE_inj_each_lat_arr),'k',"Marker","o","MarkerFaceColor","auto",'LineWidth',3) %total
hold off
grid on
set(gca,'FontSize',12)
ylabel('Injection rate [Tg/yr]')
xticklabels([])
xlabel('<- r.m.s. land P-E                       r.m.s land T ->')
legend('60N','30N','15N','EQ','15S','30S','60S','total','Location','eastoutside')
title({'Injection Rates at Different Latitudes', 'r.m.s. land T vs r.m.s land P-E'})

%% plot injection rates throughout the pareto front rms landT vs landPE normalized per standard deviation
figure('Renderer', 'painters', 'Position', [100 200 760 420])
clf
ax = axes;
hold on
plot(1:11,landT_landPE_norm_inj_each_lat_arr(1,:),"Color",colors(6,:),"Marker","square","MarkerFaceColor","auto",'LineWidth',2) %60N
plot(1:11,landT_landPE_norm_inj_each_lat_arr(2,:),"Color",colors(5,:),"Marker","square","MarkerFaceColor","auto",'LineWidth',2) %30N
plot(1:11,landT_landPE_norm_inj_each_lat_arr(3,:),"Color",colors(4,:),"Marker","square","MarkerFaceColor","auto",'LineWidth',2) %15N
plot(1:11,landT_landPE_norm_inj_each_lat_arr(4,:),"Color",colors(3,:),"Marker","square","MarkerFaceColor","auto",'LineWidth',2) %EQ
plot(1:11,landT_landPE_norm_inj_each_lat_arr(5,:),'--',"Color",colors(4,:),"Marker","^","MarkerFaceColor","auto",'LineWidth',2) %15S
plot(1:11,landT_landPE_norm_inj_each_lat_arr(6,:),'--',"Color",colors(5,:),"Marker","^","MarkerFaceColor","auto",'LineWidth',2) %30S
plot(1:11,landT_landPE_norm_inj_each_lat_arr(7,:),'--',"Color",colors(6,:),"Marker","^","MarkerFaceColor","auto",'LineWidth',2) %60S
plot(1:11,sum(landT_landPE_norm_inj_each_lat_arr),'k',"Marker","o","MarkerFaceColor","auto",'LineWidth',3) %total
hold off
grid on
set(gca,'FontSize',12)
ylabel('Injection rate [Tg/yr]')
xticklabels([])
xlabel('<- r.m.s. land P-E                       r.m.s land T ->')
legend('60N','30N','15N','EQ','15S','30S','60S','total','Location','eastoutside')
title({'Injection Rates at Different Latitudes', 'r.m.s. land T vs r.m.s land P-E normalized per std'})

%% plot injection rates throughout the pareto front rms T vs P NO EQ vs w EQ
figure('Renderer', 'painters', 'Position', [100 200 1160 560])
clf
ax1 = subplot(1,2,1);
ax1.Position = [.07 .1 .38 .75];
hold on
plot(1:11,inj_each_lat_arr_noEQ(1,:),"Color",colors(6,:),"Marker","square","MarkerFaceColor","auto",'LineWidth',2) %60N
plot(1:11,inj_each_lat_arr_noEQ(2,:),"Color",colors(5,:),"Marker","square","MarkerFaceColor","auto",'LineWidth',2) %30N
plot(1:11,inj_each_lat_arr_noEQ(3,:),"Color",colors(4,:),"Marker","square","MarkerFaceColor","auto",'LineWidth',2) %15N
plot(1:11,inj_each_lat_arr_noEQ(4,:),'--',"Color",colors(4,:),"Marker","^","MarkerFaceColor","auto",'LineWidth',2) %15S
plot(1:11,inj_each_lat_arr_noEQ(5,:),'--',"Color",colors(5,:),"Marker","^","MarkerFaceColor","auto",'LineWidth',2) %30S
plot(1:11,inj_each_lat_arr_noEQ(6,:),'--',"Color",colors(6,:),"Marker","^","MarkerFaceColor","auto",'LineWidth',2) %60S
plot(1:11,sum(inj_each_lat_arr_noEQ),'k',"Marker","o","MarkerFaceColor","auto",'LineWidth',3) %total
hold off
grid on
axis([0 12 0 18])
set(gca,'FontSize',12)
ylabel('Injection rate [Tg/yr]')
xticklabels([])
xlabel('<- r.m.s. P                       r.m.s T ->')
title('no EQ')

ax2 = subplot(1,2,2);
ax2.Position = [.52 .1 .46 .75];
hold on
plot(1:11,TvsP_inj_each_lat_arr(1,:),"Color",colors(6,:),"Marker","square","MarkerFaceColor","auto",'LineWidth',2) %60N
plot(1:11,TvsP_inj_each_lat_arr(2,:),"Color",colors(5,:),"Marker","square","MarkerFaceColor","auto",'LineWidth',2) %30N
plot(1:11,TvsP_inj_each_lat_arr(3,:),"Color",colors(4,:),"Marker","square","MarkerFaceColor","auto",'LineWidth',2) %15N
plot(1:11,TvsP_inj_each_lat_arr(4,:),"Color",colors(3,:),"Marker","square","MarkerFaceColor","auto",'LineWidth',2) %EQ
plot(1:11,TvsP_inj_each_lat_arr(5,:),'--',"Color",colors(4,:),"Marker","^","MarkerFaceColor","auto",'LineWidth',2) %15S
plot(1:11,TvsP_inj_each_lat_arr(6,:),'--',"Color",colors(5,:),"Marker","^","MarkerFaceColor","auto",'LineWidth',2) %30S
plot(1:11,TvsP_inj_each_lat_arr(7,:),'--',"Color",colors(6,:),"Marker","^","MarkerFaceColor","auto",'LineWidth',2) %60S
plot(1:11,sum(TvsP_inj_each_lat_arr),'k',"Marker","o","MarkerFaceColor","auto",'LineWidth',3) %total
hold off
grid on
axis([0 12 0 18])
set(gca,'FontSize',12)
ylabel('Injection rate [Tg/yr]')
xticklabels([])
xlabel('<- r.m.s. P                       r.m.s T ->')
legend('60N','30N','15N','EQ','15S','30S','60S','total','Location','eastoutside')
title('with EQ')

sgtitle('Injection Rates at Different Latitudes')