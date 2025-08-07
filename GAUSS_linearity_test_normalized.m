%% Simulation information
n_base = 3; %Number of baseline (SSP2-4.5) ensemble members
n_def = 3;
n_EQ = 3;
n_15N_15S = 3;
n_30N_30S = 3;
n_60N_60S = 3;
n_30N = 3;
n_30S = 3;
n_15N = 3;
n_15S = 3;
n_60N = 1;

start_base = 2000; %Number of baseline (SSP2-4.5) ensemble members
start_def = 2035;
start_EQ = 2035;
start_15N_15S = 2035;
start_30N_30S = 2035;
start_60N_60S = 2035;
start_30N = 2035;
start_30S = 2035;
start_15N = 2035;
start_15S = 2035;
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
%scenarios = {'base','def','EQ','15N_15S','30N_30S','60N_60S','30N','30S','15N','15S','60N'};
scenarios = {'base','def','15N_15S','30N_30S','30N','30S','15N','15S'};

colorscheme = [ 0  0  0;
               .6  0 .8;
               .8 .1  0;
               .9 .6  0;
               .3 .6 .2;
                0  0 .8;
                0  1  0;
                0 .8  1]; % for plotting the different strategies/scenarios

%% load data
load('C:\Users\ezrab\Documents\Cornell\Research\MATLAB data\GAUSS.mat');
load('C:\Users\ezrab\Documents\Cornell\Research\MATLAB data\GAUSS_V_zm.mat');
scenarios = {'base','def','15N_15S','30N_30S','30N','30S','15N','15S'};

%% MY standard error
addpath('C:\Users\ezrab\Documents\Cornell\Research\Code\SAI_optimization_paper-main\SAI_optimization_paper-main')

PE_std_err = calc_std_error(squeeze(mean(reshape(PRECT_base(:,:,97:336,:),288,192,12,20,[]),3))-squeeze(mean(reshape(QFLX_base(:,:,97:336,:),288,192,12,20,[]),3)));
P_std_err = calc_std_error(squeeze(mean(reshape(PRECT_base(:,:,97:336,:),288,192,12,20,[]),3)));
T_std_err = calc_std_error(squeeze(mean(reshape(TREFHT_base(:,:,97:336,:),288,192,12,20,[]),3)));

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
varnames = {'PRECT', 'TREFHT', 'QFLX'}; % Edit this line to choose which variable to process

for i = 1:length(scenarios)
    for k=1:length(varnames)
        eval(['[temporal_mean_' varnames{k} '_' scenarios{i} ', ea_temporal_mean_' varnames{k} '_' scenarios{i} ',~] = take_temporal_means(' varnames{k} '_' scenarios{i} ', start_' scenarios{i} ', comp_start, comp_end);']);
    end
end

%for reference period
for k=1:length(varnames)
    eval(['[temporal_mean_' varnames{k} '_ref, ea_temporal_mean_' varnames{k} '_ref,~] = take_temporal_means(' varnames{k} '_base, start_base, ref_start, ref_end);']);
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
    eval(['ITCZ_' scenarios{k} ' = zeros(35,2);'])
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

%% calculate alphas
%alpha is defined as the change in a variable due to GHG warming
%per degree of GHG warming
for k = 1:length(varnames)
    eval(['alpha_' varnames{k} ' = zeros(289,192,n_base);'])
    for i = 1:n_base
        eval(['alpha_' varnames{k} '(:,:,i) = (temporal_mean_' varnames{k} '_base(:,:,i) - temporal_mean_' varnames{k} '_ref(:,:,i))/delta_T0(i);'])
    end
end
alpha_SSI = (ea_temporal_mean_SSI_base-ea_temporal_mean_SSI_ref)/mean(SSI_delta_T0);
alpha_T1 = (temporal_mean_T1_base-temporal_mean_T1_ref)./delta_T0;
alpha_T2 = (temporal_mean_T2_base-temporal_mean_T2_ref)./delta_T0;
alpha_ITCZ = (temporal_mean_ITCZ_base-temporal_mean_ITCZ_ref)./delta_T0;

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
        end
    end
end
for j = 2:length(scenarios)
eval(['n = n_' scenarios{j} ';'])
    eval(['cooling_' scenarios{j} ' = zeros(n,1);'])
    eval(['mu_T1_' scenarios{j} ' = zeros(n,1);'])
    eval(['mu_T2_' scenarios{j} ' = zeros(n,1);'])
    eval(['mu_ITCZ_' scenarios{j} ' = zeros(n,1);'])
    for i = 1:n
        eval(['cooling_' scenarios{j} '(i) = mean(global_annual_mean_TREFHT_base(comp_start-start_base+1:comp_end-start_base+1,i)) - mean(global_annual_mean_TREFHT_' scenarios{j} '(comp_start-start_' scenarios{j} '+1:comp_end-start_' scenarios{j} '+1,i));'])
        eval(['mu_T1_' scenarios{j} '(i) = (temporal_mean_T1_' scenarios{j} '(i) - temporal_mean_T1_base(i))/cooling_' scenarios{j} '(i);'])
        eval(['mu_T2_' scenarios{j} '(i) = (temporal_mean_T2_' scenarios{j} '(i) - temporal_mean_T2_base(i))/cooling_' scenarios{j} '(i);'])
        eval(['mu_ITCZ_' scenarios{j} '(i) = (temporal_mean_ITCZ_' scenarios{j} '(i) - temporal_mean_ITCZ_base(i))/cooling_' scenarios{j} '(i);'])
    end
    eval(['ea_mu_SSI_' scenarios{j} ' = ((ea_temporal_mean_SSI_' scenarios{j} ' - ea_temporal_mean_SSI_ref) - alpha_SSI*mean(delta_T0))/mean(cooling_' scenarios{j} ');'])
end
%% get injection rates
%inj_rate_arr is 35(years) by 5(strategies)
% the order of strategies is:
% 0N_1.0, 15N_15S, 30N_30S, 60N_60S_1.0, Global+1C
load('C:\Users\ezrab\Documents\Cornell\Research\MATLAB data\yearly_injection_rate.mat')

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

ea_injections_def = (injections_def_001 + injections_def_002 + injections_def_003)/3;
ea_inj_def_end = mean(ea_injections_def(16:35,:));
ea_inj_EQ_end = mean(inj_rate_arr(16:35,1));
ea_inj_15N_15S_end = mean(inj_rate_arr(16:35,2));
ea_inj_30N_30S_end = mean(inj_rate_arr(16:35,3));
ea_inj_60N_60S_end = mean(inj_rate_arr(16:35,4));

%% Predict Multi-objective change in variables using linear combination of 15N_15S, 30N_30S, 30N, 30S, 15N, 15S
tot = sum(ea_inj_def_end);
a_15N_15S = (ea_inj_def_end(2)+ea_inj_def_end(3))/(tot); % fraction of response from 15N_15S
a_30N_30S = (ea_inj_def_end(1)+ea_inj_def_end(4))/(tot);
a_30N = ea_inj_def_end(4)/tot - a_30N_30S/2;
a_30S = ea_inj_def_end(1)/tot - a_30N_30S/2;
a_15N = ea_inj_def_end(3)/tot - a_15N_15S/2;
a_15S = ea_inj_def_end(2)/tot - a_15N_15S/2;

allVars = {'PRECT', 'TREFHT', 'QFLX'}; %2D vars
for k=1:length(allVars)
    eval(['mu_' allVars{k} '_predicted = a_15N_15S*mean(mu_' allVars{k} '_15N_15S,3) + a_30N_30S*mean(mu_' allVars{k} '_30N_30S,3) +' ...
        'a_30N*mean(mu_' allVars{k} '_30N,3) + a_30S*mean(mu_' allVars{k} '_30S,3) +' ...
        'a_15N*mean(mu_' allVars{k} '_15N,3) + a_15S*mean(mu_' allVars{k} '_15S,3);'])
end

allVars = {'T1', 'T2','ITCZ'}; %scalar vars
for k=1:length(allVars)
    eval(['mu_' allVars{k} '_predicted = a_15N_15S*mean(mu_' allVars{k} '_15N_15S) + a_30N_30S*mean(mu_' allVars{k} '_30N_30S) +' ...
        'a_30N*mean(mu_' allVars{k} '_30N) + a_30S*mean(mu_' allVars{k} '_30S) +' ...
        'a_15N*mean(mu_' allVars{k} '_15N) + a_15S*mean(mu_' allVars{k} '_15S);'])
end
ea_mu_SSI_predicted = a_15N_15S*ea_mu_SSI_15N_15S + a_30N_30S*ea_mu_SSI_30N_30S + ...
        a_30N*ea_mu_SSI_30N + a_30S*ea_mu_SSI_30S + ...
        a_15N*ea_mu_SSI_15N + a_15S*ea_mu_SSI_15S;

save('mu_predicted.mat','mu_TREFHT_predicted','mu_PRECT_predicted','mu_QFLX_predicted','ea_mu_SSI_predicted','mu_T1_predicted','mu_T2_predicted','mu_ITCZ_predicted')
area_weight=AREA/sum(AREA,'all');
rms_T_error = squeeze(sqrt(sum(area_weight.*((mu_TREFHT_predicted(1:288,:)-mean(mu_TREFHT_def(1:288,:,:),3)).^2),[1 2])))
rms_P_error = squeeze(sqrt(sum(area_weight.*((mu_PRECT_predicted(1:288,:)-mean(mu_PRECT_def(1:288,:,:),3)).^2),[1 2])))
%% print ensemble-averaged predicted mu values
%ea_mu_T0_predicted = squeeze(sum(mu_TREFHT_predicted(1:288,:).*AREA,[1 2])/AE)
ea_mu_T1_predicted = mean(mu_T1_predicted)
ea_mu_T2_predicted = mean(mu_T2_predicted)
ea_mu_SSI_predicted
ea_mu_P0_predicted = squeeze(sum(mu_PRECT_predicted(1:288,:).*AREA,[1 2]))/AE
ea_mu_ITCZ_predicted = mean(mu_ITCZ_predicted)

ea_mu_T1_def = mean(mu_T1_def)
ea_mu_T2_def = mean(mu_T2_def)
ea_mu_SSI_def
ea_mu_P0_def = squeeze(sum(mean(mu_PRECT_def(1:288,:,:),3).*AREA,[1 2]))/AE
ea_mu_ITCZ_def = mean(mu_ITCZ_def)

%% create T and P maps for predicted, actual, and difference for multi-objective strategy
Tmin = -3;
Tmax = 3;
Pmin = -2;
Pmax = 2;
[T_contours, T_color_map] = custom_color_map(Tmin,Tmax,25,'*RdBu');
[P_contours, P_color_map] = custom_color_map(Pmin,Pmax,11,'BrBG');

figure('Renderer', 'painters', 'Position', [20 60 800 700])
clf

ax1=subplot(3,2,1); % T predicted
ax1.Position = [.07 .7 .42 .25];
hold on
worldmap([-90 90],[-180 180])
setm(ax1, 'meridianlabel', 'off', 'parallellabel', 'off')
box on
%H = ttest2(temporal_mean_TREFHT_R1, temporal_mean_TREFHT_base, 'Dim',3);
%sH = interp2(lon2,lat,H',slon,slat');
contourfm(lat,lon2,mu_TREFHT_predicted',T_contours,'LineStyle','none')
%[II,JJ]=find(~sH');
%scatterm(slat(JJ),slon(II),1,'k')
load coastlines
geoshow(coastlat,coastlon,'Color','k')
caxis([Tmin Tmax])
hold off
title('\DeltaT Predicted');
colormap(ax1, T_color_map);
add_abc_tofig(1,0.04,-0.02,12);

ax2=subplot(3,2,2); % P predicted
ax2.Position = [.53 .7 .42 .25];
hold on
worldmap([-90 90],[-180 180])
setm(ax2, 'meridianlabel', 'off', 'parallellabel', 'off')
box on
contourfm(lat,lon2,mu_PRECT_predicted',P_contours,'LineStyle','none')
load coastlines
geoshow(coastlat,coastlon,'Color','k')
caxis([Pmin Pmax])
hold off
title('\DeltaP Predicted');
colormap(ax2, P_color_map);
add_abc_tofig(2,0.04,-0.02,12);

ax3=subplot(3,2,3); % T actual
ax3.Position = [.07 .41 .42 .25];
hold on
worldmap([-90 90],[-180 180])
setm(ax3, 'meridianlabel', 'off', 'parallellabel', 'off')
box on
contourfm(lat,lon2,mean(mu_TREFHT_def,3)',T_contours,'LineStyle','none')
load coastlines
geoshow(coastlat,coastlon,'Color','k')
caxis([Tmin Tmax])
hold off
title('\DeltaT Actual');
colormap(ax3, T_color_map);
add_abc_tofig(3,0.04,-0.02,12);

ax4=subplot(3,2,4); % P actual
ax4.Position = [.53 .41 .42 .25];
hold on
worldmap([-90 90],[-180 180])
setm(ax4, 'meridianlabel', 'off', 'parallellabel', 'off')
box on
contourfm(lat,lon2,mean(mu_PRECT_def,3)',P_contours,'LineStyle','none')
load coastlines
geoshow(coastlat,coastlon,'Color','k')
caxis([Pmin Pmax])
hold off
title('\DeltaP Actual');
colormap(ax4, P_color_map);
add_abc_tofig(4,0.04,-0.02,12);

% error check
% temperature
Z=abs((mu_TREFHT_predicted - mean(mu_TREFHT_def,3))./(2*[T_std_err; T_std_err(1,:)]));
[lon_null,lat_null]=find(Z<1);
lon_null=(lon_null-1)/287*360;
lat_null=(lat_null-1)/191*180-90;

ax5=subplot(3,2,5); % T diff
ax5.Position = [.07 .02 .42 .35];
hold on
worldmap([-90 90],[-180 180])
setm(ax5, 'meridianlabel', 'off', 'parallellabel', 'off')
box on
contourfm(lat,lon2,mu_TREFHT_predicted' - mean(mu_TREFHT_def,3)',T_contours,'LineStyle','none')
scatterm(lat_null,lon_null,.2,'k','filled');
load coastlines
geoshow(coastlat,coastlon,'Color','k')
caxis([Tmin Tmax])
cb=colorbar;
ylabel(cb,'K','FontSize',10,'Rotation',0,'FontWeight','bold')
cb.Location = 'southoutside';
hold off
title('\DeltaT Error');
colormap(ax5, T_color_map);
add_abc_tofig(5,0.04,-0.02,12);

% precip
Z=abs((mu_PRECT_predicted - mean(mu_PRECT_def,3))./(2*[P_std_err; P_std_err(1,:)]));
[lon_null,lat_null]=find(Z<1);
lon_null=(lon_null-1)/287*360;
lat_null=(lat_null-1)/191*180-90;

ax6=subplot(3,2,6); % P diff
ax6.Position = [.53 .02 .42 .35];
hold on
worldmap([-90 90],[-180 180])
setm(ax6, 'meridianlabel', 'off', 'parallellabel', 'off')
box on
contourfm(lat,lon2,mu_PRECT_predicted' - mean(mu_PRECT_def,3)',P_contours,'LineStyle','none')
scatterm(lat_null,lon_null,.2,'k','filled');
load coastlines
geoshow(coastlat,coastlon,'Color','k')
caxis([Pmin Pmax])
cb=colorbar;
ylabel(cb,'mm/day','FontSize',10,'Rotation',0,'FontWeight','bold')
cb.Location = 'southoutside';
hold off
title('\DeltaP Error');
colormap(ax6, P_color_map);
add_abc_tofig(6,0.04,-0.02,12);