clear;clc;

runstuff.num_sims = 50;  % number of simulations per choice of params
runstuff.maxDays=100; % how many days each sim
runstuff.seed=123456; % random number seed

big_multi = cell(runstuff.num_sims,5); %big cells to add 

params.qu_days = 3; %quaraintine days

params.popA_size=5000; %population size in location A
params.popB_size=5000; %population size in location B
%alpha = params.alpha; %recovery rate from state I to state R (not include hospitalization)
params.gamma=20; %wanning rate from state R to state S
%h = params.h; %hospitalization rate for both location A and B 
%r_A = params.r_A; %recovery rate from hospital in A
%r_B = params.r_B; %recovery rate from hospital in B
params.n=30; %number of infectious in the population to start the pandemic

params.beta_muA_baseline=0.0001; % the basic beta, location A
params.beta_muB_baseline=0.0002; % the basic beta, location B
params.beta_k_baseline=0.1; % dispersion
params.vuln=0.05; %vulnerable individuals

params.dh_A=0.05; %death rate in A (hospital)
params.dh_B=0.01; %death rate in B (hospital)

protoc.mig_A=0.15; %migration rate from A to B
protoc.mig_B=0.01; %migration rate from B to A
protoc.eps_A=0.1; %error rate for false negative test in A, or even testing
protoc.eps_B=0.2; %error rate for false negative test in B, or even testing

[multi_stats,multi_plot]=multi_plots_stats(runstuff,params,protoc);
                for k=1:runstuff.num_sims

                %multi_stats(k).setting=setting;
                multi_stats(k).popA_size=params.popA_size;
                multi_stats(k).popB_size=params.popB_size;
                multi_stats(k).beta_muA_baseline=params.beta_muA_baseline;
                multi_stats(k).beta_muB_baseline=params.beta_muB_baseline;
                end
           for k=1:runstuff.num_sims
                I = multi_plot{k}.I;
                num_inf = sum(I);
                I_A=multi_plot{k}.IA;
                num_infA=sum(I_A);
                I_B=multi_plot{k}.IB;
                num_infB=sum(I_B);
                big_multi{k,1}=num_inf;
                big_multi{k,2}=multi_stats(k);
                big_multi{k,3}=multi_plot{k};
                big_multi{k,4}=num_infA;
                big_multi{k,5}=num_infB;
           end

save('scen1.mat','big_multi')



