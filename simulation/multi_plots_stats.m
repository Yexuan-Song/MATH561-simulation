function [multi_run_stats,multi_run_plots]=multi_plots_stats(runstuff,params,protoc)

% how many simulations
num_sims=runstuff.num_sims;
% set random number seed
rng(runstuff.seed)

% intialize and make space
multi_run_plots=cell(num_sims,1);


% intialize and make space
multi_run_stats.total_infected=zeros(num_sims,1);
multi_run_stats.total_dead=zeros(num_sims,1);
multi_run_stats.total_tested=zeros(num_sims,1);
multi_run_stats.total_inf_locA=zeros(num_sims,1);
multi_run_stats.total_inf_locB=zeros(num_sims,1);
multi_run_stats.total_hos_locA=zeros(num_sims,1);
multi_run_stats.total_hos_locB=zeros(num_sims,1);
multi_run_stats.time=zeros(num_sims,1);

% for each simulation
for kk=1:num_sims
  
  % simulate class
  [stats,plotdata]=simulation(runstuff,params,protoc);    
  % store results
  multi_run_plots{kk}=plotdata;
  multi_run_stats(kk,1)=stats;
end
end

      
