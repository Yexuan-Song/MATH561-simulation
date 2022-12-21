%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [stats,plotdata]=simulation(runstuff,params,protoc)

num_days = runstuff.maxDays; %number of days 
qu_days = params.qu_days; %quaraintine days
popA_size = params.popA_size; %population size in location A
popB_size = params.popB_size; %population size in location B
gamma = params.gamma; %wanning rate from state R to state S
n = params.n; %number of infectious in the population to start the pandemic

beta_muA=params.beta_muA_baseline; % the basic beta, location A
beta_muB=params.beta_muB_baseline; % the basic beta, location B
beta_k=params.beta_k_baseline; % dispersion
vuln = params.vuln; %vulnerable individuals

dh_A = params.dh_A; %death rate in A (hospital)
dh_B = params.dh_B; %death rate in B (hospital)

mig_A = protoc.mig_A; %migration rate from A to B
mig_B = protoc.mig_B; %migration rate from B to A
eps_A = protoc.eps_A; %error rate for false negative test in A, or even testing
eps_B = protoc.eps_B; %error rate for false negative test in B, or even testing

huge = 1e8; %to aviod confusion when comparing inequalities

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%initializing 
%set up the number of cases 
%SIR model
%initializing the infected states 
%A_S,B_S;A_SQ,B_SQ S class
%A_I,B_I;A_IQ,B_IQ I class
%A_R,B_R;A_RQ,B_RQ Q class
%A_H,B_H Hospital class

%intuition make the vector of size popA+popB in that vector size will
%change due to migration events, but the total is not changing

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%S class
size = popA_size+popB_size;
%This helps to track which individual is in location A or location B
location = zeros(popA_size+popB_size,1);
pos = popA_size+1:popA_size+popB_size;
location(pos)=1;  %0 means location A, 1 means location B


%S class
S=ones(popA_size+popB_size,1);
SQ = zeros(popA_size+popB_size,1);
time_SQ=huge*ones(popA_size+popB_size,1);

%I class
I = zeros(popA_size+popB_size,1);
IQ = zeros(popA_size+popB_size,1);
time_I=huge*ones(popA_size+popB_size,1);
time_IQ=huge*ones(popA_size+popB_size,1);

%R class
R = zeros(popA_size+popB_size,1);
RQ = zeros(popA_size+popB_size,1);
time_R=huge*ones(popA_size+popB_size,1);
time_RQ=huge*ones(popA_size+popB_size,1);

%H class
H = zeros(popA_size+popB_size,1);
time_H=huge*ones(popA_size+popB_size,1);

%quarantine time
Q=zeros(popA_size+popB_size,1);

%flag if tested positive in the quanrantine
flag=zeros(popA_size+popB_size,1);
tt_test=zeros(popA_size+popB_size,1);

%death
D=zeros(popA_size+popB_size,1);

% this keeps track of what generation of infection is
% -1 is not infected yet
generation=-1*ones(popA_size+popB_size,1);

%set beta for each individual, need to be careful about the position since
%infecitious class is the total number of people but beta is location, 
%need 1+popA_szie to line the index for location B
beta_individualA=gamma_rate(beta_muA,beta_k,popA_size); 
beta_individualB=gamma_rate(beta_muB,beta_k,popB_size);
%combine beta
beta_individual=[beta_individualA;beta_individualB];

re_rate = gamma_rate(gamma,5,popA_size+popB_size);

%setup infectious period for each individual
%%%%%%%%%%%%%%%
infectiousperiod=mygamma(5,2,popA_size+popB_size);
m_A=mygamma(mig_A,0.1,popA_size);
m_B=mygamma(mig_B,0.1,popB_size);
migration=[m_A;m_B];

%who is vulnerable, vulnerbale individual (elder) go to
%the hospital if they get infected 
vuln=rand(popA_size+popB_size,1)<vuln;

was_infected=zeros(popA_size+popB_size,1); 
was_infected_locA=zeros(popA_size+popB_size,1);
was_infected_locB=zeros(popA_size+popB_size,1);
was_hospitalized_locA=zeros(popA_size+popB_size,1);
was_hospitalized_locB=zeros(popA_size+popB_size,1);
who_infected=zeros(popA_size+popB_size,1); 

ind_case = randperm(popA_size,n); %find out infectee 
vuln(ind_case,1)=0; %vuln not the index cases


%REMEMBER when consider location B, the starting point is
%from popA_size+1 to popA_size+popB_size

%initial infection 
%consider the pandemic happens in location A
S(ind_case,1)=0;
I(ind_case,1)=1;
generation(ind_case,1)=0;
time_I(ind_case)=0;
was_infected(ind_case,1)=1;
was_infected_locA(ind_case,1)=1;
who_infected(ind_case,1)=-1;
I_locA=zeros(size,1);
I_locB=zeros(size,1);
I_locA(ind_case,1)=1;
hos_A=zeros(size,1);
hos_B=zeros(size,1);

numsteps = num_days;
S_mat=zeros(size,numsteps+1);
SQ_mat=zeros(size,numsteps+1);
I_mat=zeros(size,numsteps+1);
I_mat_locA=zeros(size,numsteps+1);
I_mat_locB=zeros(size,numsteps+1);
IQ_mat=zeros(size,numsteps+1);
R_mat=zeros(size,numsteps+1);
RQ_mat=zeros(size,numsteps+1);
H_mat=zeros(size,numsteps+1);
l_mat=zeros(size,numsteps+1);
hosA_mat=zeros(size,numsteps+1);
hosB_mat=zeros(size,numsteps+1);

%matrix to determine how long for quarantine, this is constant not a
%distribtion
Q_mat=zeros(size,numsteps+1);

%initial state
S_mat(:,1)=S;
SQ_mat(:,1)=SQ;
I_mat(:,1)=I;
IQ_mat(:,1)=IQ;
R_mat(:,1)=R;
RQ_mat(:,1)=RQ;
H_mat(:,1)=H;
Q_mat(:,1)=Q;
l_mat(:,1)=location;
I_mat_locA(:,1)=I_locA;
I_mat_locB(:,1)=I_locB;
hosA_mat(:,1)=hos_A;
hosB_mat(:,1)=hos_B;
%start the loooooooooooooop
time=0; % check the time when disease gose to location B

for kk=1:numsteps
    %loop over time
    currentTime=kk;
    %loop all individuals might get infected  
    for k=1:size
            %random number to decide if migration to other place or not
            %now loop individuals who are not t
            rand_mig=rand;
            %loop individuals who might infect others 
            if D(k,1)==0
                %check location
                loc=location(k,1);
             for jj=1:size
                %random number to decide if get infected or not
                rand_num=rand;
                %check the location

                %is k susceptible
                if S(k,1)==1
                    %check if going to migrate
                    if  rand_mig<migration(k,1)
                        S(k,1)=0;
                        SQ(k,1)=1;
                        time_SQ(k,1)=currentTime;
                    else
                        %if not going to migrate, loop over infectious
                        if (location(jj,1)==loc && I(jj,1)==1)
                            beta = beta_individual(jj);
                        else
                            beta = 0;
                        end
                            if rand_num < beta
                                S(k,1)=0;
                                I(k,1)=1;
                                time_I(k,1)=currentTime;
                                generation(k,1)=generation(jj,1)+1;
                                was_infected(k,1)=1;
                                who_infected(k,1)=jj;
                                if loc==0
                                    was_infected_locA(k,1)=1;
                                    I_locA(k,1)=1;
                                elseif loc==1
                                    was_infected_locB(k,1)=1;
                                    I_locB(k,1)=1;
                                end
                            end
                    end
                end
            end

            if SQ(k,1)==1
                if (currentTime > time_SQ(k,1)+qu_days)
                %are we past the quarantine period?
                SQ(k,1)=0;
                S(k,1)=1;
                        if loc==0
                            location(k,1)=1;
                        elseif loc==1
                            location(k,1)=0;
                        end
                %migrate back by less chance
                migration(k,1)=migration(k,1)*0.5;
                end
            end

            if I(k,1)==1
                if rand_mig<migration(k,1)
                    I(k,1)=0;
                    IQ(k,1)=1;
                    time_IQ(k,1)=currentTime;
                else
                    if currentTime>time_I(k,1)+infectiousperiod(k,1)
                        I(k,1)=0;
                        R(k,1)=1;
                        time_R(k,1)=currentTime;
                            I_locA(k,1)=0;
                            I_locB(k,1)=0;
                    else
                        %check vunerable
                        if vuln(k,1)==1
                            I(k,1)=0;
                            H(k,1)=1;
                            time_H(k,1)=currentTime;
                            I_locA(k,1)=0;
                            I_locB(k,1)=0;
                            if loc==0
                                was_hospitalized_locA(k,1)=1;
                                hos_A(k,1)=1;
                            else
                                was_hospitalized_locB(k,1)=1;
                                hos_B(k,1)=1;
                            end
                        end
                    end
                end
            end

            if IQ(k,1)==1
                %first test if recovered or not 
                %if recovered before testing, then good
                %if not, testing may result in more Q time 
                if vuln(k,1)==1
                        IQ(k,1)=0;
                        H(k,1)=1;
                        time_H(k,1)=currentTime;
                        if loc==0
                            was_hospitalized_locB(k,1)=1;
                            hos_B(k,1)=1;
                        else
                            was_hospitalized_locA(k,1)=1;
                            hos_A(k,1)=1;
                        end
                end
                if currentTime>time_I(k,1)+infectiousperiod(k,1)
                    IQ(k,1)=0;
                    R(k,1)=1;
                        I_locB(k,1)=0;
                        I_locA(k,1)=0;
                        if loc==0
                            location(k,1)=1;
                        elseif loc==1
                            location(k,1)=0;
                        end
                else
                    test_err=rand;
                    if currentTime > time_IQ(k,1)+qu_days && IQ(k,1)==1
                    if flag(k,1)==0
                        if loc==0 %at location A, will go to locB
                            if test_err < eps_B
                                IQ(k,1)=0;
                                I(k,1)=1;
                                location(k,1)=1;
                                I_locB(k,1)=1;
                                I_locA(k,1)=0;
                            else
                            %tested positive, stay in the hotel
                                flag(k,1)=1;
                                tt_test(k,1)=1;
                            end
                        elseif loc==1 %at location B, will go to locA
                            if test_err < eps_A
                                IQ(k,1)=0;
                                I(k,1)=1;
                                location(k,1)=0;
                                I_locB(k,1)=0;
                                I_locA(k,1)=1;
                            else
                                %tested positive, stay in the hotel
                                flag(k,1)=1;
                                tt_test(k,1)=1;
                            end
                        end
                migration(k,1)=migration(k,1)*0.5;
                    end
                    end
                end
            end

            if R(k,1)==1
                if rand_mig <migration(k,1)
                    R(k,1)=0;
                    RQ(k,1)=1;
                    time_RQ(k,1)=currentTime;
                else 
                    if currentTime > time_H(k,1)+re_rate(k,1)
                        %wanning immunity, return S
                        R(k,1)=0;
                        S(k,1)=1;
                    end
                end
            end

            if RQ(k,1)==1
                if currentTime > time_RQ(k,1)+qu_days
                    RQ(k,1)=0;
                    R(k,1)=1;
                    if loc==0
                        %change location 
                        location(k,1)=1;
                    elseif loc==1
                        location(k,1)=0;
                    end
                   migration(k,1)=migration(k,1)*0.5;
                end
            end

            if H(k,1)==1
                rad=rand;
                if currentTime > time_H(k,1)+infectiousperiod(k,1)
                    H(k,1)=0;
                    R(k,1)=1;
                    time_R(k,1)=currentTime;
                    hos_A(k,1)=0;
                    hos_B(k,1)=0;
                else
                if loc==0
                    if rad<dh_A
                        H(k,1)=0;
                        D(k,1)=1;
                        hos_A(k,1)=0;
                    end
                elseif loc==1
                    if rad<dh_B
                    H(k,1)=0;
                    D(k,1)=1;
                    hos_B(k,1)=0;
                    end
                end
                end
            end
            end
    end
    if sum(I_locB)>0
        time=time+1;
    end
    S_mat(:,kk+1)=S;
    SQ_mat(:,kk+1)=SQ;
    I_mat(:,kk+1)=I;
    IQ_mat(:,kk+1)=IQ;
    R_mat(:,kk+1)=R;
    RQ_mat(:,kk+1)=RQ;
    H_mat(:,kk+1)=H;
    Q_mat(:,kk+1)=Q;
    l_mat(:,kk+1)=location;
    I_mat_locA(:,kk+1)=I_locA;
    I_mat_locB(:,kk+1)=I_locB;
    hosA_mat(:,kk+1)=hos_A;
    hosB_mat(:,kk+1)=hos_B;
end

who_died = sum(D);
tot_infected=sum(was_infected);
tot_tested=sum(tt_test);
tot_inf_locA=sum(was_infected_locA);
tot_inf_locB=sum(was_infected_locB);
tot_hos_locA=sum(was_hospitalized_locA);
tot_hos_locB=sum(was_hospitalized_locB);

stats.total_infected=tot_infected;
stats.total_dead=who_died;
stats.total_tested=tot_tested;
stats.total_inf_locA=tot_inf_locA;
stats.total_inf_locB=tot_inf_locB;
stats.total_hos_locA=tot_hos_locA;
stats.total_hos_locB=tot_hos_locB;
stats.time=time;

plotdata.S=S_mat;
plotdata.SQ=SQ_mat;
plotdata.I=I_mat;
plotdata.IQ=IQ_mat;
plotdata.R=R_mat;
plotdata.RQ=RQ_mat;
plotdata.H=H_mat;
plotdata.l=l_mat;
plotdata.IA=I_mat_locA;
plotdata.IB=I_mat_locB;
plotdata.hA=hosA_mat;
plotdata.hB=hosB_mat;

end
