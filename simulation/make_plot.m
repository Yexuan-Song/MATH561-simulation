clear;clc;
S1=load('plot7.mat');

ndays=101;
ntry=50;


%ci = nonpar(S1,ndays,ntry);
mu = meaninf(S1,ndays,ntry);
mu_inf = num_inf(S1,ntry);
mu_death=num_death(S1,ntry);
mu_locA=num_inf_A(S1,ndays,ntry);
mu_locB=num_inf_B(S1,ndays,ntry);


x = 1:ndays;

%days = 1:9:ndays-1;
%labels=1:90;
figure 
plot(x,mu_locA)

title(sprintf('Mean total infected %.2f, mean death %.2f',mu_inf,mu_death))
xlabel('days')
ylabel('cases')

%set(gca,'XTick',days);
%set(gca,'XTickLabel',labels);
hold on
%plot(x,mu_locA);
plot(x,mu_locB);
% plot(x, yCI95+yMean,'-r')                           % Plot 95% Confidence Intervals Of All Experiments
%patch([x, fliplr(x)], [ci(:,1)' fliplr(ci(:,2)')], 'b', 'EdgeColor','none', 'FaceAlpha',0.25)
hold off





function ci = nonpar(data,ndays,ntry)
    ci = zeros(ndays,2);
    for j = 1:ndays
    vec = [];
        for i = 1:ntry
        vec(i) = data.big_multi{i}(j);
        end
    sorted = sort(vec);
    lower = round(0.025*ntry+1);
    upper = round(0.975*ntry-1);
    ci(j,1) = sorted(lower);
    ci(j,2) = sorted(upper);
    end  
end

function mu = meaninf(data,ndays,ntry)
    mu = [];
    for j = 1:ndays
    val = 0;
    for i = 1:ntry
    val = val + data.big_multi{i}(j);
    end
    mu(j) = val/ntry;
    end
end


function mu_inf = num_inf(data,ntry)
    mu_inf = 0;

    for i=1:ntry
    ss  = data.big_multi(i,2);
    ss = ss{1};
    mu_inf = mu_inf+ss.total_infected;
    end
    mu_inf = mu_inf/ntry;
end


function mu_death = num_death(data,ntry)
    mu_death = 0;

    for i=1:ntry
    ss  = data.big_multi(i,2);
    ss = ss{1};
    mu_death = mu_death+ss.total_dead;
    end
    mu_death = mu_death/ntry;
end

function mu_locA = num_inf_A(data,ndays,ntry)
    mu_locA=[];
    for j=1:ndays
        val=0;
        for i=1:ntry
            ss=data.big_multi(i,4);
            ss=ss{1};
            val=val+ss(j);
        end
        mu_locA(j)=val/ntry;
    end
end

function mu_locB = num_inf_B(data,ndays,ntry)
    mu_locB=[];
    for j=1:ndays
        val=0;
        for i=1:ntry
            ss=data.big_multi(i,5);
            ss=ss{1};
            val=val+ss(j);
        end
        mu_locB(j)=val/ntry;
    end
end
