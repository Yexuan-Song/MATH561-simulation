% generate gamma distributed random variable with a given mean and
% dispersion k
% mea= mean of rate
% k = dispersion of rate
% n = number of draws
function out=gamma_rate(mea,k,n)

if k~=Inf
  a=k;
  b=mea/k;

  out=gamrnd(a,b,n,1);
else
  out=mea*ones(n,1);
end

end
