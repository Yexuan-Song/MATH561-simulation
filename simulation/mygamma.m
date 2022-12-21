% generate gamma distributed random variable with a given mean and standard
% deviation
function out=mygamma(mea,std,n)

if std~=0
  a=(mea/std)^2;
  b=mea/a;

  out=gamrnd(a,b,n,1);
else
  out=mea*ones(n,1);
end

end
