function [reLambda] = DispersRel(tau, a, b)

[f,g] = funs(tau,a,b,10);

values = roots(f,g,'marchingsquares');
real_componnents = values(:,1);

reLambda = max(real_componnents);


end

