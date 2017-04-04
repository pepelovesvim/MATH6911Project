%% the actual approximation of any input via rounding to nearest node
function f = optionPrice(S0,K,T,r,sigma,N_time,N_space,xmax)
i = ceil((1/(2*xmax/(N_space))*((K/S0)+xmax)));
A = Asian_CN(T,r,sigma,N_time,N_space,xmax);
f = S0*A(i+1);
end 