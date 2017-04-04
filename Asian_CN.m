function f=Asian_CN(T, r,sigma,N_time,N_space,xmax)
%**********************************************************
%C-N Scheme for Equation(1) with the transform T-t=tau.
%This guarentees 
%**********************************************************
delta_t=T/N_time; %step in time variable
delta_x=2*xmax/N_space; %step in space variable
%**********************************************************
phi_old=zeros(N_space-1,1);
phi_new=phi_old;

%**********************************************************
% initial condition
for i = 0:N_space-2
    phi_old(i+1) = max(-(-xmax + i*delta_x), 0);
end
%**********************************************************
% creating generator matrix L: matrix corresponding to discretized operator
L=zeros(N_space-1);
for i=1:(N_space-1)
    L(i,i)=(-sigma^2*(-xmax+i*delta_x)^2)/(delta_x)^2;
    if (i==1)
        L(i,i+1)=-(1/(2*delta_x))*(1/T + r*(-xmax+i*delta_x))+ 0.5*sigma^2*(-xmax + i*delta_x)^2/(delta_x)^2;
    elseif (i==N_space-1)
        L(i,i-1)=(1/(2*delta_x))*(1/T + r*(-xmax+i*delta_x))+ 0.5*sigma^2*(-xmax + i*delta_x)^2/(delta_x)^2;
    else
        L(i,i-1)=(1/(2*delta_x))*(1/T + r*(-xmax+i*delta_x))+ 0.5*sigma^2*(-xmax + i*delta_x)^2/(delta_x)^2;
        L(i,i+1)=-(1/(2*delta_x))*(1/T + r*(-xmax+i*delta_x))+ 0.5*sigma^2*(-xmax + i*delta_x)^2/(delta_x)^2;
    end;
end;
%**********************************************************
D=sparse(2*eye(N_space-1)+delta_t*L); 
C=sparse(2*eye(N_space-1)-delta_t*L);
A=inv(C)*D; %iteration matrix in Crank-Nicolson scheme 
%**********************************************************
% iterating in time
for j=1:N_time
    phi_new=A*phi_old;
    % update
    phi_old=phi_new;
end;
%**********************************************************
% return the result
f=phi_new;