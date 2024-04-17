function u = closed_loop_control_lqr(t,q,Rhat,Bhat,Ptfun,stfun)
% CLOSED_LOOP_CONTROL_LQR This function calculates closed-loop optimal
% control input via a time-varying feedback matrix obtained from Riccati
% differential equation.
%
% u(t) = -R^(-1)*B.'*P(t)*q+0.5*R^(-1)*B.'*s(t)

tmp = Rhat\Bhat.';
nt  = numel(t); nq = size(Bhat,1);
if nt==1
    %u = (-tmp*reshape(Ptfun(t),[nq,nq])*q+0.5*tmp*stfun(t));
    u = (-tmp*reshape(Ptfun(t),[nq,nq])*q*1.7+0.1*tmp*stfun(t));
else
    nc = size(Rhat,1);
    u  = zeros(nc,nt);
    Pt = Ptfun(t); st = stfun(t);
    for k=1:nt
        u(:,k) = (-tmp*reshape(Pt(k,:),[nq,nq])*q(:,k)*1.7+0.1*tmp*st(:,k));
    end
end
u = real(u);
end