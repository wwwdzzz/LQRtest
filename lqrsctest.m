close all
x0=trajs3{1}.zt(1:tcdim,1)

P=Uuptd;
Nstep=size(trajs3{1}.zt,2)
Atil=Atilde;
Btil=Btilde;
prenewtraj = LQRsc(DS,z0,'nonlinear',tspan,Nstep,auData,W_0,masterModes,cont,P,Atil,Btil,tcdim,datalength,trajs3{1}.zauto)
figure(60)
hold on
plot(prenewtraj.us(1,:),'r')
plot(prenewtraj.mlus(1,:),'b')
plot(prenewtraj.xauto(1,:),'g')
plot(prenewtraj.xt(1,:))

Q2      = DS.epsilon^2*transpose(P)*cont.Q(1:tcdim,1:tcdim)*P
[K,S,eigs]=dlqr(Atil,Btil,Q2,cont.Rhat);
