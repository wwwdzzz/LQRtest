function [Atil,Btil]=DMDcuknB(X,tcdim,datalength,datastart,truncatedim)


% tcdim=100;
% datalength=120;
% datastart=360;
p=tcdim+2;
r=truncatedim;

X=trajs3{1}.zt(1:tcdim,:)-trajs3{1}.zauto(1:tcdim,:);
Xhalf=X(:,datastart:datastart+datalength);
%X=trajs3{1}.zt(1:tcdim,:);
u=0.1*trajs3{1}.ut(:,datastart:datastart+datalength-1);
N=size(Xhalf,2)-1;
timei2=linspace(datastart,datastart+N,N);




X1 = Xhalf(:,1:end-1); 
X2 = Xhalf(:,2:end);
Gama=u(:,1:end);
Omega=[X1;Gama];

[U,S,V] = svd(Omega,'econ');

Utilde = U(:,1:p); 
Stilde = S(1:p,1:p);
Vtilde = V(:,1:p);

U1tilde=Utilde(1:size(X1,1),:);
U2tilde=Utilde(size(X1,1)+1:size(Utilde,1), :);

A_=X2*Vtilde*pinv(Stilde)*(U1tilde.')
B_=X2*Vtilde*pinv(Stilde)*(U2tilde.')%3.17

[Uup,Sup,Vup]=svd(X2,'econ');

 Uuptd = Uup(:,1:r); 
 Suptd = Sup(1:r,1:r);
 Vuptd = Vup(:,1:r);

Atilde=(Uuptd.')*A_*Uuptd%3.18
Btilde=(Uuptd.')*B_
eig(Atilde)

[W,eigs] = eig(Atilde);
phi=X2*Vtilde*inv(Stilde)*(U1tilde.')*Uuptd*W; %3.21
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Atd=A_;
Btd=B_;
Ctd=C;
Atrk=Atilde;
Btrk=Btilde;
Ctrk=C;

xt=zeros(tcdim,N);
yt=zeros(1,N);
ut=zeros(2,N);

epsino=zeros(r,N);
ytrk=zeros(1,N);
utrk=zeros(2,N);

X0=X(1:tcdim,datastart);
xt(:,1)=X0;
epsino(:,1)=pinv(Uuptd)*X0;
for k =1:N
    ut(:,k)=u(:,k);
    xt(:,k+1)=Atd*xt(:,k)+Btd*ut(:,k);
    yt(:,k)=xt(1,k);
    utrk(:,k)=u(:,k);
    epsino(:,k+1)=Atrk*epsino(:,k)+Btrk*ut(:,k);
   
end
xtrkre=Uuptd*epsino;
ytrkre=xtrkre(1,1:datalength);

y=X2(1,:);

figure(20);
plot(timei2,y,'b','DisplayName','origin','Color',[0,0,1]);
hold on;
figure(20);
plot(timei2,yt,'o','DisplayName','DMD','Color',[1,0,0]);
plot(timei2,ytrkre,'b','DisplayName','origin','Color',[0,1,1]);

Atil=Atilde;
Btil=Btilde;

end