function traj = scLQR0auto(P,Atil,Btil,tcdim,datalength,X0,cont,type)

Nstep=5000;%datalength;



%need to add Qf ,DS.btilde, DS.atilde
if type==1%for simple problem
    Q = eye(2);%cont.Q(1:tcdim,1:tcdim);
    Rhat = eye(1);%cont.Rhat;
    Qf = eye(2);%cont.Qf(1:tcdim,1:tcdim);
    ep   =1;
else
    Q = cont.Q(1:50,1:50);
    Rhat = 0.0001*cont.Rhat;
    Qf = cont.Qf(1:50,1:50);
    ep   =1;
end
%tspan = 

% auto part

xauto=zeros(tcdim,Nstep);
%xauto=xauto(1:tcdim,1:Nstep);


%% The process to get P matrix
%P = [];%%%%%%%%Util
% z0 = Pξ0

%xi0   = P'*z0tr;%%%%%%%
xi0   = P'*X0;


%% calculate offline matrix
bQ_N    = 2*ep*transpose(Qf*P)*xauto(:,end);%7
bQ_k    = 2*ep*transpose(Q*P)*xauto;
Q2      = ep^2*transpose(P)*Q*P; %50x50
Qfhat   = ep^2*transpose(P)*Qf*P; % 50x50

%Atil = DS.Atil;
%Btil = DS.Btil;

%%matlab LQR 
usmatlab       = zeros(size(Btil,2),Nstep);
ximatlab      = zeros(size(xi0,1), Nstep);
ximatlab(:,1) = xi0;
[K,P1,~]=dlqr(Atil,Btil,Q2,Rhat);
for j = 1:Nstep-1
    usmatlab(:,j) = - K* ximatlab(:,j);%%%8
    ximatlab(:,j+1) = Atil * ximatlab(:,j) + Btil * usmatlab(:,j);
end






M_Nm1    = (Rhat + Btil' * Qfhat * Btil)^-1;
h_Nm1    = -bQ_N' * Btil * M_Nm1 * Btil' * Qfhat' * Atil;%0
K_Nm1    = -Atil' * Qfhat * Btil * M_Nm1 * Btil' * Qfhat' * Atil;%%%%%%%%%%%%有问题
Qtil_Nm1 = Q2 + Atil' * Qfhat * Atil;%%%%%%%%%%%%%%
bN_Nm1   = bQ_k(:,end-1)' + bQ_N' * Atil;%0

Mks    = cell(1,Nstep);
hks    = cell(1,Nstep);
Kks    = cell(1,Nstep);
Qtilks = cell(1,Nstep);
bNks   = cell(1,Nstep);

ans11=cell(1,Nstep);

Mks{Nstep}    = M_Nm1;
hks{Nstep}    = h_Nm1;
Kks{Nstep}    = K_Nm1;%%%%%%%%%%%%%%%%
Qtilks{Nstep} = Qtil_Nm1;%%%%%%%%%%%%%%
bNks{Nstep}   = bN_Nm1;
ans11{Nstep}   = K_Nm1+Qtil_Nm1;


for k = 5000-1:-1:2
    % Mk h(k) Kk Qtilk bNk
    % whether they will converge ?
    Mks{k}   = (Rhat + Btil' * (Qtilks{k+1} + Kks{k+1}) * Btil)^-1;%sl
    hks{k}   = -(bNks{k+1}+hks{k+1})*Btil*Mks{k}*Btil'*(Qtilks{k+1}+Kks{k+1})'*Atil;
    Kks{k}   = -Atil'*(Qtilks{k+1}+Kks{k+1})*Btil*Mks{k}*Btil'*(Qtilks{k+1}+Kks{k+1})'*Atil;
    Qtilks{k}= Q2 + Atil'*(Qtilks{k+1}+Kks{k+1})*Atil;
    bNks{k}  = bQ_k(:,k)' + (bNks{k+1}+hks{k+1}) * Atil;

    % Mks{k}   = (Rhat + Btil' * (P1) * Btil)^-1;%sl
    % hks{k}   = -(bNks{k+1}+hks{k+1})*Btil*Mks{k}*Btil'*(P1)'*Atil;
    % Kks{k}   = -Atil'*(P1)*Btil*Mks{k}*Btil'*(P1)'*Atil;
    % Qtilks{k}= Q2 + Atil'*(P1)*Atil;
    % bNks{k}  = bQ_k(:,k)' + (bNks{k+1}+hks{k+1}) * Atil;
end
k = k - 1;
%Mks{k} = (Rhat + Btil' * (Qtilks{k+1} + Kks{k+1}) * Btil)^-1;   % M0
Mks{k} = (Rhat + Btil' * (P1) * Btil)^-1;

%% calculate control policy
us       = zeros(size(Btil,2),Nstep);
xis      = zeros(size(xi0,1), Nstep);
xis(:,1) = xi0;
Kesc=cell(1,Nstep);

for j = 1:Nstep-1

    %us(:,j) = - Mks{j} * ((0.5*(bNks{j+1}+hks{j+1})*Btil + xis(:,j)'*Atil'*(Qtilks{j+1}+Kks{j+1})*Btil)');%%%8
    %us(:,j) = - Mks{j} * Btil'*(Qtilks{j+1}+Kks{j+1})'*Atil*xis(:,j);
   
    ans11{j}=Qtilks{j+1}+Kks{j+1};
    Kesc{j}=Mks{j} * Btil'*(ans11{j})'*Atil;%多了个转置？
    us(:,j)=-Kesc{j}*xis(:,j);
    xis(:,j+1) = Atil * xis(:,j) + Btil * us(:,j);
end
j = j + 1;
us(:,j) = -Mks{j}*((0.5*bQ_N'*Btil + xis(:,j)'*Atil'*Qfhat*Btil)');%30







%% evalute state of original system
%
% xt     = xauto+ep*real(P*xis);
xt     =real(P*xis);
outdof = tcdim;
xout   = xt%(outdof,:);
% invariance PDE error
% res   = invariance_PDE_residual(DS,tspan,zt,ut);
% ratio = ratio_ext2int(DS,tspan,zt,ut);
%% output
traj      = struct();

traj.us   = us;
%traj.etat = etat;
%traj.xit  = xit;
traj.xt   = xout;
traj.xauto=xauto;
traj.mlus=usmatlab;
traj.mlxt=real(P*ximatlab);

end


