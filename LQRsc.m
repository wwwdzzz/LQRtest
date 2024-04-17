function traj = LQRsc(DS,x0,proj,tspan,Nstep,autData,Wmap,masterModes,cont,P,Atil,Btil,tcdim,datalength,autopart,X0)

%need to add Qf ,DS.btilde, DS.atilde
Q = cont.Q(1:tcdim,1:tcdim);
Rhat = 0.1*cont.Rhat;
Qf = cont.Qf(1:tcdim,1:tcdim);
ep   = DS.epsilon;
n    = size(DS.M,1);
%tspan = 

% auto part
eta0   = get_initial_p0(DS,masterModes,x0,proj,autData,Wmap);
z0   = (x0-reduced_to_full(eta0,Wmap,[],0))/ep;

nonode = @(t,eta) auto_red_dyn(eta,autData);
option = odeset('RelTol',1e-8,'AbsTol',1e-10);


%forward simulation of reduced dynamics
[~,etat] = ode45(nonode,tspan,eta0,option);
etat    = transpose(etat);
xauto = reduced_to_full_traj([],etat,Wmap); % u=0
%xauto=zeros(tcdim,Nstep);
xauto=xauto(1:tcdim,1:Nstep);


% z0=autopart(:,1);
% xauto=autopart(1:tcdim,1:Nstep);
% %xauto=zeros(tcdim,Nstep);

%% The process to get P matrix
%P = [];%%%%%%%%Util
% z0 = Pξ0

z0tr=z0(1:tcdim,:);
xi0   = P'*z0tr;%%%%%%%

%% calculate offline matrix
bQ_N    = 2*ep*transpose(Qf*P)*xauto(:,end);%7
bQ_k    = 2*ep*transpose(Q*P)*xauto;
Q2      = ep^2*transpose(P)*Q*P; %50x50
Qfhat   = ep^2*transpose(P)*Qf*P; % 50x50

%Atil = DS.Atil;
%Btil = DS.Btil;

M_Nm1    = (Rhat + Btil' * Qfhat * Btil)^-1;
h_Nm1    = -bQ_N' * Btil * M_Nm1 * Btil' * Qfhat' * Atil;%0
K_Nm1    = -Atil' * Qfhat * Btil * M_Nm1 * Btil' * Qfhat' * Atil;%1-5
Qtil_Nm1 = Q2 + Atil' * Qfhat * Atil;
bN_Nm1   = bQ_k(:,end-1)' + bQ_N' * Atil;%0

Mks    = cell(1,Nstep);
hks    = cell(1,Nstep);
Kks    = cell(1,Nstep);
Qtilks = cell(1,Nstep);
bNks   = cell(1,Nstep);

Mks{Nstep}    = M_Nm1;
hks{Nstep}    = h_Nm1;
Kks{Nstep}    = K_Nm1;
Qtilks{Nstep} = Qtil_Nm1;
bNks{Nstep}   = bN_Nm1;



for k = Nstep-1:-1:2
    % Mk h(k) Kk Qtilk bNk
    % whether they will converge ?
    Mks{k}   = (Rhat + Btil' * (Qtilks{k+1} + Kks{k+1}) * Btil)^-1;%sl
    hks{k}   = -(bNks{k+1}+hks{k+1})*Btil*Mks{k}*Btil'*(Qtilks{k+1}+Kks{k+1})'*Atil;
    Kks{k}   = -Atil'*(Qtilks{k+1}+Kks{k+1})*Btil*Mks{k}*Btil'*(Qtilks{k+1}+Kks{k+1})'*Atil;
    Qtilks{k}= Q2 + Atil'*(Qtilks{k+1}+Kks{k+1})*Atil;
    %bNks{k}  = bQ_k(:,k+1)' + bQ_k(:,k)' * Atil;
    bNks{k}  = bQ_k(:,k)' + (bNks{k+1}+hks{k+1}) * Atil;
end
k = k - 1;
Mks{k} = (Rhat + Btil' * (Qtilks{k+1} + Kks{k+1}) * Btil)^-1;   % M0

%% calculate control policy
us       = zeros(size(Btil,2),Nstep);
xis      = zeros(size(xi0,1), Nstep);
xis(:,1) = xi0;

for j = 1:Nstep-1
    us(:,j) = - Mks{j} * ((0.5*(bNks{j+1}+hks{j+1})*Btil + xis(:,j)'*Atil'*(Qtilks{j+1}+Kks{j+1})*Btil)');%%%8
    xis(:,j+1) = Atil * xis(:,j) + Btil * us(:,j);
end
j = j + 1;
us(:,j) = -Mks{j}*((0.5*bQ_N'*Btil + xis(:,j)'*Atil'*Qfhat*Btil)');%30

%%matlab LQR 
usmatlab       = zeros(size(Btil,2),Nstep);
ximatlab      = zeros(size(xi0,1), Nstep);
ximatlab(:,1) = xi0;
[K,~,~]=dlqr(Atil,Btil,Q2,Rhat);
for j = 1:Nstep-1
    usmatlab(:,j) = - K* ximatlab(:,j);%%%8
    ximatlab(:,j+1) = Atil * ximatlab(:,j) + Btil * usmatlab(:,j);
end


%% evalute state of original system
%
 xt     = xauto+ep*real(P*xis);
%xt     =real(P*xis);
outdof = tcdim;
xout   = xt%(outdof,:);
% invariance PDE error
% res   = invariance_PDE_residual(DS,tspan,zt,ut);
% ratio = ratio_ext2int(DS,tspan,zt,ut);
%% output
traj      = struct();
traj.time = tspan;
traj.us   = us;
%traj.etat = etat;
%traj.xit  = xit;
traj.xt   = xout;
traj.xauto=xauto;
traj.mlus=usmatlab;
traj.mlxt=real(P*ximatlab);

end


