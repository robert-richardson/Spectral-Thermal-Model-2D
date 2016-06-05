function [A,B,C,E,rr,zz,Psi,Te,ye] = func_spectral_galerkin(p, N)
% This function implements the spectral-Galerkin method described in the
% paper. It's outputs consist of the state matrices for the resulting
% low-order model and additional variables storing 2-D distributed
% temperature data.

% Copyright (c) 2016 by Robert Richardson, Shi Zhao, David Howey
% and The Chancellor, Masters and Scholars of the University of Oxford.
% See the licence file LICENCE.txt for more information.

%% Construct basis functions (which satisfy the homogeneous BCs)
% Roughly corresponds to Section 2.1: Chebyshev-Galerkin approximation

% For convenience, define:
H = p.z2 - p.z1;                                        % cell height
R = p.r2 - p.r1;                                        % inner radius - outer radius

% Boundary conditions in required format (Eq. 3):
am = -p.hl/p.kr;    bm = 2/R;
ap = p.hr/p.kr;     bp = 2/R;
cm = -p.hb/p.kz;    dm = 2/H;
cp = p.ht/p.kz;     dp = 2/H;

% Grid points for numerical quadrature
Ng = 41;                                                % number of grid points (note: this is just for the numerical quadrature - it is different from N, the number of basis functions in each dimension)
Ng_centre = round(Ng/2);                                % grid point of centre node

% Evaluate radial/axial basis functions at grid points
% rhat/zhat =   radial/axial grid coordinates in domain [-1:1]
% Phi_r/Phi_z = radial/axial basis function values at grid points
[rhat, Phi_r] = func_basis_fn(Ng-1,am,bm,ap,bp,N-1);    % (Eq. 13)
[zhat, Phi_z] = func_basis_fn(Ng-1,cm,dm,cp,dp,N-1);    % (Eq. 16)

% Psi (column vector of basis functions, eq. 30)
N_s = N^2;                                              % number of model states
Psi = zeros(N_s,(Ng)*(Ng));                             % initialise Psi
k = 0;
for i = 1:1:N
    for j = 1:1:N
        Psi(k+1,:) = kron(Phi_r(i,:),Phi_z(j,:));       
        k = k+1;
        
    end
end


%% Construct state matrices
% Roughly corresponds to Section 2.3: Solution Algorthm

% Parameters for numerical quadrature
[~,w] = clencurt(Ng-1);                     % Clenshaw-curtis quadrature weights
wr = repelem(w,1,Ng)';                      % Clenshaw-curtis quadrature weights (radial)
wz = repmat(w,1,Ng)';                       % Clenshaw-curtis quadrature weights (axial)
[~,Dm] = chebdif(Ng,1);                     % Chebyshev differentiation matrix
D2m = Dm^2;                                 % Chebyshev differentiation matrix (2nd order)
Dr = 2/R*Dm;                                % scaling factor
D2r = (2/R)^2*D2m;                          % scaling factor 
D2z = (2/H)^2*D2m;                          % scaling factor
I = eye(size(D2m));                         % Identity matrix
r = p.r1+(rhat+1)*R/2;                      % radial coordinates (physical space)
z = p.z1+(zhat+1)*H/2;                      % axial coordinates (physical space)
rr = repelem(r,1,Ng)';                      % radial coordinate (kronecker form)
zz = repmat(z,1,Ng)';                       % axial coordinate (kronecker space) 

% Initialise matrices
E = zeros(N^2,N^2);
A = zeros(N^2,N^2); 
B = zeros(N^2,2);
C = zeros(4,N^2);

% Assign matrix values
for i = 0:N_s-1
    psi = Psi(i+1,:)';
    psir = (kron(Dr,I)*psi);
    psirr = (kron(D2r,I)*psi);
    psizz = (kron(I,D2z)*psi);
    for j = 0:N_s-1
        psj = Psi(j+1,:)';
        E(i+1,j+1) = p.rho*p.cp*sum(rr.*psi.*psj.*wr.*wz)*R/2*H/2;
        A(i+1,j+1) = sum((p.kr*rr.*psirr + p.kr*psir + p.kz*rr.*psizz).*psj.*wr.*wz)*R/2*H/2;
    end;
    B(i+1,1) = sum(rr.*psi.*wr.*wz)*R/2*H/2;
    psisq = reshape(psi,Ng,Ng);
    C(1,i+1) = psisq(end,Ng_centre);
    C(2,i+1) = psisq(1,Ng_centre);
    C(3,i+1) = psisq(Ng_centre,end);
    C(4,i+1) = psisq(Ng_centre,1);
end;


%% Homogenization of the boundary conditions
% Roughly corresponds to Section 2.2

% P matrices (Eq. 26)
P_t = zeros(N,N);
P_b = zeros(N,N);
P_r = zeros(N,N);
P_l = zeros(N,N);
for k = 1:N
    phi_rk = Phi_r(k,:);
    phi_zk = Phi_z(k,:);
    for i = 1:N
        phi_ri = Phi_r(i,:);
        phi_zi = Phi_z(i,:);
        P_b(k,i) = sum(r.*phi_rk.*phi_ri.*w)*R/2;
        P_t(k,i) = sum(r.*phi_rk.*phi_ri.*w)*R/2;
        P_l(k,i) = sum(p.r1*phi_zk.*phi_zi.*w)*H/2;
        P_r(k,i) = sum(p.r2*phi_zk.*phi_zi.*w)*H/2;
    end
end

% s vectors (Eq. 27)
s_b = zeros(N,1);
s_t = zeros(N,1);
s_l = zeros(N,1);
s_r = zeros(N,1);
for i = 1:N
    phi_ri = Phi_r(i,:);
    phi_zi = Phi_z(i,:);
    s_b(i,1) = sum(r.*phi_ri.*w)*R/2;
    s_t(i,1) = sum(r.*phi_ri.*w)*R/2;
    s_l(i,1) = sum(p.r1*phi_zi.*w)*H/2;
    s_r(i,1) = sum(p.r2*phi_zi.*w)*H/2;
end

% k/j values (Eq. 25)
k1 = cp + dp;       k2 = cp + 2*dp;
k3 = dm - cm;       k4 = cm - 2*dm;
j1 = ap + bp;       j2 = ap + 2*bp;
j3 = bm - am;       j4 = am - 2*bm;

% BCs (Eq. 3):
e_l = am*p.Tinfl;
e_r = ap*p.Tinfr;
e_b = cm*p.Tinfb;
e_t = cp*p.Tinft;

% d vectors (Eq. 24):
dI =    (k4/(k1*k4-k2*k3))*(P_b\(e_t*s_b) - (k2/k4)*(P_t\(e_b*s_t)));
dII =   (k1/(k1*k4-k2*k3))*(P_t\(e_b*s_t) - (k3/k1)*(P_b\(e_t*s_b)));
dIII =  (j4/(j1*j4-j2*j3))*(P_l\(e_r*s_l) - (j2/j4)*(P_r\(e_l*s_r)));
dIV =   (j1/(j1*j4-j2*j3))*(P_r\(e_l*s_r) - (j3/j1)*(P_l\(e_r*s_l)));

% Auxilary function (Eq. 19)
Te = kron(Phi_r'*dI,zhat') + kron(Phi_r'*dII,(zhat').^2) + ...
    kron(rhat',Phi_z'*dIII) + kron((rhat').^2,Phi_z'*dIV);

T_esq = reshape(Te,Ng,Ng);
ye(1) = T_esq(end,Ng_centre);
ye(2) = T_esq(1,Ng_centre);
ye(3) = T_esq(Ng_centre,end);
ye(4) = T_esq(Ng_centre,1);
% ye(1) = ue(end);  ye(2) = ue(Ng); ye(3) = ue(end-Ng+1); ye(4) = ue(1);

% uerr and uezz
Ter = kron(Dr,I)*Te;
Terr = kron(D2r,I)*Te;
Tezz = kron(I,D2z)*Te;

for i = 0:N_s - 1;
    psi = Psi(i+1,:)';
    B(i+1,2) = sum((p.kr*rr.*Terr + p.kz*rr.*Tezz + p.kr.*Ter).*psi.*wr.*wz)*R/2*H/2;
end


end




