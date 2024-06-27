function D_H=Compute_DH(Phi,Zeta,J,lambda1,mu1)
DW = 1; DH = 1; [nely,nelx] = size(Phi); Material_num = size(lambda1,2); 
%% Initial data: get volume fraction according to representative points' level set function value
nodenrs=reshape(1:nelx*nely,nely,nelx);        
edofVec=reshape(2*nodenrs(1:end-1,1:end-1)-1,(nelx-1)*(nely-1),1);   
edofMat=repmat(edofVec,1,8)+repmat([0 1 2*(nely-1)+[2 3 4 5] 2 3],(nelx-1)*(nely-1),1);
EleNodesID=edofMat(:,2:2:8)./2;
H = ones(nely,nelx);
for i = 1 : Material_num - 1
    Phi_max1 = Phi - Zeta(i);
    H(Phi_max1 >= 0) = i+1;
end
x_pre = sum( H(EleNodesID), 2 ) / 4;
x = reshape(x_pre,nely-1,nelx-1); x = round(x); x=flipud(x);
%% FEM data initialization 
nelx = nelx-1;  nely = nely-1;
nele = nelx*nely;
EW = DW / nelx;
EH = DH / nely;
%% Preparation FE analysis
nodenrs=reshape(1:(1+nelx)*(1+nely),1+nely,1+nelx);
edofVec=reshape(2*nodenrs(1:end-1,1:end-1)-1,nelx*nely,1);
edofMat=repmat(edofVec,1,8)+repmat([0 1 2*nely+[2 3 4 5] 2 3],nelx*nely,1);
iK=kron(edofMat,ones(8,1))';
jK=kron(edofMat,ones(1,8))';
lambda = zeros(size(x)); mu = lambda;  % Material properties in the different elements
for i = 1 : Material_num
    lambda = lambda + lambda1(i)*(x==i);
    mu = mu + mu1(i)*(x==i);
end
%% Stiffness, force matrix and necessary partial information of a unit cell
[keLambda, keMu, feLambda, feMu] = elementMats(EW/2,EH/2,inv(J));
%% Define loads and supports
alldofs = 1:2*(nely+1)*(nelx+1);
fixeddofs11 = 1:1:2*(nely+1);
fixeddofs12 = 2*nelx*(nely+1)+1:1:2*(nelx+1)*(nely+1);
fixeddofs21 = union( 1:2*(nely+1):2*nelx*(nely+1)+1 , 2:2*(nely+1):2*nelx*(nely+1)+2 );
fixeddofs22 = union( 2*(nely+1)-1:2*(nely+1):2*(nelx+1)*(nely+1)-1 , 2*(nely+1):2*(nely+1):2*(nelx+1)*(nely+1) );
freedofs = setdiff(alldofs,fixeddofs12);
freedofs = setdiff(freedofs,fixeddofs22);
freedofs = setdiff(freedofs,[1,2]);
%% Assemble the total stiffness matrix and the load vector
U = zeros(2*(nely+1)*(nelx+1),3);
sK = keLambda(:)*lambda(:).' + keMu(:)*mu(:).'; % The corresponding stiffness matrix entries
K  = sparse(iK(:), jK(:), sK(:)); K = (K+K')/2;
sF = feLambda(:)*lambda(:).'+feMu(:)*mu(:).'; % Assembly three load cases corresponding to the three strain cases
iF = repmat(edofMat',3,1);
jF = [ones(8,nele); 2*ones(8,nele); 3*ones(8,nele)];
F  = sparse(iF(:), jF(:), sF(:), 2*(nelx+1)*(nely+1), 3);
%% Apply periodic boundary conditions and solve
K(fixeddofs11,alldofs) = K(fixeddofs11,alldofs)+K(fixeddofs12,alldofs); 
K(alldofs,fixeddofs11) = K(alldofs,fixeddofs11)+K(alldofs,fixeddofs12);
F(fixeddofs11,:) = F(fixeddofs11,:)+F(fixeddofs12,:);
K(fixeddofs21,alldofs) = K(fixeddofs21,alldofs)+K(fixeddofs22,alldofs);
K(alldofs,fixeddofs21) = K(alldofs,fixeddofs21)+K(alldofs,fixeddofs22);
F(fixeddofs21,:) = F(fixeddofs21,:)+F(fixeddofs22,:);
U(freedofs,:) = K(freedofs,freedofs)\F(freedofs,:); % solve KU=F
U(fixeddofs12,:) = U(fixeddofs11,:); U(fixeddofs22,:) = U(fixeddofs21,:); 
%% Integration to calculate homogenized elastic moduli
feLambda1 = feLambda'; feMu1 = feMu';
temp_fe = feLambda1(:)*lambda(:)'+feMu1(:)*mu(:)';
temp_total=temp_fe(:);
D_equ=reshape(temp_total(:),3,8*nele)*U(edofMat',:);
De = repmat(lambda(:),1,9) .* repmat([1 1 0 1 1 0 0 0 0],nele,1) + repmat(mu(:),1,9) .* repmat([2 0 0 0 2 0 0 0 1],nele,1);
D_H = (reshape(sum(De),3,3)-D_equ/EW/EH*det(J))/nele;
disp('****** Homogenised elasticity tensor ******'); disp(D_H)
%% COMPUTE ELEMENT STIFFNESS MATRIX AND FORCE VECTOR
function [keLambda, keMu, feLambda, feMu] = elementMats(a,b,J1)
CMu = diag([2 2 1]); CLambda = zeros(3); CLambda(1:2,1:2) = 1; % Constitutive matrix contributions
poi_gau_x = [ -1/sqrt(3); 1/sqrt(3) ]; poi_gau_y = poi_gau_x; ceo_gau = [ 1; 1 ]; % Gauss point
keLambda = zeros(8,8); keMu = zeros(8,8);
feLambda = zeros(8,3); feMu = zeros(8,3);
L = [1 0 0 0;0 0 0 1;0 1 1 0];
for ii=1:length(poi_gau_x)
  for jj=1:length(poi_gau_y)
    JN = (1/4)*[-(1-poi_gau_y(jj)) (1-poi_gau_y(jj)) (1+poi_gau_y(jj)) -(1+poi_gau_y(jj));...
                -(1-poi_gau_x(ii)) -(1+poi_gau_x(ii)) (1+poi_gau_x(ii)) (1-poi_gau_x(ii))]; % Differentiated shape functions
    coo = J1 * [-a a a -a; -b -b b b]; J = JN * coo'; JJ = inv(J); % Jacobian
    weight = ceo_gau(ii)*ceo_gau(jj)*det(J); % Weight factor at this point
    G = kron(eye(2),JJ); dN = zeros(4,8); % Strain-displacement matrix
    dN(1,1:2:8) = JN(1,:); dN(2,1:2:8) = JN(2,:); dN(3,2:2:8) = JN(1,:); dN(4,2:2:8) = JN(2,:);
    B = L*G*dN;
    keLambda = keLambda + weight*(B' * CLambda * B); % Element matrices
    keMu = keMu + weight*(B' * CMu * B);
    feLambda = feLambda + weight*(B' * CLambda * diag([1 1 1])); % Element loads
    feMu = feMu + weight*(B' * CMu * diag([1 1 1]));
  end
end
%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%
%  Input:
%    Phi - Level set function value for the matrix cell, default size 101*101
%    Zeta - Level set height for the matrix cell
%    J - Jacobi matrix for the matrix cell deformation 
%    lambda1 - Lame's constant [0.3/(1-0.3^2)]
%    mu1 - Lame's constant
%  Output:
%    D_H - Homogenized elastic moduli
%  ---------------------------------------
%  Please feel free to contact us with any questions! 
%  - Chuang Ma, Dalian University of Technology
%  - Yichao Zhu, Dalian University of Technology
%  - Xu Guo, Dalian University of Technology
%  - chuangma@mail.dlut.edu.cn / yichaozhu@dlut.edu.cn / guoxu@dlut.edu.cn
%  [D_H]=Compute_DH(Phi_rectangle,[0.3 0.4],[1 0;0 1],[1e-10 0.329670329670330 0.659340659340659],[1e-10 0.384615384615385 0.769230769230769]);
%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%