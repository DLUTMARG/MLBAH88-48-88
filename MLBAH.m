function Com = MLBAH(DW,DH,nelx,nely,J_0,variable_f)
%% Initial parameters for macroscopic design
num_ele = nelx*nely;  % mesh number along x and y axis
den(1:nely,1:nelx) = 1; % phi value of the elements (used to describe the shape of the design domain)
EW = DW / nelx;  % length of element
EH = DH / nely;    % width of element
[ x ,y ] = meshgrid( linspace(0,DW,nelx+1) , linspace(0,DH,nely+1) );% coordinates of nodes
%% Calculating derivatives of the mapping function (Jacobi matrices)
[coo_repx , coo_repy] = meshgrid(linspace(0,DW,nelx) , linspace(0,DH,nely));
coo_rep = [coo_repx(:)'+ EW/2 ; coo_repy(:)'+ EH/2];
J = [ J_0(1),J_0(5),J_0(6),J_0(11),J_0(12),J_0(12),J_0(13);
    J_0(3),J_0(8),J_0(9),J_0(15),J_0(16),J_0(16),J_0(17);
    J_0(2),J_0(6),J_0(7),J_0(12),J_0(13),J_0(13),J_0(14);
    J_0(4),J_0(9),J_0(10),J_0(16),J_0(17),J_0(17),J_0(18) ]*...
    [ ones(1,size(coo_rep(1,:),2));coo_rep(1,:);coo_rep(2,:);...
    coo_rep(1,:).^2;coo_rep(1,:).*coo_rep(2,:);coo_rep(1,:).*coo_rep(2,:);coo_rep(2,:).^2 ];
frac = variable_f(1).*coo_rep(1,:).*coo_rep(1,:) + variable_f(2).*coo_rep(2,:).*coo_rep(2,:) ...
    + variable_f(3).*coo_rep(1,:).*coo_rep(2,:) + variable_f(4).*coo_rep(1,:) + variable_f(5).*coo_rep(2,:) + variable_f(6);
J_norm = J./sqrt(J(1,:).*J(4,:)-J(2,:).*J(3,:));
% gradedstructure(J_0,Phi_max,phi,variable_f);
%% Preparation FE analysis
fixeddofs = 2*(nely+1)*nelx+1:2*(nely+1)*(nelx+1); % Define loads and supports
alldofs     = 1:2*(nely+1)*(nelx+1);
freedofs    = setdiff(alldofs,fixeddofs);
F = sparse(2*(nely+1)*(nelx+1),1);
loaddof = 3:2:2*(nely+1)-3;
F(loaddof) = 1/nely; F([1 2*(nely+1)-1]) = 0.5/nely;
%Preparation FE analysis
nodenrs=reshape(1:(1+nelx)*(1+nely),1+nely,1+nelx);        
edofVec=reshape(2*nodenrs(1:end-1,1:end-1)-1,nelx*nely,1);
edofMat=repmat(edofVec,1,8)+repmat([0 1 2*nely+[2 3 4 5] 2 3],nelx*nely,1);
EleNodesID=edofMat(:,2:2:8)./2;
iK=kron(edofMat,ones(8,1))'; 
jK=kron(edofMat,ones(1,8))';
%% Assembled stiffness matrix
KE = zeros(8,8); U = zeros(2*(nely+1)*(nelx+1),1); sK = zeros(64,num_ele);
[D_H] = ComputewithNN_DH(J_norm,frac);
for i=1:num_ele
    KE = BasicKe(D_H(:,i));  % homogenized stiffness matrix
    sK(:,i)=KE(:);
end
K = sparse(iK(:),jK(:),sK(:));
U(freedofs,:) = K(freedofs,freedofs)\F(freedofs,:);
Com = F'*U;
%% Post-processing and visualisation
Figt = 1; % 1 : x-direction displacement; 2 : y-direction displacement
p = patch('Faces',EleNodesID,'Vertices',[x(:),y(:)],'FaceColor','interp' ); 
cdata = U(Figt:2:end); 
cdata = cdata - (max(cdata)-min(cdata))/2; 
Tick = round(linspace(min(U(Figt:2:end)),max(U(Figt:2:end)),5)); hcb = colorbar; 
set(hcb,'Ticks',linspace(min(cdata),max(cdata),5),'FontSize',20,'TickLabels',{num2str(Tick(1)),num2str(Tick(2)),num2str(Tick(3)),num2str(Tick(4)),num2str(Tick(5))}); set(p,'FaceColor','interp','FaceVertexCData',cdata); set(p,'EdgeColor','none'); 
colormap(darkb2r(min(cdata),max(cdata))); axis equal; axis off;
end
%% -----------------------Subfunctions----------------------
%% Element stiffness
function Ke = BasicKe(D)
D = [D(1) D(2) D(3);D(2) D(4) D(5);D(3) D(5) D(6)]; a = 1 ; b = 1 ;
r11 =  [ b^2*D(1) , a*b*D(3) , a*b*D(7) , a^2*D(9) ]/16; k11 = [  16/3 ;  4 ;  4 ;  16/3 ]; k21 = [ -16/3 ;  4 ; -4 ;   8/3 ];
r12 =  [ b^2*D(7) , a*b*D(9) , a*b*D(4) , a^2*D(6) ]/16; k12 = [ -16/3 ; -4 ;  4 ;   8/3 ]; k22 = [  16/3 ; -4 ; -4 ;  16/3 ];
r21 =  [ b^2*D(3) , a*b*D(2) , a*b*D(9) , a^2*D(8) ]/16; k13 = [  -8/3 ; -4 ; -4 ;  -8/3 ]; k23 = [   8/3 ; -4 ;  4 ; -16/3 ];
r22 =  [ b^2*D(9) , a*b*D(8) , a*b*D(6) , a^2*D(5) ]/16; k14 = [   8/3 ;  4 ; -4 ; -16/3 ]; k24 = [  -8/3 ;  4 ;  4 ;  -8/3 ];
k31 = [  -8/3 ; -4 ; -4 ;  -8/3 ]; k41 = [   8/3 ; -4 ;  4 ; -16/3 ];
k32 = [   8/3 ;  4 ; -4 ; -16/3 ]; k42 = [  -8/3 ;  4 ;  4 ;  -8/3 ];
k33 = [  16/3 ;  4 ;  4 ;  16/3 ]; k43 = [ -16/3 ;  4 ; -4 ;   8/3 ];
k34 = [ -16/3 ; -4 ;  4 ;   8/3 ]; k44 = [  16/3 ; -4 ; -4 ;  16/3 ];
Ke1 = reshape([ r11 ; r21 ; r12 ; r22 ] * [ k11 , k12 , k13 , k14 ],2,8);
Ke2 = reshape([ r11 ; r21 ; r12 ; r22 ] * [ k21 , k22 , k23 , k24 ],2,8);
Ke3 = reshape([ r11 ; r21 ; r12 ; r22 ] * [ k31 , k32 , k33 , k34 ],2,8);
Ke4 = reshape([ r11 ; r21 ; r12 ; r22 ] * [ k41 , k42 , k43 , k44 ],2,8);
Ke = [ Ke1 ; Ke2 ; Ke3 ; Ke4 ];
Ke = Ke/(a*b);
end
function D_H = ComputewithNN_DH(J_norm,frac)
for i = 1:6
    load(['D:\Study_on_AH_and_opt\paper4\Publishing Code\NNtrain\VF_net\net',num2str(i),'\net']);
    load(['D:\Study_on_AH_and_opt\paper4\Publishing Code\NNtrain\VF_net\net',num2str(i),'\ps_i']);
    load(['D:\Study_on_AH_and_opt\paper4\Publishing Code\NNtrain\VF_net\net',num2str(i),'\ps_t']);  
    A1 = mapminmax('apply',[J_norm;frac],ps_i);
    anewn = sim(net,A1);
    anew(i,:) = mapminmax('reverse',anewn,ps_t);
end
D_H = zeros(6,size(J_norm,2));
for i=1:size(J_norm,2)
    R = [anew(1,i) anew(2,i) anew(4,i);0 anew(3,i) anew(5,i);0 0 anew(6,i)];
    D_H1 = R'*R;
    D_H(:,i) = [D_H1(1);D_H1(2);D_H1(3);D_H1(5);D_H1(6);D_H1(9)];
end
end
%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%
%  Input:
%    DW -  Width of design space
%    DH - Height of design space
%    nelx - Number of transverse elements 
%    nely - Number of vertical elements 
%    J_0 - Coefficients of the mapping function
%    variable_f - Coefficients of the level-set height function
%  Output:
%    Com - Compliance of the interesting structure
%  ---------------------------------------
%  Please feel free to contact us with any questions! 
%  - Chuang Ma, Dalian University of Technology
%  - Yichao Zhu, Dalian University of Technology
%  - Xu Guo, Dalian University of Technology
%  - chuangma@mail.dlut.edu.cn / yichaozhu@dlut.edu.cn / guoxu@dlut.edu.cn
%  MLBAH(2,1,200,100,[ 1; 0; 0; 1; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0 ],[ -0.06; -0.08; 0.07; -0.06; 0.19; 0.18 ]);
%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%