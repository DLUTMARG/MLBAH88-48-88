function []=gradedstructure(J_0,Phi_max,fun)
%% J of each cell 
% 11, 21, 12, 22 ***(1,2,3,4)**** := 1 , x , y , x^2 , xy , xy , y^2
J = [ J_0(1),J_0(5),J_0(6),J_0(11),J_0(12),J_0(12),J_0(13);
    J_0(3),J_0(8),J_0(9),J_0(15),J_0(16),J_0(16),J_0(17);
    J_0(2),J_0(6),J_0(7),J_0(12),J_0(13),J_0(13),J_0(14);
    J_0(4),J_0(9),J_0(10),J_0(16),J_0(17),J_0(17),J_0(18) ];

%% draw picture: phi( Jx /h)
alf = 0.0;
DW = 2; DH = 1;
nelx = 1600; nely = 800;
EW = DW / nelx;  % length of element
EH = DH / nely;  % width of element
[ x ,y ] = meshgrid( EW * ( 0 :  nelx) , EH * ( 0 : nely) );

enlarge_num=0.1; clf;
xx = x(:)'; yy = y(:)';

pit_y = [ J(1), J(3); J(2), J(4) ]*[xx;yy] + ...
    [ J(5), J(7); J(6), J(8) ]*[xx.*xx;xx.*yy]/2 + ...
    [ J(9), J(11); J(10), J(12) ]*[yy.*xx;yy.*yy]/2 + ...
    [ J(13), J(15); J(14), J(16) ]*[xx.*xx.*xx;xx.*xx.*yy]/3 + ...
    [ J(17), J(19); J(18), J(20) ]*[xx.*yy.*xx;xx.*yy.*yy]/3 + ...
    [ J(21), J(23); J(22), J(24) ]*[xx.*yy.*xx;xx.*yy.*yy]/3 + ...
    [ J(25), J(27); J(26), J(28) ]*[yy.*yy.*xx;yy.*yy.*yy]/3;

pit_y = pit_y/enlarge_num;
coo = pit_y-floor(pit_y);
[ xx_i ,yy_i ] = meshgrid( 0.01 * ( 0 :  100) , 0.01 * ( 0 : 100) );
Phi_Jx_de = interp2(xx_i,yy_i,reshape(Phi_max,101,101),reshape(coo(1,:),nely+1,nelx+1),reshape(coo(2,:),nely+1,nelx+1));
% Equipotential line
frac = (fun(1)*x.*x/DW/DW + fun(2)*y.*y/DH/DH + fun(3)*x.*y/DW/DH + ...
    fun(4)*x/DW + fun(5)*y/DH + fun(6));
% regular
frac(frac<0) = 0;
frac(frac>0.33) = 0.33;
hold on;
contourf(x,  y, Phi_Jx_de-frac, [alf,alf]);
axis equal
axis off
end