% Crank Nicolson Adamas Bashfort scheme stability region

beta = 0:0.01:2;
alpha = 0:0.01:10;
alpha = -alpha;
[X,Y] = meshgrid(alpha,beta);

Z1 = ((-X/2 - (3*1i*Y)/2 - 1)+sqrt((X/2 + (3*1i*Y)/2 + 1).^2 -...
     1i*(3 - X).*Y))./(3 - X) ;
 
Z2 = ((-X/2 - (3*1i*Y)/2 - 1)-sqrt((X/2 + (3*1i*Y)/2 + 1).^2 -...
     1i*(3 - X).*Y))./(3 - X) ;
 
 Z1 = abs(Z1);
 Z2 = abs(Z2);
 Z = zeros(size(Z1));
 Z(1:end) = max(Z1(1:end),Z2(1:end));
 V = 0:0.1:1;
 contour(X,Y,Z,V,'linewidth',3)
 set(gca, 'fontsize',14)
 xlabel('\alpha')
 ylabel('\beta')
 title('Stability Contour for CNAB')