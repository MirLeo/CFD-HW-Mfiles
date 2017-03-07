%% Stokes Eddies
%
% solves steady Navier Stokes equations for driven cavity problem using
% streamfunction/vorticity formulation
%
%   depends on: 
%   SystemMat
%   ConstructRhsSteadyNS
    close all;  clear; clc;
%
% Problem parameters:
    U = -1;
    Rmax = 1;
    alpha = pi/2;
    Re = -70;
%   
% Set up finite difference grid
    M =100; dr = Rmax/(M-1);
    N =100; dth = alpha/(N-1);
    [rg, thg] = meshgrid(0: dr :Rmax, ...
                         alpha: -dth: 0);                     
%
% Unknowns and numbering
    numUn = M*N;
    nP = reshape(1:numUn, size(rg));
    nO = reshape(numUn+1:2*numUn, size(rg));
    numUn = 2*numUn;
%
% Build system matrices and rhs
    PsiOmSys = SystemMat(numUn, nP, nO, M, N, alpha, dr, dth);
    PsiOmSys = sparse(PsiOmSys);
    A = PsiOmSys;
    [LL, UU, PP, QQ, RR] = lu(A);

    spparms('spumoni', 0)  % write out information about sparse algorithms
    spy(PsiOmSys)
    drawnow
    %rhs = ConstructRhs(numUn, nP, nO, M, N, Rmax, dr, U);
%%
%
% Solve
%     tic
%     psivort = A \ rhs ;  % sparse solve
%     t = toc;
%     disp(['Time taken for linear system solve = ', num2str(t)]);
%     
    %iterative solve
    sol0 = zeros(numUn,1);
    rhs = ConstructRhsSteadyNS(numUn, nP, nO, M, N, Rmax, dr, dth, U, sol0, Re);
    iter = 1;
    iterMax = 1000;
    Tol = 1e-10;
    err = 1;
    
    while abs(err) > Tol
        
        psivort = QQ*(UU\(LL\(PP*(RR\rhs))));
        rhs = ConstructRhsSteadyNS(numUn, nP, nO, M, N, Rmax, dr, dth, U, psivort, Re);
        iter = iter +1;
        err = norm(sol0 -psivort)/norm(psivort)
        sol0 = psivort;
        error(iter)=err;
        if iter == iterMax
            break
        end
    end
    iter
    figure()
    plot(log(error))
    title('error decay, Re=70')
    ylabel('log(error)')
    xlabel('iterations')
    %
% Plot
    psi = reshape(psivort(1:numUn/2), size(rg));
    omega = reshape(psivort(numUn/2+1:numUn), size(rg));
    figure()
    subplot(1, 2, 1)
        pcolor(rg.*cos(thg), rg.*sin(thg), psi); colorbar;
        shading flat;  colormap(jet);  
        xlabel('x')
        ylabel('y')
        title('Streamfunction, Re=40')
        axis([0 Rmax 0 Rmax])
        axis square
    subplot(1, 2, 2)
        pcolor(rg.*cos(thg), rg.*sin(thg), omega); colorbar;
        shading flat;  colormap(jet);  
        xlabel('x')
        ylabel('y')
        title('Vorticity, Re=1')
        axis([0 Rmax 0 Rmax])
        axis square
%%
%  Look for Eddies!
    figure()
    subplot(1, 2, 1)
        contour(rg.*cos(thg), rg.*sin(thg), psi, [0 0],'k','LineWidth',2); 
        hold on;
        contour(rg.*cos(thg), rg.*sin(thg), psi, 40); 
%         c = linspace(0, 4.4d-4, 40);
%         contour(rg.*cos(thg), rg.*sin(thg), psi, c); 
        shading flat;  colormap(jet);  
        xlabel('x')
        ylabel('y')
        title('Streamfunction')
        axis([0 Rmax 0 Rmax])
        axis square
    subplot(1, 2, 2)
        c = linspace(-20, 20, 40);
        contour(rg.*cos(thg), rg.*sin(thg), omega, c);
        shading flat;  colormap(jet);  
        xlabel('r')
        ylabel('\theta')
        title('Vorticity')
        axis([0 Rmax 0 Rmax])
        axis square
