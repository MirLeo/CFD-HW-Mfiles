%% Stokes Eddies
%
% solves Navier Stokes equations for driven cavity problem using
% streamfunction/vorticity formulation
%System is solved using SBDF scheme
%   depends on:
%   ConstructRhsNS_IMEX
%   SystemMatNS_IMEX

    close all;  %clear; clc;
%
% Problem parameters:
    U = 1;
    Rmax = 1;
    alpha = 3*pi/2;
    Re = 1000;
    dt = 0.002;
    beta = 2*dt/(3*Re);
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
    PsiOmSys = SystemMatNS_IMEX(numUn, nP, nO, M, N, alpha, dr, dth, beta);
    PsiOmSys = sparse(PsiOmSys);
    A = PsiOmSys;
    [LL, UU, PP, QQ, RR] = lu(A);

     spparms('spumoni', 0)  % write out information about sparse algorithms
 %   spy(PsiOmSys)
 %   drawnow
    %rhs = ConstructRhk= 1:Tstepss(numUn, nP, nO, M, N, Rmax, dr, U);
    %time integration
    tic
    sol = zeros(numUn,1);
    solOld  = zeros(numUn,1);
    rhs = ConstructRhsNS_IMEX(numUn, nP, nO, M, N, Rmax, dr, dth, dt, U, sol, solOld);
%     Tfinal = 10;
%     Tsteps = Tfinal/dt;
    step = 1;
    err = 1;
    while abs(err)>1e-8
        psivort = QQ*(UU\(LL\(PP*(RR\rhs))));
        err = norm(sol -psivort)/norm(psivort)
        solOld = sol;
        sol = psivort;
        rhs = ConstructRhsNS_IMEX(numUn, nP, nO, M, N, Rmax, dr, dth, dt, U, sol, solOld);
        time = step * dt;
        if (mod(time,5)==0 && time<30)
            psi = reshape(psivort(1:numUn/2), size(rg));
            omega = reshape(psivort(numUn/2+1:numUn), size(rg));
            figure()
            subplot(1, 2, 1)
            contour(rg.*cos(thg), rg.*sin(thg), psi,40); colorbar;
            shading flat;  colormap(jet);  
            xlabel('x')
            ylabel('y')
            title(sprintf('Streamfunction, Re = %d, time = %d',Re, time),'fontsize',14)
            axis([-Rmax Rmax -Rmax Rmax])
            axis square
            subplot(1, 2, 2)
            contour(rg.*cos(thg), rg.*sin(thg), omega,80); colorbar;
            shading flat;  colormap(jet);  
            xlabel('x')
            ylabel('y')
            title('Vorticity','fontsize',14)
            axis([-Rmax Rmax -Rmax Rmax])
            axis square
        end
        step = step +1;
        error(step) = err;
       
    end
Tfinal = dt * step

% Plot
    psi = reshape(psivort(1:numUn/2), size(rg));
    omega = reshape(psivort(numUn/2+1:numUn), size(rg));
    figure()
    subplot(1, 2, 1)
        contour(rg.*cos(thg), rg.*sin(thg), psi,40); colorbar;
        shading flat;  colormap(jet);  
        xlabel('x')
        ylabel('y')
        title(sprintf('Streamfunction, Re = %d, Tfinal = %d',Re, Tfinal),'fontsize',14)
        axis([-Rmax Rmax -Rmax Rmax])
        axis square
    subplot(1, 2, 2)
        contour(rg.*cos(thg), rg.*sin(thg), omega,80); colorbar;
        shading flat;  colormap(jet);  
        xlabel('x')
        ylabel('y')
        title('Vorticity','fontsize',14)
        axis([-Rmax Rmax -Rmax Rmax])
        axis square
%%
%  Look for Eddies!
    figure()
    subplot(1, 2, 1)
        contour(rg.*cos(thg), rg.*sin(thg), psi, [0 0],'k','LineWidth',2); 
        hold on;
        contour(rg.*cos(thg), rg.*sin(thg), psi, 80); 
%         c = linspace(0, 4.4d-4, 40);
%         contour(rg.*cos(thg), rg.*sin(thg), psi, c); 
        shading flat;  colormap(jet);  
        xlabel('x')
        ylabel('y')
        title('Streamfunction')
        axis([0 Rmax 0 Rmax])
        axis square
    subplot(1, 2, 2)
        c = linspace(-20, 20, 80);
        contour(rg.*cos(thg), rg.*sin(thg), omega, c);
        shading flat;  colormap(jet);  
        xlabel('r')
        ylabel('\theta')
        title('Vorticity')
        axis([-Rmax Rmax -Rmax Rmax])
        axis square
