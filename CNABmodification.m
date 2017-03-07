%Modification in AdvectDiffuse.m   
case 5  %CNAB
            %
            %system matrix for remaining steps
            diag1 = (1+sigma)*ones(1,M);
            up1 = -sigma/2 * ones(1,M-1);
            low1 = -sigma/2 * ones(1,M-1);
            sysMatrix1 = diag(diag1) + diag(up1, 1) + diag(low1, -1);
            sysMatrix1(1, end) = -sigma/2;
            sysMatrix1(end, 1) = -sigma/2;
            %
            %system matrix for remaining steps
            diagonal = (1+sigma)*ones(1,M);
            up = -sigma/2 * ones(1,M-1);
            low = -sigma/2 * ones(1,M-1);
            sysMatrix = diag(diagonal) + diag(up, 1) + diag(low, -1);
            sysMatrix(1, end) = -sigma/2;
            sysMatrix(end, 1) = -sigma/2;