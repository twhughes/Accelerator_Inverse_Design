function [DEX,DEY,DHX,DHY] = yeeder(NGRID,RES,BC,kinc)

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Calculate yee grid derivative operators
        %
        % Note: for normalized grid, use this function as follows:
        % (DEX,DEY,DHX,DHY) = yeeder(NGRID,k0*RES,BC,kinc/k0);
        %
        % NGRID = [Nx,Ny]   grid dimensions
        % RES = [dx,dy]     grid resolution of 1x grid
        % BC = [xbc,ybc]    boundary conditions
        %       bc = 0     -> Dirichlet (noes not require kinc)
        %       bc = -1    -> Periodic  (requires kinc)
        %       bc = -2    -> Periodic  (requires kinc)
        % kinc = [kx, ky]   incident wave vector (only required for periodic BC
        %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        if nargin < 4
            kinc = [0,0];
        end
        
        if (BC(1) == -1)
           kinc = [0,0];
           BC(1) = -2;
        end
        if (BC(2) == -1)
           kinc = [0,0];
           BC(2) = -2;
        end
        
        Nx = NGRID(1);
        Ny = NGRID(2);
        
        dx = RES(1);
        dy = RES(2);
        
        xbc = BC(1);
        ybc = BC(2);

        kx = kinc(1);
        ky = kinc(2);

        N = Nx*Ny;

        Lambda_x = dx*Nx;
        Lambda_y = dy*Ny;

        %DEX = -1*speye(N,N);
        %DEY = -1*speye(N,N);

        % Construct DE operators

        if (Nx > 1)
        DEX = spdiags([-1*ones(N,1),ones(N,1)],[0,1],N,N);
        else
        DEX = spdiags([-1*zeros(N,1),ones(N,1)],[0,1],N,N);
        end
        if (Ny > 1)
        DEY = spdiags([-1*ones(N,1),ones(N,1)],[0,Nx],N,N);
        else
        DEY = spdiags([-1*zeros(N,1),ones(N,1)],[0,Nx],N,N);        
        end

        % fix DEX boundary conditions for Dirichlet boundary conditions
        for i = (Nx:Nx:N-Nx)
            DEX(i,i+1) = 0;
        end

        % if periodic boundary conditions in x

        if (xbc == -2)
            for i = (Nx:Nx:N)
               DEX(i,i-Nx+1) = exp(1i*Lambda_x*kx);
            end
        end

        if (ybc == -2)
           Z = spdiags(exp(1i*Lambda_y*ky)*ones(N,1),N-Nx,N,N);
           DEY = DEY+transpose(Z(:,1:N));
        end

        if (Nx == 1 && xbc == -2)
            DEX = (1i*kx).*speye(N);
        end
        if (Ny == 1 && ybc == -2)
            DEY = (1i*ky).*speye(N);
        end

        DEX = DEX./dx;
        DEY = DEY./dy;

        DHX = -DEX';
        DHY = -DEY';
end
