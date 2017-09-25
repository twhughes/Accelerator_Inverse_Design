function [sx,sy] = calcpml2d(NGRID,NPML)

        eta0 = 376.73;     % free space impedence
        a_max = 3.0;       % 0 <= a_max <= 5
        p = 3.0;           % 3 <= p <= 5
        sigma_p_max = 1.0; % sigma_prime ~ 1

        Nx = NGRID(1);        Ny = NGRID(2);

        Nxlo = NPML(1); Nxhi = NPML(2); Nylo = NPML(3); Nyhi = NPML(4);


        sx = ones(Nx,Ny);
        sy = ones(Nx,Ny);

        for nx = (1:Nxlo)
                xp = nx/Nxlo;
                sig_p_x = sigma_p_max*(sin(pi*xp/2)^2);
                a = 1 + a_max*(xp)^p;
                sx(Nxlo-nx+1,:) = a*(1+1i*eta0*sig_p_x);
        end
        for nx = (1:Nxhi)
                xp = nx/Nxhi;
                sig_p_x = sigma_p_max*(sin(pi*xp/2)^2);
                a = 1 + a_max*(xp)^p;
                sx(Nx-Nxhi+nx,:) = a*(1+1i*eta0*sig_p_x);
        end
        for ny = (1:Nylo)
                yp = ny/Nylo;
                sig_p_y = sigma_p_max*(sin(pi*yp/2)^2);
                a = 1 + a_max*(yp)^p;
                sy(:,Nylo-ny+1) = a*(1+1i*eta0*sig_p_y);
        end      
        for ny = (1:Nyhi)
                yp = ny/Nyhi;
                sig_p_y = sigma_p_max*(sin(pi*yp/2)^2);
                a = 1 + a_max*(yp)^p;
                sy(:,Ny-Nyhi+ny) = a*(1+1i*eta0*sig_p_y);
        end
end
