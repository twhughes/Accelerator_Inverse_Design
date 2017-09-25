function [fields,extra] = FDFD(ER,MuR,RES,NPML,BC,lambda0,Pol,b,AF)
    
    extra = {};
    fields = {};
    k0 = 2*pi/lambda0;
    kinc = [1,0];
    dx = RES(1);
    dy = RES(2);

    NGRID = size(ER);
    NGRID2 = size(MuR);
    %assert(NGRID==NGRID2,'mu and epsilon matrices are not the same size')
    [Nx,Ny] = size(ER);
    %compute PML
    [sx,sy] = calcpml2d(NGRID,NPML);
    %
    ERxx = ER./sx.*sy;
    ERyy = ER.*sx./sy;
    ERzz = ER.*sx.*sy;
    MuRxx = MuR./sx.*sy;
    MuRyy = MuR.*sx./sy;
    MuRzz = MuR.*sx.*sy;

    %compute material matrices
    N = Nx*Ny;
    ERxx_inv = spdiags(1./ERxx(:),0,N,N);
    ERyy_inv = spdiags(1./ERyy(:),0,N,N);
    ERzz = spdiags(ERzz(:),0,N,N);
    MuRxx_inv = spdiags(1./MuRxx(:),0,N,N);
    MuRyy_inv = spdiags(1./MuRyy(:),0,N,N);
    MuRzz = spdiags(MuRzz(:),0,N,N);
    %compute derivative matrices
    [DEX,DEY,DHX,DHY] = yeeder(NGRID,k0*RES,BC,kinc/k0);
    %construct system matrix A
    
    if (nargin<9)
        if (strcmp(Pol,'Ez'))
            A = DHX*MuRyy_inv*DEX + DHY*MuRxx_inv*DEY + ERzz;
        elseif (strcmp(Pol,'Hz'))
            A = DEX*ERyy_inv*DHX + DEY*ERxx_inv*DHY + MuRzz;
        else
            exit(1)
        end
        AF = factorize(A);
    end
    b = b(:);
    f = AF\b;
    
    %f = mldivide(A,b);
   
    SX = spdiags(sx(:),0,Nx*Ny, Nx*Ny);
    SY = spdiags(sy(:),0,Nx*Ny, Nx*Ny);
    Hz = reshape(f,[Nx,Ny]);
    Ex = reshape(1i*ERxx_inv*DHY*f,[Nx,Ny]);
    Ey = reshape(-1i*ERyy_inv*DHX*f,[Nx,Ny]);
    Ez = Hz;
    Hx = Ex;
    Hy = Ey;
    
    x = [1i*ERxx_inv*DHY*f; -1i*ERyy_inv*DHX*f];
    
    fields.Ex = Ex;
    fields.Ey = Ey;
    fields.Ez = Ez;
    fields.Hx = Hx;
    fields.Hy = Hy;
    fields.Hz = Hz;
    fields.x = x;
    
    extra.eps = ER;
    extra.derivatives = {};    
    extra.derivatives.DEX = DEX;
    extra.derivatives.DEY = DEY;
    extra.derivatives.DHX = DHX;
    extra.derivatives.DHY = DHY;
    extra.AF = AF;
    
end