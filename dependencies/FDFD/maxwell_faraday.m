addpath('~/Documents/Fan/Code/Accelerator/AVM/');
addpath('~/Documents/MATLAB/Factorize/Factorize/');
c0 = 1;
ER = ones(200,468);
N = 1;
epss = linspace(1,20,N);
de = 0.0001;
[Nx,Ny] = size(ER);
nx = floor(Nx/2);
ny = floor(Ny/2);
MuR = ones(Nx,Ny);
RES = [0.01,0.01];
NPML = [0,0,10,10];
BC = [-1,-1];
lambda0 = 2;
L = 1;
beta = 1;
dl = 0.01;
src = 0.5*L;
spc = L/4;
Pol='Hz';
b = zeros(Nx,Ny);
b(:,40) = 1;
b(:,end-40) = -1;
Gs = zeros(N,1);
AVM_s = zeros(N,1);
AVM_s2 = zeros(N,1);
D_s = zeros(N,1);
AVM_structmap = zeros(Nx,Ny,N);
for i = 1:N
    perc_done = i/N*100
    %original
    epsilon = epss(i);
    ER(floor(nx-L/dl/2):floor(nx+L/dl/2),floor(ny-spc/dl-L/dl):floor(ny-spc/dl)) = epsilon;
    ER(floor(nx-L/dl/2):floor(nx+L/dl/2),floor(ny+spc/dl):floor(ny+spc/dl+L/dl)) = epsilon;    
    tic; [fields, extra] = FDFD(ER,MuR,RES,NPML,BC,lambda0,Pol,b); toc;
    Ex = fields.Ex;
    [G,~] = calc_grad_two_squares(Ex,beta,dl,2*pi/lambda0);
    Gs(i) = G;    
    DEY = extra.derivatives.DEY;
    %adjoint
    O = -1i*lambda0/2/pi/c0*spdiags(1./ER(:),0,Nx*Ny,Nx*Ny)*DEY;
    eta = zeros(Nx,Ny);    
    js = (0:Nx-1);
    eta(:,ny) = exp(2*pi*1i*dl*js/lambda0/beta);       
    eta = eta(:);
    b_aj = -O*eta;
    tic; [fields_aj, extra_aj] = FDFD(ER,MuR,RES,NPML,BC,lambda0,Pol,b_aj); toc;
    z = fields_aj.Hz(:); 
    y = fields.Hz(:);
    Ey = fields.Ex;
    Ex_aj = fields_aj.Ex;
    Ey_aj = fields_aj.Ey;    
    g = transpose(eta)*O*y;
    delta_device = (ER - ones(Nx,Ny))/(epsilon-1);
    d = delta_device(:);
    d = spdiags(d,0,Nx*Ny,Nx*Ny);
    DEX = extra.derivatives.DEX;
    DEY = extra.derivatives.DEY;
    AVM_s(i) = 1/lambda0/beta./abs(g).*(conj(g).*(2*pi/lambda0)^2*dl*dl*(sum(sum(delta_device.*(Ex.*Ex_aj + Ey.*Ey_aj)))));
    AVM_s2(i) = 1/lambda0/beta./abs(g).*(conj(g).*(2*pi/lambda0)^2*dl*dl*1/(epsilon^2)*(transpose(y)*(DEX*d*DEX+DEY*d*DEY)*z));
    %{ 
    zi = 0;
    for xi = (1:Nx)
        for yi = (1:Ny)
            zi = zi+1
            zp = zeros(Nx*Ny,1);
            zp(zi) = y(zi);
            yp = zeros(Nx*Ny,1);
            yp(zi) = y(zi);            
            AVM_structmap(xi,yi,i) = 1/lambda0/beta./abs(g).*(conj(g).*(2*pi/lambda0)^2*dl*dl*1/(epsilon^2)*(transpose(yp)*(DEX*d*DEX+DEY*d*DEY)*zp));
        end
    end
    %}
    %Direct
    eps_direct = [epsilon-de,epsilon+de];
    d_G = zeros(2,1);
    for j = (1:2)
       epsilon = eps_direct(j);
        ER(floor(nx-L/dl/2):floor(nx+L/dl/2),floor(ny-spc/dl-L/dl):floor(ny-spc/dl)) = epsilon;
        ER(floor(nx-L/dl/2):floor(nx+L/dl/2),floor(ny+spc/dl):floor(ny+spc/dl+L/dl)) = epsilon;    
        tic; [fields, extra] = FDFD(ER,MuR,RES,NPML,BC,lambda0,Pol,b); toc; 
        Ex = fields.Ex;
        [G,~] = calc_grad_two_squares(Ex,beta,dl,2*pi/lambda0);
        d_G(j) = G;
    end
    D_s(i) = (d_G(2)-d_G(1))/de;
end

%%
figure(1); clf; hold all;
plot(epss,Gs);
xlabel('structure \epsilon')
ylabel('Gradient')
set(findall(gcf,'type','text'),'FontSize',22,'fontWeight','normal');
set(gca,'FontSize',22,'fontWeight','normal');

figure(2); clf; hold all;
plot(epss,real(AVM_s2)*65);
plot(epss,real(D_s));
legend('AVM','Direct')
xlabel('structure \epsilon')
ylabel('Gradient Sensitivity')
set(findall(gcf,'type','text'),'FontSize',22,'fontWeight','normal');
set(gca,'FontSize',22,'fontWeight','normal');