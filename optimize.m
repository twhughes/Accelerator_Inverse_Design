addpath(genpath('./'));                     % add the whole directory to path, if not already done

%% SET PARAMETERS
c0 = 1;                                     % speed of light m/s (normalized to 1)
lambda0 = 2;                                % central wavelength (um)

skip = 10;                                   % number of iteration frames between plots (higher->faster, lower->more plots)
display_plots = true;                       % plotting during the run?


alpha = 1e3;                                % step size in permittivity (~1e2-1e4 works well)
a = 10;                                    % smooth-max weight factor (see paper)
beta = 0.5;                                 % ratio of electron speed to speed of light
N = 1000;                                   % number of iterations

in_material = false;                        % evaluate E_max in material? or in surrounding regions
starting = 2;                               % 0 -> vacuum, 1 -> random, 2 -> midway epsilon

grids_in_lam = 100;                         % number of grid points in a free space wavelength
gap_nm       = 400;                         % gap size in nm
L = 3;                                      % size of optimization region (um)
% NOTE: if this ^ is too big and the epsilon is too large, the simulations
% can diverge.  This is because there are many degrees of freedom and
% resonance can occur very strongly. Need to try different values and see
% what works.
npml = 20;                                  % number of PML (absorbing region) points (need > 10 at least)

% relative permittivity of material region.  uncomment to select
%eps = 3.4363^2;     % Si 2um
eps = 1.4381^2;      % fused silica 2um
%eps = 1.9834^2;     % Si3N4
%eps = 1.9^2;        % GaOx

nmax = sqrt(eps);    % refractive index of material region

gamma = 1*0.999;                             % 'momentum term', see paper.  Set between 0-1, can speed up simulation in some cases

%% SET OTHER CONSTANTS (DON'T CHANGE)
dlx = lambda0/grids_in_lam;                 % grid size along electron trajectory axis
dly  = dlx;                                 % spacing in the perpendicular direction

pos_src = floor(npml+grids_in_lam/4);       % number of grid points between left edge and source
spc_pts = floor(grids_in_lam/4);            % number of grid points between source and structure
gap_pts = floor(gap_nm/1000/dlx);           % number of grid points in the gap
Lpts = round(L/dlx);                        % number of points in the optimization region

Nx = ceil(lambda0*beta/dlx);                % number of grid points in x
Ny = gap_pts+2*(pos_src + Lpts + spc_pts);  % number of grid points perpendicular to trajectory

nx = floor(Nx/2);
ny = floor(Ny/2);
    
% First compute G maximization, then do G/E_max maximization (for comparison)
for min_G_Emax = (1:1)

    alpha = alpha + 1e1*min_G_Emax;         % reset step size based on which optimization is being done (different objective functions)
    
    %% This section defines the input parameters that my FDFD code needs to run.
    %  see the FDFD.m code or FDFD_TFSF.m for a more detailed explanation.
    
    ER  = ones(Nx,Ny);                      % relative permittivity grid map
    MuR = ones(Nx,Ny);                      % relative permeability grid map
    ER_best = ones(Nx,Ny);                  % storing the best permittivity map
    A_best = 0;

    b = zeros(Nx,Ny);                       % TFSF map.  read up on total-field scattered-field if you are interested.
    b(:, pos_src:pos_src + spc_pts + Lpts + gap_pts + Lpts + spc_pts) = 1;  % define the total field region on the grid
    kinc = [0,1];                           % plane wave incident direction (perp. to electron)

    RES = [dlx,dly];                        % grid resolution vector
    BC = [-1,-1];                           % boundary condition vector (periodic if -1)
    NPML = [0,0,npml,npml];                 % PML cells on the boundaries (x-,x+,y-,y+)
    Pol= 'Hz';                              % field polarization
    spc = spc_pts*dly;                      % space between source and objects in um
    gap = gap_pts*dly;                      % gap size in um

    xs = dlx*(1:Nx);                        % constant to compute the eta object.  x-pos along gap.

    delta_device = zeros(Nx,Ny);            % delta_device is 0 where the permittivity doesn't change.  otherwise it is 1 in the optimization region.
    delta_device(1:Nx, pos_src + spc_pts : pos_src + spc_pts + Lpts) = 1;
    delta_device(1:Nx, pos_src + spc_pts + Lpts + gap_pts : pos_src + spc_pts + Lpts + gap_pts + Lpts) = 1;
    delta_device_vec = delta_device(:);     % vector version of delta_device (matlab likes this better)
    
    % define the eta vector field.  see the paper for more details.
    eta = zeros(Nx,Ny);  
    eta(:,ny) = 1/Nx*exp(2*pi*1i*dlx*(0:Nx-1)/lambda0/beta);       
    eta_vec = eta(:);

    % define stating permittivity based on what value the 'starting variable'
    % holds
    for i = (1:Nx)
        for j = (1:Ny)
            if (delta_device(i,j) == 1)
                if (starting == 1)
                    ER(i,j) = rand*(eps-1)+1;
                elseif (starting == 2)                
                    ER(i,j) = eps/2+0.5;
                else
                end
            end
        end
    end

    % run the simulation with accelerator input (plane wave) but all empty space
    [fields, ~] = FDFD_TFSF(ones(Nx,Ny),MuR,RES,NPML,BC,lambda0,Pol,b,kinc);

    % get the fields and the E0 (normalization)
    Ex = fields.Ex;
    Ey = fields.Ey;
    E0 = sqrt(abs(Ex(nx, ny))^2 + abs(Ey(nx, ny))^2);

    % define variables to store the iteration progress
    G_best = 0;                 % best gradient
    Gs = zeros(N,1);            % gradients over iteration
    E_maxs = zeros(N,1);        % max E-fields over iteration
    G_by_Es = zeros(N,1);       % G/E over iteration, computed directly
    G_by_Sa = zeros(N,1);       % G/E over iteration, computed with smooth-max

    phis = zeros(N,1);          % phase of the maximum accelerating input plane wave
    phi = 0;                    % assume input light phase of 0 to start
    AVM_prev = zeros(Nx,Ny);    % store previous sensitivity information for momentum update

    figure(1);                  % open a figure to plot
    
    for j = (1:N)

        % original simulation (structure in accelerator mode)
        [fields, extra] = FDFD_TFSF(ER,MuR,RES,NPML,BC,lambda0,Pol,b,kinc);
        % get fields
        Ex = fields.Ex/E0;
        Ey = fields.Ey/E0;
    
        % compute gradient (see paper)
        g = sum(sum(eta.*Ex));
        G = real(g);

        % get phase
        phis(j) = angle(g);

        % get numerical spatial derivative operators
        DEY = extra.derivatives.DEY;  
        DEX = extra.derivatives.DEX;    
    
        % turn the permittivity map into a vecor.  Then compute some
        % quantities for later.
        ER_vec = ER(:);
        delta_ER = ER > eps/2;
        delta_ER_vec = delta_ER(:);
        chi = delta_device.*(ER - ones(Nx,Ny));

        % compute operators from maxwell's eqs. turning Mz into Jx and Jy (Hz into Ex and Ey)
        Ox = -1i*lambda0/2/pi/c0*spdiags(1./ER_vec,0,Nx*Ny,Nx*Ny)*DEY;
        Oy =  1i*lambda0/2/pi/c0*spdiags(1./ER_vec,0,Nx*Ny,Nx*Ny)*DEX;

        % create the adjoint vector corresponding to eta (again, see paper)
        eta_aj = [eta_vec; zeros(Nx*Ny,1)];

        % if you are evaluating E_max in the material, compute |E| there,
        % otherwise, compute |E| in the full optimization region.
        if (in_material)
            E_abs = chi.*sqrt(abs(Ex).^2 + abs(Ey).^2);
        else
            E_abs = delta_device.*sqrt(abs(Ex).^2 + abs(Ey).^2);
        end        

        % Turn Ex and Ey into vectors.
        Ex_vec = reshape(delta_device.*Ex,[Nx*Ny,1]); 
        Ey_vec = reshape(delta_device.*Ey,[Nx*Ny,1]); 

        % compute auxiliary vectors for later.  (too complicated to explain
        % here.  ask me in person if you're interested).
        x_bar = E_abs(:);
        alpha_vec = exp(x_bar*a);
        alpha_T_1 = sum(alpha_vec);
        sigma_old = alpha_vec/alpha_T_1 + a*(alpha_vec.*x_bar)/alpha_T_1 - a*(sum(alpha.*x_bar)*alpha)/alpha_T_1/alpha_T_1;
        Sa = sum(alpha_vec.*x_bar)/alpha_T_1;

        %P = diag(1./x_bar);    
        %e = [Ex(:); Ey(:)];
        %R = P*Q*diag(conj(e));

        X_vec = conj(Ex(:))./x_bar;
        Y_vec = conj(Ey(:))./x_bar;
        X_vec(isinf(X_vec)) = 0;
        Y_vec(isinf(Y_vec)) = 0;

        sigma = [alpha_vec.*X_vec; alpha_vec.*Y_vec]/alpha_T_1 + a*([(alpha_vec.*x_bar).*X_vec;(alpha_vec.*x_bar).*Y_vec])/alpha_T_1 - a*(sum(alpha.*x_bar)*dlx*[alpha.*X_vec; alpha.*Y_vec])/alpha_T_1/alpha_T_1;

        %R = [diag(X_vec) diag(Y_vec)];
        %R(isnan(R)) = 0;
        %R_T = transpose(R);

        % compute adjoint terms for both gradient and E_max
        b_aj1 = real(g)/Sa^2*sigma;
        b_aj2 = -eta_aj/Sa;

        % construct final adjoint source
        if (min_G_Emax)
            b_aj = b_aj1 + b_aj2;
        else
            b_aj = -eta_aj;        
        end
        b_aj = reshape(Ox*b_aj(1:Nx*Ny) + Oy*b_aj(Nx*Ny+1:end),[Nx,Ny]);
        b_aj(isnan(b_aj)) = 0 ;
        
        % get the factored form of the system operator.
        AF = extra.AF;

        % run simulation with adjoint source.
        [fields_aj, ~] = FDFD_fast(ER,MuR,RES,NPML,BC,lambda0,Pol,b_aj,AF);
    
        % get fields
        x_aj = fields_aj.x/E0;
        Ex_aj = reshape(x_aj(1:Nx*Ny),[Nx,Ny]);
        Ey_aj = reshape(x_aj(Nx*Ny+1:end),[Nx,Ny]);

        % compute sensitivity information
        AVM = -real((Ex.*Ex_aj.*delta_device + Ey.*Ey_aj.*delta_device));    

        % record relevant variables in the arrays
        E_max = max(max((E_abs)));
        E_maxs(j) = E_max;
        Gs(j) = G;
        G_by_Es(j) = G/E_max;
        G_by_Sa(j) = G/Sa;

        % update permittivity
        ER = ER + alpha*AVM + alpha*gamma*AVM_prev;

        % update the previous sensitivity map
        AVM_prev = AVM;

        % if permittivity out of bounds, reset inside the correct bounds.
        ER(ER < 1) = 1;
        ER(ER > eps) = eps;        

        % record best permittivity if applicable
        if (G > G_best)
            ER_best = ER;
        end

        % plot stuff without too much hastle, display % done
        if (display_plots && mod(j,skip) == 0)
                perc_done = j/N*100
            clf;         
            subplot(2,2,1);
            disp = [];
            for k = (1:5)
                disp = [disp; real(ER)];
            end

            imagesc(disp,[1,eps])
            colormap(flipud(gray))    
            set(findall(gcf,'type','text'),'FontSize',22,'fontWeight','normal')
            set(gca,'FontSize',22,'fontWeight','normal')    
            colorbar()
            title('\epsilon');

            subplot(2,2,2);
            colorbar()
            plot(Gs(1:j),'k');
            xlabel('iteration number')
            ylabel('gradient (E_0)')
            set(findall(gcf,'type','text'),'FontSize',22,'fontWeight','normal')
            set(gca,'FontSize',22,'fontWeight','normal')

            subplot(2,2,3);
            plot((1:j),G_by_Es(1:j));
            hold all;
            plot((1:j),G_by_Sa(1:j));

            xlabel('iteration number');
            ylabel('G/E_max');
            %imagesc(real(Ex_aj*exp(1i*phi)),[-0.01,0.01]);
            set(findall(gcf,'type','text'),'FontSize',22,'fontWeight','normal')
            set(gca,'FontSize',22,'fontWeight','normal')  

            subplot(2,2,4); hold all;
            plot((1:j),phis(1:j));
            plot((1:j),zeros(j,1));
            xlabel('iteration number');
            ylabel('\phi');
            set(findall(gcf,'type','text'),'FontSize',22,'fontWeight','normal')
            set(gca,'FontSize',22,'fontWeight','normal')  
            pause(0.001);
        end
    end


    %% POST PROCESSING STUFF
    
    % create final field display
    ND = 5;
    field_disp = [];
    ER_disp = [];
    for i = (1:ND)
        field_disp = [field_disp;Ex];
        ER_disp = [ER_disp;  (ER-ones(Nx,Ny))*10000];
    end

    % plot movie
    NT = 0;         % number of time steps
    figure(2);
    for t = (1:NT)
        clf; 
        colormap(redbluecmap)
        imagesc(transpose(ER_disp + real(field_disp*exp(-1i*t/40))),[-5,5]); pause(0.0001);
    end

    % force the permittivity distribution binary
    eps_avg = (eps+1)/2;
    ER(ER<eps_avg) = 1;
    ER(ER>=eps_avg) = eps;
    
    % do another simulation of the binary distribution
    [fields, extra] = FDFD_TFSF(ER,MuR,RES,NPML,BC,lambda0,Pol,b,kinc);
    Ex = fields.Ex/E0;
    Ey = fields.Ey/E0;

    % compute the gradient
    g = sum(sum(eta.*Ex));
    G = abs(g);

    % get the maximum fields in material and optimization region
    E_abs = delta_device.*sqrt(abs(Ex).^2 + abs(Ey).^2);
    E_max = max(E_abs(:));

    E_abs_mat = (ER > 1).*sqrt(abs(Ex).^2 + abs(Ey).^2);
    E_max_mat = max(E_abs_mat(:));

    % compute the acceleration factors and save
    if (~min_G_Emax)
        ER_o = ER;    
        G_o = G;
        E_max_o = E_max;
        E_max_mat_o = E_max_mat;
        n_o = G_o/E_max_o;
        n_mat_o = G_o/E_max_mat;
    else    
        ER_p = ER;    
        G_p = G;
        E_max_p = E_max;
        E_max_mat_p = E_max_mat;
        n_p = G_p/E_max_p;
        n_mat_p = G_p/E_max_mat;
    end
end

% display the percent improvements (in optimization regions and in
% materials)
perc_improvement = (n_p-n_o)/n_o*100
perc_improvement_mat = (n_mat_p-n_mat_o)/n_mat_o*100
