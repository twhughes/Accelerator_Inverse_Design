solveropts.returnDiv = false;      %return an E field divergence operator for bound charge calculation below
inspect_only = false;
addpath(genpath('~/Documents/MATLAB/maxwellfdfd-master'));

N = 50;
ns = linspace(400,100000,N);
single_times=zeros(N,1);
double_times=zeros(N,1);
for i = (1:N)
    n = floor(sqrt(ns(i))/2);
    X = n;
    Y = n;
    x = floor(X/2);
    y = floor(Y/2);

    B = Box([-x,x; -y,y; 0,1]);
    tic;
    [E1, H1, obj_array, src_array, extra, epsTyler] = maxwell_run(...
        'OSC', 1, lam, ...
        'DOM', {'vacuum', 'none', 1.0}, [-X X; -Y Y; 0 1], [1,1,1] , BC.p, [0 0 0], ...
        'OBJ', ...
            {'box','k',eps}, ...
                B,...                            
        'SRCM', PointSrc(Axis.z,[0,0,0],1), ...
        solveropts, inspect_only);
    single_times(i) = toc();
    tic;
    [E2, H2, obj_array, src_array, extra, epsTyler] = maxwell_run(...
        'OSC', 1, lam, ...
        'DOM', {'vacuum', 'none', 1.0}, [-X X; -Y Y; 0 1], [1,1,1] , BC.p, [0 0 0], ...
        'OBJ', ...
            {'box','k',eps}, ...
                B,...                            
        'SRCM', PointSrc(Axis.z,[0,0,0],1), ...
        'SRCJ', PointSrc(Axis.z,[0.5,0.5,0.5],1), ...    
        solveropts, inspect_only);
    double_times(i) = toc();
%{
    Ex1 = E1{Axis.x}.data_original;
    Ey1 = E1{Axis.y}.data_original;
    Ez1 = E1{Axis.z}.data_original;
    Hx1 = H1{Axis.x}.data_original;
    Hy1 = H1{Axis.y}.data_original;
    Hz1 = H1{Axis.z}.data_original;

    Ex2 = E2{Axis.x}.data_original;
    Ey2 = E2{Axis.y}.data_original;
    Ez2 = E2{Axis.z}.data_original;
    Hx2 = H2{Axis.x}.data_original;
    Hy2 = H2{Axis.y}.data_original;
    Hz2 = H2{Axis.z}.data_original;
    %}
end
%%
clf; hold all;
scatter(ns./1000,single_times,'filled');
scatter(ns./1000,double_times,'filled');
xlabel('number of grid points (thousands)')
ylabel('simulation time (sec)')
legend('M_z source','M_z and J_z source');
set(findall(gcf,'type','text'),'FontSize',25,'fontWeight','normal')
set(gca,'FontSize',25,'fontWeight','normal')