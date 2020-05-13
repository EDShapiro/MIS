%function(  grid , kldiv) =  surconstruct(dim, level, type grid? Implement options to choose from on the type of grid)
dim = 2;
level = 4;
grid = sgpp.Grid.createPolyGrid(dim,8);
grid.getGenerator().regular(level);
gridStorage = grid.getStorage();

%%
%Retrieveing Linear CLenshaw Curtis Grid Points
y = zeros(gridStorage.getSize(),dim);
    for i = 0:gridStorage.getSize()-1
           gp = gridStorage.getPoint(i);
           y(i+1,:) = [gp.getStandardCoordinate(0), gp.getStandardCoordinate(1)];
    end
    gridsize(level) =length(y);
%%
% Sparse grid, grid point conversion for sparse grid surrogate construction
% of Vlasov-Poisson Nonlinear Landau Damping Problem 
% %Convert [0,1] to [a,b] for each random variable (f^(-1) = x*(b-a) +a)
% lower1 = -1; upper1 = 1;
% lower2= 0.0 ; upper2 = pi/2;
% xconv = zeros(length(y),2);
% xconv(:,1) = y(:,1)*(upper1 - lower1) + lower1; %Scaling sparse grid points to be from [lower, upper]
% xconv(:,2) = y(:,2)*(upper2 - lower2) + lower2;
% eps1d = xconv(:,1);
% eps2d = xconv(:,2);
% samples = length(xconv(:,1));
% ynew = cat(2, eps1d, eps2d);
%%
%Modified script to incorporate uncertainty in background field and
%instability
%eps=0.5; % Amplitude of perturbation, 0.05 for linear, 0.5 for nonlinear
kx=0.5;    % Wave vector
L=2*pi/kx; % length of domain
qm=-1;    % negative charge to mass = electrons, $\frac{q}{m}=-1$

% Initial condition
f0 =@(x,v,eps1d, eps2d)   (1+eps1d.*cos((kx)*x))./sqrt(2*pi).*exp(-0.5.*(v + eps2d).^2);
 
% Background to be subtracted from the phase space plot
background=@(x,v, mu) 1./sqrt(2*pi).*exp(-0.5.*(v-mu).^2);
% External electric field
E_ext=@(x,t) 0.*x+0.*t;

dt=0.1; % Time step
tmax=50;
rungekutta_order=3; %Order of the runge Kutta time integrator


Nx=32; % Number of cells in spatial direction
Nv=32; % Number of cells in velocity direciton
vmax=4.5;
vmin=-4.5;
%%
v=linspace(vmin,vmax,Nv+1).';v=v(1:end-1);
x=linspace(0,L,Nx+1).';x=x(1:end-1);

[XX,VV]=ndgrid(x,v); 
%%
switch rungekutta_order
    case 1
        rksd=[1, 0]; %Symplectic Euler
        rksc=[0, 1];
    case 2
        rksd=[0.5, 0.5 ];
        rksc=[0, 1];
    case 3
        rksd=[2/3, -2/3, 1  ];
        rksc=[ 7/24, 3/4, -1/24];
    case 4
        rk4sx=real((2^(1/3) +2^(-1/3)-1)/6);
        rksd=[ 2*rk4sx +1 , -4*rk4sx-1, 2*rk4sx+1, 0];
        rksc=[ rk4sx + 0.5 , -rk4sx, -rk4sx, rk4sx +0.5];
end

%Build Matrix for spatial advection and all stages
XSHIFT=zeros(Nx,Nv,length(rksd));
m=fftshift(((1:Nx)-Nx/2-1));
for jdx=1:length(rksd)
   XSHIFT(:,:,jdx)=exp((-1j*dt*rksd(jdx)*2*pi/L)*(v*m).');
end
%Preallocate matrix for velocity advection
VSHIFT=zeros(Nx,Nv);

%diagnostic
numt=ceil(tmax/dt);
time = linspace(0,tmax,numt);
fieldenergy=zeros(numt,1);
kineticenergyi=zeros(numt,gridStorage.getSize());

EE=zeros(Nx,numt);


% Set initial condition
%%
%Defining a 3-D Array the initial particle distributions based off of the
%perturbed data.


for i = 1:gridStorage.getSize()
   
fi(:,:,i)=f0(XX,VV,y(i,1), y(i,2));
fi(:,:,i)=fft(fi(:,:,i),[],1);
 
end

% for i = 1:samples
%     b(:,:,i) = background(
% end

%%
%Unit Test for UQ with Sparse Grid Points
%figure('Name', 'Phase Space Density', 'Numbertitle','off');
%Adding loop to consider perturbations in the electric field
%Unit Test for UQ
%figure('Name', 'Phase Space Density', 'Numbertitle','off');
%Adding loop to consider perturbations in the electric field
for i =1:gridStorage.getSize()
    
    
    for tdx=1:numt
    for sdx=1:length(rksd)
        rho=sum(fi(:,:,i),2).*(vmax-vmin)/Nv; %rho(x) -> integrate over v
        rho(1)=0; %remove (ion) background
        E=rho./(-1j*(2*pi/L)*fftshift(((1:Nx)-Nx/2-1).'));
        % remove constant fourier mode
        E(1)=0;

        if (sdx==length(rksd))
            %fieldenergy(tdx)=sum(permute(E,[2 1])*E)/2; %L2 norm of Electric field
            %EE(:,tdx)=E;
        end

        E=ifft(E,'symmetric');
        E=E+E_ext(x,tdx*dt);
        %plot(x,ifft(E,'symmetric'));
            fi(:,:,i)=ifft(fi(:,:,i),[],1,'symmetric'); %back transform
        % f is fully backtransformed, calculate diagnostics
        if (sdx==length(rksd))
            
                kineticenergyi(tdx,i)=sum(fi(:,:,i),1)*v.^2*(vmax-vmin)/Nv;
                %sh_entropy(tdx)=sum(sum(f1.*log(abs(f1))))*L/Nx*(vmax-vmin)/Nv;
            %kl_entropy(tdx)=sum(sum(f1.*log(abs(f1./f0(XX,VV)))))...
                %*L/Nx*(vmax-vmin)/Nv;
            
        end

        % Fourier transform in v
        fi(:,:,i)=fft(fi(:,:,i),[],2);
        % Build Matrix for spatial advection
         m=fftshift((1:Nv)  -Nv/2-1);
         VSHIFT=exp((1j*dt*rksc(sdx)*2*pi/(vmax-vmin))*E.*m);

        fi(:,:,i)=fi(:,:,i).*VSHIFT;
        fi(:,:,i)=ifft(fi(:,:,i),[],2,'symmetric');
        %Advection in X
        fi(:,:,i)=fft(fi(:,:,i),[],1);
        fi(:,:,i)=fi(:,:,i).*XSHIFT(:,:,sdx);
    end
        
    end
    %Plot 1 kinetic energy sample from a given pertubration

end


%% Finding the time centered average around the 400th time point of the simulation.

keqoi = kineticenergyi(390:410,:); %Setting up an array for the quantity of interest at a time point (t=40),
% when the system is in a steady state, across all sparse grid points.
mkeqoi = zeros(gridStorage.getSize(),1);
for i = 1: gridStorage.getSize()
    mkeqoi(i) = mean(keqoi(:,i)); 
end
%%
% Creating a coefficients vector for the interpolant - Running more than
% once with different data causes the script to crash
alpha = sgpp.DataVector(gridStorage.getSize());
alpha.setAll(0);
for i = 0:gridStorage.getSize()-1
    gp = gridStorage.getPoint(i);
    alpha.set(i, mkeqoi(i+1));
end
%Hierarchisation of coefficients of the polynomial interpolant - I am not
%sure what is going on here.
sgpp.createOperationHierarchisation(grid).doHierarchisation(alpha);


%%
% Mapping random sample points to grid values, and then creating a 2-D grid
% of input values for the interpolant.
load('RandLD1.mat', 'epsm', 'KE769') %Loading the necessary data from brute force simulation
%% Turn this into a function or program - Find the KDE for one
%Finding the mean kinetic energy for every perturbation in first dim, fixed 1
%point for second dim.
mkeqoi1 =zeros(length(epsm49(:,1)),1);
for i = 1:length(epsm49(:,1))
    mkeqoi1(i) = mean(KE49(390:410,i));
end





%% Test to ensure interpolant is setup properly by evaluating the interpolant at the initial points. Unit Test
for i = 1:gridStorage.getSize()
    p = sgpp.DataVector(dim);
    p.set(0, y(i,1));
    p.set(1, y(i,2));
    opEval = sgpp.createOperationEvalNaive(grid);
    a = opEval.eval(alpha,p)
if a ~= mkeqoi
    fprintf('Function not interpolated properly')
    break 
end
end
%% MSE Calculation and Vector Setup for KDE
u1 = zeros(length(epsm49(:,1)),1);

%Evaluating interpolant at a set of perturbation points from the brute
%force simulation.
for i = 1:length(epsm49(:,1))
p = sgpp.DataVector(dim);
p.set(0, epsm49(i,1));
p.set(1, epsm49(i,2));
opEval = sgpp.createOperationEvalNaive(grid);
u1(i) = opEval.eval(alpha,p);
clear p opeval 
end
%%
save('dist.mat', 'u1', 'mkeqoi1')
%%
%Calculating interpolant error by evaluating MSE between sparse grid and
% brute force results.
for i = 1:length(epsm49(:,1))
    e(i) = (u1(i) - mkeqoi1(i))^2;%
end
    ef1 = sqrt(sum(e))/length(epsm49(:,1))
%% Creating KDEs
pd6 = fitdist(mkeqoi1, 'Kernel')
x6 = min(mkeqoi1)-2:1:max(mkeqoi1) + 2;
p6 = pdf(pd6,x6);
pds6 = fitdist(u1, 'Kernel')
xs6 = min(u1)-2:1:max(u1) + 2;
ps6 = pdf(pds6,xs6);
%%
%% Kullback-Liebler Divergence Calculation 
%Choosing the min/max values of the data to setup a region of integration
if min(u1) <= min(mkeqoi1)
    minv = min(u1);
else
    minv = min(mkeqoi1);
end

if max(u1) <= max(mkeqoi1)
    maxv = max(mkeqoi);
else
    maxv = max(u1);
end
minv = fix(minv);
maxv = ceil(maxv);
dx = 0.001
range = minv:dx:maxv;
%Trapezoidal rule for integration
P = pdf(pd6, range);
Q = pdf(pds6, range);
klint = P.*log(P./Q);
kldiv = (klint(1) + klint(length(klint))+ 2*sum(klint(2:length(klint)-1)))*dx/2;
%%
save('pdf.mat', 'klms', '-append')
%%
%load('pdf.mat', 'x7', 'p7','x','p')
figure(1)
hold all

plot(x6,p6, '--')
plot(xs6,ps6)

title('KDE Comparison - 49 MC Samples') % of Time Centered Average Kinetic Energy Around t = 400 for 100 R.S. of the for Wave Amplitude, Fixed KE AVG')
xlabel('KDE Support')
ylabel('Density')
legend('Monte Carlo KDE 49 Samples', 'Level 4 Surrogate 49 Samples')

%%
load('pdf.mat', 'x','p','p1','x1', 'p2','x2', 'p3', 'x3', 'p4', 'x4', 'p5', 'x5', 'p6', 'x6', 'p7','x7')
%%
figure(2)
hold all

%x = min(mkeqoi1)-2:1:max(mkeqoi1) + 2;

plot(x,p, '-*', x1,p1, x2,p2, x3,p3, x4, p4, x5, p6, x6, p6, x7, p7)

%pd = fitdist(u1, 'Kernel')
% x = min(u1)-2:1:max(u1) + 2;
% p = pdf(pd,x);
% plot(x,p)
title('KDE Convergence Graph Level 1-7 Surrogates vs Monte Carlo Simulation') % of Time Centered Average Kinetic Energy Around t = 400 for 100 R.S. of the for Wave Amplitude, Fixed KE AVG')
legend('Monte Carlo KDE')
ylabel('Density')
xlabel('KDE Support')
%% Error analysis graph plotting theoretical polynomial interpolation error vs sparse grid level, and MSE vs sparse grid level
%%
hold all
pd = fitdist(mkeqoi2', 'Kernel')
x = 73:1:105;
p = pdf(pd,x);
plot(x,p)

pd = fitdist(u2, 'Kernel')
x = 73:1:100;
p = pdf(pd,x);
plot(x,p)

%%
figure(3)
x = linspace(1,7,7)
plot(x, klv, '*')
title('Plot of Kulbach-Liebler Divergence vs Sparse Grid Level')
xlabel('Surrogate Level')
ylabel('KL Divergence')
figure(4)
plot(gridsize, klv, '*')
title('Plot of K-L Divergence Vs. Grid Size (Number of Grid Points)')
xlabel('Sparse Grid Size')
ylabel('K-L Divergence')
%% Plotting theoretical error vs sparse grid level and actual error vs sparse grid level
%Find the normalized error
hold all
x = linspace(1,7, 7);
norme = error./rssq(error);
etheory = 2.^(-2*x);
normet = etheory./rssq(etheory);
let = log2(normet);
le = log2(norme);
figure(1)
title( ' Interpolation Normalized MSE Comparison')
plot(x, norme, '*', x, normet, 'o')
legend('Experimental Error', 'Theoretical Error')
title('Theoretical & Experimential Interpolation Error')
xlabel('Sparse Grid Level')
ylabel('Error')

figure(2)
plot(x, le, '*', x, let, 'o')
title( 'Normalized Semi-Log Error Graph')
title('Theoretical & Experimential Interpolation Error')
xlabel('Sparse Grid Level')
ylabel('Error')
%%
figure('Name','Kinetic Energy','Numbertitle','off');
for i = 1:10
    semilogy(time, kineticenergyi(:,i));
    xlabel('time'); %grid on;
    ylabel('kinetic energy');
    hold all
end
%%
x = 1:1:7
plot(x,error,'*')