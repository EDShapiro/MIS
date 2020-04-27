dim = 2;
grid = sgpp.Grid.createLinearClenshawCurtisGrid(dim);
level = 3;
grid.getGenerator().regular(level);
gridStorage = grid.getStorage();

%%
y = zeros(gridStorage.getSize(),dim);
    for i = 1:gridStorage.getSize() -1
           gp = gridStorage.getPoint(i);
           y(i,:) = [gp.getStandardCoordinate(0), gp.getStandardCoordinate(1)];
    end
%%
% Sparse grid, grid point conversion for sparse grid surrogate construction
% of Vlasov-Poisson Nonlinear Landau Damping Problem 
%Convert [0,1] to [a,b] for each random variable (f^(-1) = x*(b-a) +a)
lower1 = -1; upper1 = 1;
lower2= 0.0 ; upper2 = pi/2;
xconv = zeros(length(y),2);
xconv(:,1) = y(:,1)*(upper1 - lower1) + lower1; %Scaling sparse grid points to be from [lower, upper]
xconv(:,2) = y(:,2)*(upper2 - lower2) + lower2;
eps1d = xconv(:,1);
eps2d = xconv(:,2);
samples1 = length(xconv(:,1));
samples2 = length(xconv(:,2));
%%
%Modified script to incorporate uncertainty in background field and
%instability
%eps=0.5; % Amplitude of perturbation, 0.05 for linear, 0.5 for nonlinear
kx=0.5;    % Wave vector
L=2*pi/kx; % length of domain
qm=-1;    % negative charge to mass = electrons, $\frac{q}{m}=-1$

% Initial condition
f0 =@(x,v,eps1d, eps2d)   (1+eps1d.*cos((kx+eps2d)*x))./sqrt(2*pi).*exp(-0.5.*(v).^2);
 
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
kineticenergy=zeros(numt,samples1,samples2);
meankineticenergy=zeros(samples1,samples2, 1);
sh_entropy=zeros(numt,1); %Shannon entropy
kl_entropy=zeros(numt,1); %Kullback-Leibler entropy

EE=zeros(Nx,numt);


% Set initial condition
%%
%Defining a matrix for all of the distributions based off of the perturbed
%fields

for i = 1:samples1
   for j =1:samples2
    f(:,:,i,j)=f0(XX,VV,eps1d(i), eps2d(j));
    f(:,:,i,j)=fft(f(:,:,i,j),[],1);
   end
end

% for i = 1:samples
%     b(:,:,i) = background(
% end

%%
%Unit Test for UQ
%figure('Name', 'Phase Space Density', 'Numbertitle','off');
%Adding loop to consider perturbations in the electric field
%Unit Test for UQ
%figure('Name', 'Phase Space Density', 'Numbertitle','off');
%Adding loop to consider perturbations in the electric field
for i =1:samples1
    for j = 1:samples2
    
    for tdx=1:numt
    for sdx=1:length(rksd)
        rho=sum(f(:,:,i,j),2).*(vmax-vmin)/Nv; %rho(x) -> integrate over v
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
            f(:,:,i,j)=ifft(f(:,:,i,j),[],1,'symmetric'); %back transform
        % f is fully backtransformed, calculate diagnostics
        if (sdx==length(rksd))
            
                kineticenergy(tdx,i,j)=sum(f(:,:,i,j),1)*v.^2*(vmax-vmin)/Nv;
                %sh_entropy(tdx)=sum(sum(f1.*log(abs(f1))))*L/Nx*(vmax-vmin)/Nv;
            %kl_entropy(tdx)=sum(sum(f1.*log(abs(f1./f0(XX,VV)))))...
                %*L/Nx*(vmax-vmin)/Nv;
            
        end

        % Fourier transform in v
        f(:,:,i,j)=fft(f(:,:,i,j),[],2);
        % Build Matrix for spatial advection
         m=fftshift((1:Nv)  -Nv/2-1);
         VSHIFT=exp((1j*dt*rksc(sdx)*2*pi/(vmax-vmin))*E.*m);

        f(:,:,i,j)=f(:,:,i,j).*VSHIFT;
        f(:,:,i,j)=ifft(f(:,:,i,j),[],2,'symmetric');
        %Advection in X
        f(:,:,i,j)=fft(f(:,:,i,j),[],1);
        f(:,:,i,j)=f(:,:,i,j).*XSHIFT(:,:,sdx);
    end
        
    end
    %Plot 1 kinetic energy sample from a given pertubration

    end
end
%%
t_vec = linspace(0,50, 500);
plot(t_vec, kineticenergy(:,1,4))
%%
keqoi = kineticenergy(400,:,:); %Setting up an array for the quantity of interest at a time point (t=50),
% when the system is in a steady state, across all sparse grid points.
%%
% Creating a coefficients vector for the interpolant - Running more than
% once with different data causes the script to crash
alpha = sgpp.DataVector(gridStorage.getSize());
alpha.setAll(0);
for k = 1:gridStorage.getSize() 
    for i =1:samples1
        for j = 1:samples2
            alpha.set(k, keqoi(i,j));
        end
    end
end

%%
% Mapping random sample points to grid values, and then creating a 2-D grid
% of input values for the interpolant.
load('RandLD.mat', 'eps1d', 'eps2d') %Loading the random variables from brute force simulation

e1 = zeros(length(eps1d)); e2 = length(eps2d);
z = cat(eps1d,eps2d);

e1 = (eps1d(i) + 1)./2;

e2 = 2/pi*eps2d;

%%
p = sgpp.DataVector(dim);
p.set(0, 0.7516);
p.set(1, 0.1265);
opEval = sgpp.createOperationEvalNaive(grid);
fprintf('u(0.52, 0.73) = %.4f\n', opEval.eval(alpha,p));
%%
save ('SGLD.mat', 'rho', 'E', 'kineticenergy', 'meankineticenergy', 'f', 'eps1d', 'eps2d', 'eps3d')
%%
figure('Name','Kinetic Energy','Numbertitle','off');
for i = 1:10
    semilogy(time, kineticenergy(:,i));
    xlabel('time'); %grid on;
    ylabel('kinetic energy');
    hold all
end
