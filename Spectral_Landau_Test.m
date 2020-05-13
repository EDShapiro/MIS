%Random variable creation for brute force evaluation of solutions to
%Vlasov-Poisson - Nonlinear Landau Damping 
samples1 = 49;
lower1=0; upper1 = 1;
dist = 'Uniform';
eps1 = makedist(dist,0,1);
%t = truncate(pd, lower, upper);

%%
samples2 =49;
dist2 = 'Uniform';
lower2 = 0.0; upper2 = 1;
eps2 = makedist(dist2, 0.0, 1);
%eps2 = truncate(eps2, lower2, upper2);

%%
eps1d = random(eps1, samples1, 1);
eps2d = random(eps2, samples2, 1);
%%
epsm49 = [eps1d,eps2d];
%%
samples3 =10;
dist3 = 'Normal'
lower3=0; upper3 = 1;
eps3 = makedist(dist,'mu', 0, 'sigma', 0.1);
eps3 = truncate(eps3, lower3, upper3);
%t = truncate(pd, lower, upper);
eps3d = random(eps3, samples3, 1);

%%
%Modified script to incorporate uncertainty in background field and
%instability
%eps=0.5; % Amplitude of perturbation, 0.05 for linear, 0.5 for nonlinear
kx=0.5;    % Wave vector
L=2*pi/kx; % length of domain
qm=-1;    % negative charge to mass = electrons, $\frac{q}{m}=-1$

% Initial condition
f0 =@(x,v,eps1d, eps2d)   (1+eps1d*cos((kx)*x))./sqrt(2*pi).*exp(-0.5.*(v+eps2d).^2);
 
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
kineticenergyr=zeros(numt,samples1);
meankineticenergy=zeros(samples1);
sh_entropy=zeros(numt,1); %Shannon entropy
kl_entropy=zeros(numt,1); %Kullback-Leibler entropy

EE=zeros(Nx,numt);


% Set initial condition
%%
%Defining a matrix for all of the distributions based off of the perturbed
%fields

for i = 1:samples1
   
    fr(:,:,i)=f0(XX,VV,epsm49(i,1), epsm49(i,2));
    fr(:,:,i)=fft(fr(:,:,i),[],1); 
end

%%

%%
%Unit Test for UQ
%figure('Name', 'Phase Space Density', 'Numbertitle','off');
%Adding loop to consider perturbations in the electric field
for i =1:samples1
    for tdx=1:numt
    for sdx=1:length(rksd)
        rho=sum(fr(:,:,i),2).*(vmax-vmin)/Nv; %rho(x) -> integrate over v
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
            fr(:,:,i)=ifft(fr(:,:,i),[],1,'symmetric'); %back transform
        % f is fully backtransformed, calculate diagnostics
        if (sdx==length(rksd))
            
                kineticenergyr(tdx,i)=sum(fr(:,:,i),1)*v.^2*(vmax-vmin)/Nv;
                %sh_entropy(tdx)=sum(sum(f1.*log(abs(f1))))*L/Nx*(vmax-vmin)/Nv;
            %kl_entropy(tdx)=sum(sum(f1.*log(abs(f1./f0(XX,VV)))))...
                %*L/Nx*(vmax-vmin)/Nv;
            
        end

        % Fourier transform in v
        fr(:,:,i)=fft(fr(:,:,i),[],2);
        % Build Matrix for spatial advection
         m=fftshift((1:Nv)  -Nv/2-1);
         VSHIFT=exp((1j*dt*rksc(sdx)*2*pi/(vmax-vmin))*E.*m);

        fr(:,:,i)=fr(:,:,i).*VSHIFT;
        fr(:,:,i)=ifft(fr(:,:,i),[],2,'symmetric');
        %Advection in X
        fr(:,:,i)=fft(fr(:,:,i),[],1);
        fr(:,:,i)=fr(:,:,i).*XSHIFT(:,:,sdx);
    end
        
    end
    %Plot 1 kinetic energy sample from a given pertubration

end
%%
KE49 = kineticenergyr;
%%
%Setting up arrays to save simulation outputs to perform statistical
%analysis on, evaluate on the sparse grid interpolant for comparison, and
%whatever else is needed.
save ('RandLD1.mat', 'KE321', 'epsm321', '-append' )

%% 
keqoir = kineticenergyr(390:410,:,:);
mkeqoir = zeros(100,1);
for i = 1:100
    mkeqoir(i) = mean(keqoir(:,i,1));
end
%%
%Creating a vector of the mean electron kinetic energy values, the
%quantitiy of interest in this UQ experiment.
meankineticenergyr = zeros(samples1,1);
  for i = 1:samples1
        meankineticenergyr(i) = mean(keqoir(1,i,:));
  end



%%

 %%
figure('Name','Kinetic Energy','Numbertitle','off');
for i = 1:10
    for j = 1:3
    semilogy(time, kineticenergyr(:,i,j));
    xlabel('time'); %grid on;
    ylabel('kinetic energy');
    hold all
    end
end
%%
veci = zeros(samples2,1);

for i = 1:samples2
    veci(i) = keqoir(1,1,i);
end

%%
figure('Name','Kinetic Energy @ t = 400','Numbertitle','off');
    plot(eps2d, v, '*');
    xlabel('Perturbation in Amplitude'); %grid on;
    ylabel('kinetic energy');
    hold all
    