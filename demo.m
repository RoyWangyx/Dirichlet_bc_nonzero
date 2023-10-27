%Number of grid points on [0,1]^2 
%i.e. uniform mesh with step h=1/(s-1)
s = 256;

%Create mesh (only needed for plotting)
[X,Y] = meshgrid(0:(1/(s-1)):1);

%Parameters of covariance C = tau^(2*alpha-2)*(-Laplacian + tau^2 I)^(-alpha)
%Note that we need alpha > d/2 (here d= 2) 
%Laplacian has zero Neumann boundry
%alpha and tau control smoothness; the bigger they are, the smoother the
%function
alpha = 2;
tau = 3;

%Generate random coefficients from N(0,C)
norm_a = GRF(alpha, tau, s);

%Exponentiate it, so that a(x) > 0
%Now a ~ Lognormal(0, C)
%This is done so that the PDE is elliptic
lognorm_a = exp(norm_a);


%Forcing function, f(x) = 0 g random iid gaussian
f = zeros(s,s);
g =ones(s,s);
g_rand = randn(s);
%Solve PDE: - div(a(x)*grad(p(x))) = f(x)
lognorm_p = solve_gwf(lognorm_a,f,g);
lognorm_p_rand = solve_gwf(lognorm_a,f,g_rand);

%Plot coefficients and solutions
subplot(2,2,1)
surf(X,Y,lognorm_a); 
view(2); 
shading interp;
colorbar;
subplot(2,2,2)
surf(X,Y,lognorm_p); 
view(2); 
shading interp;
colorbar;
subplot(2,2,3)
surf(X,Y,g_rand); 
view(2); 
shading interp;
colorbar;
subplot(2,2,4)
surf(X,Y,lognorm_p_rand); 
view(2); 
shading interp;
colorbar;
