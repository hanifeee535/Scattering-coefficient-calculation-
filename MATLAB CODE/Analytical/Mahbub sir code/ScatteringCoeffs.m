clc
clear all;
close all;

% Setup
% fsz = 20; %for the oral presentation
% fsz=11;  % for the proceedings
lw = 2.8;
fs = 20;

% Working directory
current_dir = cd;
cd(current_dir);

global k1 k3 K1 K3 a rho_i rho_o
f  = 10 .^([-2:.01:7]);   
nf = length(f);
                   
a=0.3e-6;               % particle radius
rho_o=997;               % density of water
rho_i=2100;              % density of particle (silica)
rho_cap=rho_i/rho_o;
vo=1497;                 %outside speed (water)
vi=5968;                 % inside speed (silica)
mu0_s=30.9e9;
mu_doubleprime= 0.000891;

eta=0;   % shear viscosity
att_o=0.023e-12;
att_i=2.6e-22;

k1=2*pi*f./vo+1i*att_o*f.^2;
K1=2*pi*f./vi+1i*att_i*f.^2;
k3=sqrt(rho_o*pi*f./mu_doubleprime)*(1+1i);
K3=sqrt(rho_i)*(2*pi*f)./sqrt(mu0_s);
% n_max=10; % maximum mode

xc=k1*a;
xs=k3*a;
Xc=K1*a;
Xs=K3*a;

% Bessel functions
[j0xc, j0pxc,h0xc,h0pxc] = SpherBess(0, xc);
[j0Xc, j0pXc,h0Xc,h0pXc] = SpherBess(0, Xc);

[j1xc, j1pxc,h1xc,h1pxc] = SpherBess(1, xc);
[j1Xc, j1pXc,h1Xc,h1pXc] = SpherBess(1, Xc);

[j1xs, j1pxs,h1xs,h1pxs] = SpherBess(1, xs);
[j1Xs, j1pXs,h1Xs,h1pXs] = SpherBess(1, Xs);

% elastic scattering coefficient: Exact analytical forms 

T0CC  = (rho_cap./Xs.^2).*(Xs.^2.*j0Xc + 4*Xc.*j0pXc).*(xc.*j0pxc);
T0CC  = T0CC - (1./xs.^2).*(xs.^2.*j0xc + 4*xc.*j0pxc).*(Xc.*j0pXc);
T0CCD = (1./xs.^2).*(xs.^2.*h0xc + 4*xc.*h0pxc).*(Xc.*j0pXc);
T0CCD = T0CCD - (rho_cap./Xs.^2).*(Xs.^2.*j0Xc + 4*Xc.*j0pXc).*(xc.*h0pxc);
T0CCe  = T0CC./T0CCD;

T1CC = (1i.*xc.^3).*(h1xs-(xs.*h1pxs)).*(rho_cap-1);
T1CCD = 3*(((4*rho_cap-7)*h1xs)+ ((1+2*rho_cap)*xs.*h1pxs));
T1CCe = T1CC./T1CCD;

T1CS = (rho_cap).*xc;
T1CSD = ((4*rho_cap-7)*h1xs)+ ((1+2*rho_cap)*xs.*h1pxs);
T1CSe = T1CS./T1CSD;

%rigit particle
T0CCr  = -j0pxc./h0pxc;  % rigid particle: Exact form

Num = (xc.*j1pxc).*((xs.*h1pxs)+h1xs);
Num = Num - (2*j1xc.*h1xs);
Den = (xc.*h1pxc).*((xs.*h1pxs)+h1xs);
Den = Den- (2*h1xc.*h1xs);
T1CCr = -(Num./Den);  %rigit particle 

Num2 = (xc.*j1pxc).*h1pxc;
Num2 = Num2- (j1xc.*xc.*h1pxc);
Den2 = (xc.*h1pxc).*((xs.*h1pxs)+h1xs);
Den2 = Den2-(2*h1xc.*h1xs);
T1CSr = -(Num2./Den2);




XT=10.^[-3, +2]; 

%T0CC
figure('NumberTitle','on', 'Name','T_0^CC');
plot (real(xs), real (T0CCe./xc.^3), 'Color', '[0 1 1]', 'LineStyle','-','LineWidth',lw);hold on;
plot (real(xs), real (T0CCr./xc.^3), 'Color', 'k', 'LineStyle','-','LineWidth',lw);
set(gca, 'XScale','log', 'YScale','log', 'FontSize',fs);
set(gca, 'XTick',XT, 'XLim',[min(XT) max(XT)]);

xlabel ('\Ree(k_{s}a)', 'FontWeight','Bold', 'FontSize',fs);
ylabel('T_{0}^{CC}', 'FontWeight','Bold', 'FontSize',fs);
legend('(T_{0}^{CCe})','(T_{0}^{CCr})');



figure('NumberTitle','on', 'Name','T_0^CC');
plot (real(xs), imag (T0CCe./xc.^3), 'Color', '[0 1 1]', 'LineStyle','-','LineWidth',lw);hold on;
plot (real(xs), imag (T0CCr./xc.^3), 'Color', 'k', 'LineStyle','-','LineWidth',lw);
set(gca, 'XScale','log', 'YScale','log', 'FontSize',fs);
set(gca, 'XTick',XT, 'XLim',[min(XT) max(XT)]);

xlabel ('\Ree(k_{s}a)', 'FontWeight','Bold', 'FontSize',fs);
ylabel('T_{0}^{CC}', 'FontWeight','Bold', 'FontSize',fs);
legend('\Imm(T_{0}^{CCe})','\Imm(T_{0}^{CCr})');

%T1CC

figure('NumberTitle','on', 'Name','T_1^CC');
plot (real(xs), real (T1CCe./xc.^3), 'Color', '[0 1 1]', 'LineStyle','-','LineWidth',lw);hold on;
plot (real(xs), real (T1CCr./xc.^3), 'Color', 'k', 'LineStyle','-','LineWidth',lw);
set(gca, 'XScale','log', 'YScale','log', 'FontSize',fs);
set(gca, 'XTick',XT, 'XLim',[min(XT) max(XT)]);

xlabel ('\Ree(k_{s}a)', 'FontWeight','Bold', 'FontSize',fs);
ylabel('T_{1}^{CC}', 'FontWeight','Bold', 'FontSize',fs);
legend('(T_{1}^{CCe})','(T_{1}^{CCr})');



figure('NumberTitle','on', 'Name','T_1^CC');
plot (real(xs), imag (T1CCe./xc.^3), 'Color', '[0 1 1]', 'LineStyle','-','LineWidth',lw);hold on;
plot (real(xs), imag (T1CCr./xc.^1), 'Color', 'k', 'LineStyle','-','LineWidth',lw);
set(gca, 'XScale','log', 'YScale','log', 'FontSize',fs);
set(gca, 'XTick',XT, 'XLim',[min(XT) max(XT)]);

xlabel ('\Ree(k_{s}a)', 'FontWeight','Bold', 'FontSize',fs);
ylabel('T_{1}^{CC}', 'FontWeight','Bold', 'FontSize',fs);
legend('\Imm(T_{1}^{CCe})','\Imm(T_{1}^{Cr})');

%T1CS

figure('NumberTitle','on', 'Name','T_1^CS');
plot (real(xs), real (T1CSe.*h1xs./xc.^3), 'Color', '[0 1 1]', 'LineStyle','-','LineWidth',lw);hold on;
plot (real(xs), real (T1CSr./xc.^1), 'Color', 'k', 'LineStyle','-','LineWidth',lw);
set(gca, 'XScale','log', 'YScale','log', 'FontSize',fs);
set(gca, 'XTick',XT, 'XLim',[min(XT) max(XT)]);

xlabel ('\Ree(k_{s}a)', 'FontWeight','Bold', 'FontSize',fs);
ylabel('T_{1}^{CS}', 'FontWeight','Bold', 'FontSize',fs);
legend('(T_{1}^{CSe})','(T_{1}^{CSr})');



figure('NumberTitle','on', 'Name','T_1^CS');
plot (real(xs), imag (T1CSe.*h1xs./xc.^3), 'Color', '[0 1 1]', 'LineStyle','-','LineWidth',lw);hold on;
plot (real(xs), imag (T1CSr./xc.^3), 'Color', 'k', 'LineStyle','-','LineWidth',lw);
set(gca, 'XScale','log', 'YScale','log', 'FontSize',fs);
set(gca, 'XTick',XT, 'XLim',[min(XT) max(XT)]);

xlabel ('\Ree(k_{s}a)', 'FontWeight','Bold', 'FontSize',fs);
ylabel('T_{1}^{CS}', 'FontWeight','Bold', 'FontSize',fs);
legend('\Imm(T_{1}^{CSe})','\Imm(T_{1}^{CSr})');



function [jn,jnprime,hn1,hn1prime]=SpherBess(n,x)
%spherical bessel and hankel function of order n and their argument x and
%their derivatives
sq=sqrt(pi ./(2*x));
jn=sq.*besselj(n+0.5,x);
jnprime=sq.*((n./x).* besselj(n+0.5,x)-besselj(n+1.5,x));

hn1=sq.* besselh(n+0.5,1,x);
hn1prime=sq.*((n./x).* besselh(n+0.5,1,x)-besselh(n+1.5,1,x));


end %end of bessel function
