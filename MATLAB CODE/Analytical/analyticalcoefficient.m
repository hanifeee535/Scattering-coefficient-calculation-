function analyticalcoefficient(plot_reim)
clc
clear all
close all

tic

%set up the plots
lw = 2;
fs =16;

%physical properties of fluid and solid
a = 300e-6;
rho_o = 997;
rho_i = 2100;
rho_cap = rho_i/rho_o;
vo = 1497;
vi = 5968;
mu0_s = 30.9e9;
mu_pp = .000891;
eta = 0;
att_0 = .023e-12;
att_i = 2.6e-22;

%freuency vector
f = 10.^([-3:.2:6.15]);
n = 0:2;

nf = length(f);
k1 = 2*pi*f./vo + 1i*att_0*f.^2;
K1 = 2*pi*f./vi + 1i*att_i*f.^2;
k3 = (1+1i)* sqrt(rho_o*2*pi*f./2*mu_pp);
K3 = sqrt(rho_i)*(2*pi*f)./sqrt(mu0_s);

c1 = k1.*a;
c3 = k3.*a;
C1 = K1.*a;
C3 = K3.*a;


%analytical calculation 

if n == 0
    [jnc, jnpc, hnc, hnpc] = SpherBess (n,c1);
    [jns, jnps, hns, hnps] = SpherBess (n,c3);
    [jnC, jnpC, hnC, hnpC] = SpherBess (n,C1);
    [jnS, jnpS, hnS, hnpS] = SpherBess (n,C3);


    M0 = -1i*((n/(2*n+1)*(rho_cap-1))-((rho_cap./C3.^2).*2*n*(n-1)));
    E0 = -1i*(2*n+1)./((n*hns)-(c3.*hnps));
    N0 = ((rho_cap./(2*n+1))-((rho_cap./C3.^2).*2*(n-1))+ ((1./C3.^2).*2*(n^2-1)/(2*n+1))).*(hns+ c3.*hnps) + (n+1)/(2*n+1).*((2*n*(n-1)./c3.^2)-1).*hns;

    T0CC = -(((n/(2*n+1))*(rho_cap -1) + (2*n*(n-1))*(rho_cap./C3.^2)).*c1.^(2*n+1))./ ((M0+E0).*N0).*(factorial(2*n)/(2^n*factorial(n))).*(factorial(2*n+1)/(2^n*factorial(n)));
    T0CS = (1i*((n/(2*n+1)).*(rho_cap -1) + (2*n*(n-1)).*(rho_cap./C3.^2)).*c1.^n)./ ((Mn+En).*Nn).* n.*((n*hns) - (c3.*hnps)).*(factorial(2*n)/(2^n*factorial(n)));
    T0CCr= (1i*(c3.*hnps)+hns-((n+1)*hns).*n.*c1.^(2*n+1))./ ((c3.*hnps+hns+(n.*hns)).*(n+1).*(factorial(2*n+1)/(2^n*factorial(n))).*(factorial(2*n)/(2^n*factorial(n))));
    T0CSr= -c1.^n./(((c3.*hnps)+hns+(n*hns)).*(n+1)*(2*n+1)*(factorial(2*n)/(2^n*factorial(n))));
    
    
    
elseif n==1
    [jnc, jnpc, hnc, hnpc] = SpherBess (n,c1);
    [jns, jnps, hns, hnps] = SpherBess (n,c3);
    [jnC, jnpC, hnC, hnpC] = SpherBess (n,C1);
    [jnS, jnpS, hnS, hnpS] = SpherBess (n,C3);


    M1 = -1i*((n/(2*n+1)*(rho_cap-1))-((rho_cap./C3.^2).*2*n*(n-1)));
    E1 = -1i*(2*n+1)./((n*hns)-(c3.*hnps));
    N1 = ((rho_cap./(2*n+1))-((rho_cap./C3.^2).*2*(n-1))+ ((1./C3.^2).*2*(n^2-1)/(2*n+1))).*(hns+ c3.*hnps) + (n+1)/(2*n+1).*((2*n*(n-1)./c3.^2)-1).*hns;

    T1CC = -(((n/(2*n+1))*(rho_cap -1) + (2*n*(n-1))*(rho_cap./C3.^2)).*c1.^(2*n+1))./ ((M1+E1).*N1).*(factorial(2*n)/(2^n*factorial(n))).*(factorial(2*n+1)/(2^n*factorial(n)));
    T1CS = (1i*((n/(2*n+1)).*(rho_cap -1) + (2*n*(n-1)).*(rho_cap./C3.^2)).*c1.^n)./ ((Mn+En).*Nn).* n.*((n*hns) - (c3.*hnps)).*(factorial(2*n)/(2^n*factorial(n)));
    T1CCr= (1i*(c3.*hnps)+hns-((n+1)*hns).*n.*c1.^(2*n+1))./ ((c3.*hnps+hns+(n.*hns)).*(n+1).*(factorial(2*n+1)/(2^n*factorial(n))).*(factorial(2*n)/(2^n*factorial(n))));
    T1CSr= -c1.^n./(((c3.*hnps)+hns+(n*hns)).*(n+1)*(2*n+1)*(factorial(2*n)/(2^n*factorial(n))));
    
    
elseif n ==2
    [jnc, jnpc, hnc, hnpc] = SpherBess (n,c1);
    [jns, jnps, hns, hnps] = SpherBess (n,c3);
    [jnC, jnpC, hnC, hnpC] = SpherBess (n,C1);
    [jnS, jnpS, hnS, hnpS] = SpherBess (n,C3);


    M2 = -1i*((n/(2*n+1)*(rho_cap-1))-((rho_cap./C3.^2).*2*n*(n-1)));
    E2 = -1i*(2*n+1)./((n*hns)-(c3.*hnps));
    N2 = ((rho_cap./(2*n+1))-((rho_cap./C3.^2).*2*(n-1))+ ((1./C3.^2).*2*(n^2-1)/(2*n+1))).*(hns+ c3.*hnps) + (n+1)/(2*n+1).*((2*n*(n-1)./c3.^2)-1).*hns;

    T2CC = -(((n/(2*n+1))*(rho_cap -1) + (2*n*(n-1))*(rho_cap./C3.^2)).*c1.^(2*n+1))./ ((M2+E2).*N2).*(factorial(2*n)/(2^n*factorial(n))).*(factorial(2*n+1)/(2^n*factorial(n)));
    T2CS = (1i*((n/(2*n+1)).*(rho_cap -1) + (2*n*(n-1)).*(rho_cap./C3.^2)).*c1.^n)./ ((Mn+En).*Nn).* n.*((n*hns) - (c3.*hnps)).*(factorial(2*n)/(2^n*factorial(n)));
    T2CCr= (1i*(c3.*hnps)+hns-((n+1)*hns).*n.*c1.^(2*n+1))./ ((c3.*hnps+hns+(n.*hns)).*(n+1).*(factorial(2*n+1)/(2^n*factorial(n))).*(factorial(2*n)/(2^n*factorial(n))));
    T2CSr= -c1.^n./(((c3.*hnps)+hns+(n*hns)).*(n+1)*(2*n+1)*(factorial(2*n)/(2^n*factorial(n))));
    
    
end

calculation_time = toc;

disp(sprintf('calculation time = %d s',round (calculation_time)));

%Range of variable

yc_min = min (real(c1));
yc_max = max (real(c1));

ys_min = min (real(c3));
ys_max = max (real(c3));

ysp_min = min (real(C3));
ysp_max = max (real(C3));
XT = 10.^[-3,2];
plot_reim = 0;

switch plot_reim  
    %T2CC
    
    
    
    case {0} %real /imaginary parts on a single figure 
        
         %T0CC
        figure('NumberTitle','off', 'Name','T_0^CC');
        plot (real(c3), real (T0CC), 'Color', '[0 1 1]', 'LineStyle','-','LineWidth',lw);hold on;
        
        plot (real(c3), imag (T0CC), 'Color', 'b', 'LineStyle','-','LineWidth',lw);
        
        plot (real(c3), real (T0CCr), 'Color', 'k', 'LineStyle','-','LineWidth',lw);
        
        plot (real(c3), imag (T0CCr), 'Color', '[1 0 1]', 'LineStyle','-','LineWidth',lw);
        set(gca, 'XScale','log', 'YScale','log', 'FontSize',fs);
        set(gca, 'XTick',XT, 'XLim',[min(XT) max(XT)]);
        xlabel ('\Ree(k_{s}a)', 'FontWeight','Bold', 'FontSize',fs);
        ylabel('T_{0}^{CC}', 'FontWeight','Bold', 'FontSize',fs);
        legend('\Ree(T_{0}^{CCe})','\Imm(T_{0}^{CCe})','\Ree(T_{0r}^{CC})','\Imm(T_{0r}^{CC})');
        
         %T1CC
        figure('NumberTitle','off', 'Name','T_1^CC');
        plot (real(c3), real (T1CC./c1.^3), 'Color', 'g', 'LineStyle','-','LineWidth',lw);hold on;
        
        plot (real(c3), imag (T1CC./c1.^3), 'Color', 'b', 'LineStyle','-','LineWidth',lw);
        
        plot (real(c3), real (T1CCr./c1.^3), 'Color', 'k', 'LineStyle','-','LineWidth',lw);
        
        plot (real(c3), imag (T1CCr./c1.^3), 'Color', 'y', 'LineStyle','-','LineWidth',lw);
        set(gca, 'XScale','log', 'YScale','log', 'FontSize',fs);
        set(gca, 'XTick',XT, 'XLim',[min(XT) max(XT)]);
        xlabel ('\Ree(k_{s}a)', 'FontWeight','Bold', 'FontSize',fs);
        ylabel('T_{1}^{CC}', 'FontWeight','Bold', 'FontSize',fs);
        legend('\Ree(T_{1}^{CCe})','\Imm(T_{1}^{CCe})','\Ree(T_{1s}^{CCr})','\Imm(T_{1s}^{CCr})');
        
         %T1Cs
        figure('NumberTitle','off', 'Name','T_1^CS');
        plot (real(c3), real (T1CS./h1s./c1), 'Color', '[0 1 1]', 'LineStyle','-.','LineWidth',lw);hold on;
        
        plot (real(c3), imag (T1CS./h1s./c1), 'Color', 'b', 'LineStyle','-.','LineWidth',lw);
        
        plot (real(c3), real (T1CSr./h1s./c1), 'Color', 'k', 'LineStyle','-','LineWidth',lw);
        
        plot (real(c3), imag (T1CSr./h1s./c1), 'Color', '[1 0 1]', 'LineStyle','-','LineWidth',lw);
        set(gca, 'XScale','log', 'YScale','log', 'FontSize',fs);
        set(gca, 'XTick',XT, 'XLim',[min(XT) max(XT)]);
        xlabel ('\Ree(k_{s}a)', 'FontWeight','Bold', 'FontSize',fs);
        ylabel('T_{1}^{CS}', 'FontWeight','Bold', 'FontSize',fs);
        legend('\Ree(T_{1}^{CSe})','\Imm(T_{1}^{CSe})','\Ree(T_{1s}^{CSr})','\Imm(T_{1s}^{CSr})');
        
         %T2CC
        figure('NumberTitle','off', 'Name','T_2N^CC');
        plot (real(c3), real (T2CC./c1.^5), 'Color', '[0 1 1]', 'LineStyle','-.','LineWidth',lw);hold on;
        
        plot (real(c3), imag (T2CC./c1.^5), 'Color', 'b', 'LineStyle','-.','LineWidth',lw);
        
        plot (real(c3), real (T2CCr./c1.^5), 'Color', 'k', 'LineStyle','-','LineWidth',lw);
        
        plot (real(c3), imag (T2CCr./c1.^5), 'Color', '[1 0 1]', 'LineStyle','-','LineWidth',lw);
        set(gca, 'XScale','log', 'YScale','log', 'FontSize',fs);
        set(gca, 'XTick',XT, 'XLim',[min(XT) max(XT)]);
        xlabel ('\Ree(k_{s}a)', 'FontWeight','Bold', 'FontSize',fs);
        ylabel('T_{2}^{CC}', 'FontWeight','Bold', 'FontSize',fs);
        legend('\Ree(T_{2}^{CCe})','\Imm(T_{2}^{CCe})','\Ree(T_{2s}^{CCe})','\Imm(T_{2s}^{CCe})');
        
        %T2CS
        figure('NumberTitle','off', 'Name','T_2N^CS');
        plot (real(c3), real (T2CS.*h2s./c1.^2), 'Color', '[0 1 1]', 'LineStyle','-.','LineWidth',lw);hold on;
        
        plot (real(c3), imag (T2CS.*h2s./c1.^2), 'Color', 'b', 'LineStyle','-.','LineWidth',lw);
        
        plot (real(c3), real (T2CSr.*h2s./c1.^2), 'Color', 'k', 'LineStyle','-','LineWidth',lw);
        
        plot (real(c3), imag (T2CSr.*h2s./c1.^2), 'Color', '[1 0 1]', 'LineStyle','-','LineWidth',lw);
        set(gca, 'XScale','log', 'YScale','log', 'FontSize',fs);
        set(gca, 'XTick',XT, 'XLim',[min(XT) max(XT)]);
        xlabel ('\Ree(k_{s}a)', 'FontWeight','Bold', 'FontSize',fs);
        ylabel('T_{2}^{Cs}', 'FontWeight','Bold', 'FontSize',fs);
        legend('\Ree(T_{2}^{CSe})','\Imm(T_{2}^{CSe})','\Ree(T_{2s}^{CSr})','\Imm(T_{2s}^{CSr})');
        case {2}
        %modulas/phase on two figures
        
        
        %T0CC
        
        figure('NumberTitle','off', 'Name',' abs T_0^CC');
        plot (real(c3), real (T0CC), 'Color', '[0 1 1]', 'LineStyle','-','LineWidth',lw);hold on;
        plot (real(c3), real (T0CCr), 'Color', 'k', 'LineStyle','-','LineWidth',lw);
        set(gca, 'XScale','log', 'YScale','log', 'FontSize',fs);
        set(gca, 'XTick',XT, 'XLim',[min(XT) max(XT)]);
        xlabel ('\Ree(k_{s}a)', 'FontWeight','Bold', 'FontSize',fs);
        ylabel(' T_{0}^{CC}', 'FontWeight','Bold', 'FontSize',fs);
        legend('(T_{0}^{CCs})','(T_{0s}^{CCr})');
        
        %arg
        figure('NumberTitle','off', 'Name',' arg T_0^CC');
        plot (real(c3), imag (T0CC), 'Color', '[0 1 1]', 'LineStyle','-','LineWidth',lw);hold on;
        plot (real(c3), imag (T0CCr), 'Color', 'k', 'LineStyle','-','LineWidth',lw);
        set(gca, 'XScale','log', 'YScale','log', 'FontSize',fs);
        set(gca, 'XTick',XT, 'XLim',[min(XT) max(XT)]);
        xlabel ('\Ree(k_{s}a)', 'FontWeight','Bold', 'FontSize',fs);
        ylabel(' T_{0}^{CC}', 'FontWeight','Bold', 'FontSize',fs);
        legend('(Arg(T_{0}^{CCe})','Arg(T_{0s}^{CCr})');
        
         %T1CC
        
        figure('NumberTitle','off', 'Name',' abs T_1^CC');
        plot (real(c3), real (T1CC./c1.^3), 'Color', '[0 1 1]', 'LineStyle','-','LineWidth',lw);hold on;
        plot (real(c3), real (T1CCr./c1.^3), 'Color', 'k', 'LineStyle','-','LineWidth',lw);
        set(gca, 'XScale','log', 'YScale','log', 'FontSize',fs);
        set(gca, 'XTick',XT, 'XLim',[min(XT) max(XT)]);
        xlabel ('\Ree(k_{s}a)', 'FontWeight','Bold', 'FontSize',fs);
        ylabel(' T_{1}^{CC}', 'FontWeight','Bold', 'FontSize',fs);
        legend('(T_{1}^{CCe})','(T_{1s}^{CCr})');
        
        %arg
        figure('NumberTitle','off', 'Name',' arg T_1^CC');
        plot (real(c3), imag (T1CC./c1.^3), 'Color', '[0 1 1]', 'LineStyle','-','LineWidth',lw);hold on;
        plot (real(c3), imag (T1CCr./c1.^3), 'Color', 'k', 'LineStyle','-','LineWidth',lw);
        set(gca, 'XScale','log', 'YScale','log', 'FontSize',fs);
        set(gca, 'XTick',XT, 'XLim',[min(XT) max(XT)]);
        xlabel ('\Ree(k_{s}a)', 'FontWeight','Bold', 'FontSize',fs);
        ylabel(' T_{1}^{CC}', 'FontWeight','Bold', 'FontSize',fs);
        legend('(Arg(T_{1}^{CCe})','Arg(T_{1s}^{CCr})');
        
        
         %T1CS
        
        figure('NumberTitle','off', 'Name',' abs T_1^CS');
        plot (real(c3), real (T1CS./h1s./c1), 'Color', '[0 1 1]', 'LineStyle','-','LineWidth',lw);hold on;
        plot (real(c3), real (T1CSr./h1s./c1), 'Color', 'k', 'LineStyle','-','LineWidth',lw);
        set(gca, 'XScale','log', 'YScale','log', 'FontSize',fs);
        set(gca, 'XTick',XT, 'XLim',[min(XT) max(XT)]);
        xlabel ('\Ree(k_{s}a)', 'FontWeight','Bold', 'FontSize',fs);
        ylabel(' T_{1}^{CS}', 'FontWeight','Bold', 'FontSize',fs);
        legend('(T_{1}^{CSe})','(T_{1s}^{CSr})');
        
        %arg
        figure('NumberTitle','off', 'Name',' arg T_1^CS');
        plot (real(c3), imag (T1CS./h1s./c1), 'Color', '[0 1 1]', 'LineStyle','-','LineWidth',lw);hold on;
        plot (real(c3), imag (T1CSr./h1s./c1), 'Color', 'k', 'LineStyle','-','LineWidth',lw);
        set(gca, 'XScale','log', 'YScale','log', 'FontSize',fs);
        set(gca, 'XTick',XT, 'XLim',[min(XT) max(XT)]);
        xlabel ('\Ree(k_{s}a)', 'FontWeight','Bold', 'FontSize',fs);
        ylabel(' T_{1}^{CS}', 'FontWeight','Bold', 'FontSize',fs);
        legend('(Arg(T_{1}^{CSe})','Arg(T_{1s}^{CSr})');
        
         
        
        
        
        
        
        
        %T2CC
        
        figure('NumberTitle','off', 'Name',' abs T_2N^CC');
        plot (real(c3), real (T2CC./c1.^5), 'Color', '[0 1 1]', 'LineStyle','-','LineWidth',lw);hold on;
        plot (real(c3), real (T2CCr./c1.^5), 'Color', 'k', 'LineStyle','-','LineWidth',lw);
        set(gca, 'XScale','log', 'YScale','log', 'FontSize',fs);
        set(gca, 'XTick',XT, 'XLim',[min(XT) max(XT)]);
        xlabel ('\Ree(k_{s}a)', 'FontWeight','Bold', 'FontSize',fs);
        ylabel(' T_{2N}^{CC}', 'FontWeight','Bold', 'FontSize',fs);
        legend('(T_{2N}^{CCe})','(T_{2Ns}^{CCr})');
        
        %arg
        figure('NumberTitle','off', 'Name',' arg T_2N^CC');
        plot (real(c3), imag (T2CC./c1.^5), 'Color', '[0 1 1]', 'LineStyle','-','LineWidth',lw);hold on;
        plot (real(c3), imag (T2CCr./c1.^5), 'Color', 'k', 'LineStyle','-','LineWidth',lw);
        set(gca, 'XScale','log', 'YScale','log', 'FontSize',fs);
        set(gca, 'XTick',XT, 'XLim',[min(XT) max(XT)]);
        xlabel ('\Ree(k_{s}a)', 'FontWeight','Bold', 'FontSize',fs);
        ylabel(' T_{2N}^{CC}', 'FontWeight','Bold', 'FontSize',fs);
        legend('(Arg(T_{2N}^{CCe})','Arg(T_{2Ns}^{CCr})');
        
        
        %T2CS
        
        figure('NumberTitle','off', 'Name',' abs T_2N^CS');
        plot (real(c3), real (T2CS.*h2s./c1.^2), 'Color', '[0 1 1]', 'LineStyle','-','LineWidth',lw);hold on;
        plot (real(c3), real (T2CSr.*h2s./c1.^2), 'Color', 'k', 'LineStyle','-','LineWidth',lw);
        set(gca, 'XScale','log', 'YScale','log', 'FontSize',fs);
        set(gca, 'XTick',XT, 'XLim',[min(XT) max(XT)]);
        xlabel ('\Ree(k_{s}a)', 'FontWeight','Bold', 'FontSize',fs);
        ylabel(' T_{2N}^{CS}', 'FontWeight','Bold', 'FontSize',fs);
        legend('(T_{2N}^{CSe})','(T_{2Ns}^{CSr})');
        
        %arg
        figure('NumberTitle','off', 'Name',' arg T_2N^CS');
        plot (real(c3), imag (T2CS.*h2s./c1.^2), 'Color', '[0 1 1]', 'LineStyle','-','LineWidth',lw);hold on;
        plot (real(c3), imag (T2CSr.*h2s./c1.^2), 'Color', 'k', 'LineStyle','-','LineWidth',lw);
        set(gca, 'XScale','log', 'YScale','log', 'FontSize',fs);
        set(gca, 'XTick',XT, 'XLim',[min(XT) max(XT)]);
        xlabel ('\Ree(k_{s}a)', 'FontWeight','Bold', 'FontSize',fs);
        ylabel(' T_{2N}^{CS}', 'FontWeight','Bold', 'FontSize',fs);
        legend('(Arg(T_{2N}^{CSe})','Arg(T_{2Ns}^{CSr})');
        
       
end
end

function [jn,jnprime,hn1,hn1prime]=SpherBess(n,x)
%spherical bessel and hankel function of order n and their argument x and
%their derivatives
sq=sqrt(pi ./(2*x));
jn=sq.*besselj(n+0.5,x);
jnprime=sq.*((n./x).* besselj(n+0.5,x)-besselj(n+1.5,x));

hn1=sq.* besselh(n+0.5,1,x);
hn1prime=sq.*((n./x).* besselh(n+0.5,1,x)-besselh(n+1.5,1,x));


end %end of bessel function