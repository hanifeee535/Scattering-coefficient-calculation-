function ScatteringCoefficients(plot_reim)
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


nf = length(f);
k1 = 2*pi*f./vo + 1i*att_0*f.^2;
K1 = 2*pi*f./vi + 1i*att_i*f.^2;
k3 = (1+1i)* sqrt(rho_o*2*pi*f./2*mu_pp);
K3 = sqrt(rho_i)*(2*pi*f)./sqrt(mu0_s);

media.rho_i = rho_i;
media.rho_o = rho_o;
order_N = 2;


for kf = 1:nf
   if (fix(kf/100)==(kf/100)), clc; disp(sprintf('%d%%',round(100*kf/nf))); end
    media.a = a;
    media.k1 = k1(kf);
    media.k3 = k3(kf);
    media.K1 = K1(kf);
    media.K3 = K3(kf);
    [T] = CoefficientThermalNeglected(order_N,media);
    
   %coefficients
   T0CC(kf) = T(1,1);
   
   T1CC(kf) = T(2,1);
   T1CS(kf) = T(2,2);
   T1SC(kf) = T(2,3);
   T1SS(kf) = T(2,4);
   
   T2CC(kf) = T(3,1);
   T2CS(kf) = T(3,2);
   T2SC(kf) = T(3,3);
   T2SS(kf) = T(3,4);
   
end

c1 = k1.*a;
c3 = k3.*a;
C1 = K1.*a;
C3 = K3.*a;

[j0c, j0pc, h0c, h0pc] = SpherBess (0,c1);
[j0s, j0ps, h0s, h0ps] = SpherBess (0,c3);
[j0C, j0pC, h0C, h0pC] = SpherBess (0,C1);
[j0S, j0pS, h0S, h0pS] = SpherBess (0,C3);

[j2c, j2pc, h2c, h2pc] = SpherBess (2,c1);
[j2s, j2ps, h2s, h2ps] = SpherBess (2,c3);
[j2C, j2pC, h2C, h2pC] = SpherBess (2,C1);
[j2S, j2pS, h2S, h2pS] = SpherBess (2,C3);

[j1c, j1pc, h1c, h1pc] = SpherBess (1,c1);
[j1s, j1ps, h1s, h1ps] = SpherBess (1,c3);
[j1C, j1pC, h1C, h1pC] = SpherBess (1,C1);
[j1S, j1pS, h1S, h1pS] = SpherBess (1,C3);


% analytical expressions
AA = ((rho_cap./(C3.^2)).*(((C3.^2).*j0C)+(4*C1.*j0pC)).*c1.*j0pc);
BB = ( (1./c3.^2).*(((c3.^2).*j0c)+(4*C1.*j0pC)).*C1.*j0pC);
T0CCs = AA./BB;

CC =  (1i*(c1.^3).*(h1s-(c3.*h1ps)).*(rho_cap-1));
DD = (3*(((4*rho_cap)-7).*h1s)+(c3.*(1+2*rho_cap).*h1ps));
T1CCs = CC./DD; 
T1CSs = -((rho_cap-1)*c1)./ (((4*rho_cap-7)*h1ps)+ (1+2*rho_cap)*c3.*h1ps);
T1SCs = (2*(rho_cap -1)*(c1.^2))./ c3.*((((4*rho_cap)-7)*h1ps)+(((2*rho_cap)+1)*c3.*h1ps)); 
T1SSs = -((((4*rho_cap)-7)*j1s) + (((2*rho_cap)+1)*c3.*j1ps))./ ((((4*rho_cap)-7)*h1s) + (((2*rho_cap)+1)*c3.*h1ps));

T2CCs = 1i*(2/135)*(c1.^5).*(c3.*h2ps - 2*h2s)./ (c3.*h2ps+3*h2s);
T2CSs = -(c1.^2)./(9*(c3.*h2ps+3*h2s));
T2SCs = -2*c1.^3./(3*c3.*(c3.*h2ps +3*h2s));
T2SSs = -(c3.*j2ps + 3*j2s)./(c3.*h2ps +3*h2s);

calculation_time = toc;

disp(sprintf('calculation time = %d s',round (calculation_time)));

%Range of variable

yc_min = min (real(c1));
yc_max = max (real(c1));

ys_min = min (real(c3));
ys_max = max (real(c3));

ysp_min = min (real(C3));
ysp_max = max (real(C3));
XT = 10.^[-3,3];

plot_reim = 0;

switch plot_reim  
    %T2CC
    
    
    
    case {0} %real /imaginary parts on a single figure 
        
         %T0CC
        figure('NumberTitle','off', 'Name','T_0^CC');
        plot (real(c3), real (T0CC), 'Color', '[0 1 1]', 'LineStyle','-','LineWidth',lw);hold on;
        
        plot (real(c3), imag (T0CC), 'Color', 'b', 'LineStyle','-','LineWidth',lw);
        
        plot (real(c3), real (T0CCs), 'Color', 'k', 'LineStyle','-','LineWidth',lw);
        
        plot (real(c3), imag (T0CCs), 'Color', '[1 0 1]', 'LineStyle','-','LineWidth',lw);
        set(gca, 'XScale','log', 'YScale','log', 'FontSize',fs);
        set(gca, 'XTick',XT, 'XLim',[min(XT) max(XT)]);
        xlabel ('\Ree(k_{s}a)', 'FontWeight','Bold', 'FontSize',fs);
        ylabel('T_{0}^{CC}', 'FontWeight','Bold', 'FontSize',fs);
        legend('\Ree(T_{0}^{CC})','\Imm(T_{0}^{CC})','\Ree(T_{0s}^{CC})','\Imm(T_{0s}^{CC})');
        
         %T1CC
        figure('NumberTitle','off', 'Name','T_1^CC');
        plot (real(c3), real (T1CC./c1.^3), 'Color', 'g', 'LineStyle','-','LineWidth',lw);hold on;
        
        plot (real(c3), imag (T1CC./c1.^3), 'Color', 'b', 'LineStyle','-','LineWidth',lw);
        
        plot (real(c3), real (T1CCs./c1.^3), 'Color', 'k', 'LineStyle','-','LineWidth',lw);
        
        plot (real(c3), imag (T1CCs./c1.^3), 'Color', 'y', 'LineStyle','-','LineWidth',lw);
        set(gca, 'XScale','log', 'YScale','log', 'FontSize',fs);
        set(gca, 'XTick',XT, 'XLim',[min(XT) max(XT)]);
        xlabel ('\Ree(k_{s}a)', 'FontWeight','Bold', 'FontSize',fs);
        ylabel('T_{1}^{CC}', 'FontWeight','Bold', 'FontSize',fs);
        %legend('\Ree(T_{1}^{CC})','\Imm(T_{1}^{CC})','\Ree(T_{1s}^{CC})','\Imm(T_{1s}^{CC})');
        
         %T1Cs
        figure('NumberTitle','off', 'Name','T_1^CS');
        plot (real(c3), real (T1CS./h1s./c1), 'Color', '[0 1 1]', 'LineStyle','-.','LineWidth',lw);hold on;
        
        plot (real(c3), imag (T1CS./h1s./c1), 'Color', 'b', 'LineStyle','-.','LineWidth',lw);
        
        plot (real(c3), real (T1CSs./h1s./c1), 'Color', 'k', 'LineStyle','-','LineWidth',lw);
        
        plot (real(c3), imag (T1CSs./h1s./c1), 'Color', '[1 0 1]', 'LineStyle','-','LineWidth',lw);
        set(gca, 'XScale','log', 'YScale','log', 'FontSize',fs);
        set(gca, 'XTick',XT, 'XLim',[min(XT) max(XT)]);
        xlabel ('\Ree(k_{s}a)', 'FontWeight','Bold', 'FontSize',fs);
        ylabel('T_{1}^{CS}', 'FontWeight','Bold', 'FontSize',fs);
        legend('\Ree(T_{1}^{CS})','\Imm(T_{1}^{CS})','\Ree(T_{1s}^{CS})','\Imm(T_{1s}^{CS})');
        
         %T1SC
        figure('NumberTitle','off', 'Name','T_1^SC');
        plot (real(c3), real (T1SC./h1s./c1.^2), 'Color', '[0 1 1]', 'LineStyle','-.','LineWidth',lw);hold on;
        
        plot (real(c3), imag (T1SC./h1s./c1.^2), 'Color', 'b', 'LineStyle','-.','LineWidth',lw);
        
        plot (real(c3), real (T1SCs./h1s./c1.^2), 'Color', 'k', 'LineStyle','-','LineWidth',lw);
        
        plot (real(c3), imag (T1SCs./h1s./c1.^2), 'Color', '[1 0 1]', 'LineStyle','-','LineWidth',lw);
        set(gca, 'XScale','log', 'YScale','log', 'FontSize',fs);
        set(gca, 'XTick',XT, 'XLim',[min(XT) max(XT)]);
        xlabel ('\Ree(k_{s}a)', 'FontWeight','Bold', 'FontSize',fs);
        ylabel('T_{1}^{SC}', 'FontWeight','Bold', 'FontSize',fs);
        legend('\Ree(T_{1}^{SC})','\Imm(T_{1}^{SC})','\Ree(T_{1s}^{SC})','\Imm(T_{1s}^{SC})');
        
         %T1SS
        figure('NumberTitle','off', 'Name','T_1^SS');
        plot (real(c3), real (T1SS./j1s./h1s), 'Color', '[0 1 1]', 'LineStyle','-.','LineWidth',lw);hold on;
        
        plot (real(c3), imag (T1SS./j1s./h1s), 'Color', 'b', 'LineStyle','-.','LineWidth',lw);
        
        plot (real(c3), real (T1SSs./j1s./h1s), 'Color', 'k', 'LineStyle','-','LineWidth',lw);
        
        plot (real(c3), imag (T1SSs./j1s./h1s), 'Color', '[1 0 1]', 'LineStyle','-','LineWidth',lw);
        set(gca, 'XScale','log', 'YScale','log', 'FontSize',fs);
        set(gca, 'XTick',XT, 'XLim',[min(XT) max(XT)]);
        xlabel ('\Ree(k_{s}a)', 'FontWeight','Bold', 'FontSize',fs);
        ylabel('T_{1}^{SS}', 'FontWeight','Bold', 'FontSize',fs);
        legend('\Ree(T_{1}^{SS})','\Imm(T_{1}^{SS})','\Ree(T_{1s}^{SS})','\Imm(T_{1s}^{SS})');
        
        
    
    %T2CC
        figure('NumberTitle','off', 'Name','T_2N^CC');
        plot (real(c3), real (T2CC./c1.^5), 'Color', '[0 1 1]', 'LineStyle','-.','LineWidth',lw);hold on;
        
        plot (real(c3), imag (T2CC./c1.^5), 'Color', 'b', 'LineStyle','-.','LineWidth',lw);
        
        plot (real(c3), real (T2CCs./c1.^5), 'Color', 'k', 'LineStyle','-','LineWidth',lw);
        
        plot (real(c3), imag (T2CCs./c1.^5), 'Color', '[1 0 1]', 'LineStyle','-','LineWidth',lw);
        set(gca, 'XScale','log', 'YScale','log', 'FontSize',fs);
        set(gca, 'XTick',XT, 'XLim',[min(XT) max(XT)]);
        xlabel ('\Ree(k_{s}a)', 'FontWeight','Bold', 'FontSize',fs);
        ylabel('T_{2}^{CC}', 'FontWeight','Bold', 'FontSize',fs);
        legend('\Ree(T_{2}^{CC})','\Imm(T_{2}^{CC})','\Ree(T_{2s}^{CC})','\Imm(T_{2s}^{CC})');
        
        %T2CS
        figure('NumberTitle','off', 'Name','T_2N^CS');
        plot (real(c3), real (T2CS.*h2s./c1.^2), 'Color', '[0 1 1]', 'LineStyle','-.','LineWidth',lw);hold on;
        
        plot (real(c3), imag (T2CS.*h2s./c1.^2), 'Color', 'b', 'LineStyle','-.','LineWidth',lw);
        
        plot (real(c3), real (T2CSs.*h2s./c1.^2), 'Color', 'k', 'LineStyle','-','LineWidth',lw);
        
        plot (real(c3), imag (T2CSs.*h2s./c1.^2), 'Color', '[1 0 1]', 'LineStyle','-','LineWidth',lw);
        set(gca, 'XScale','log', 'YScale','log', 'FontSize',fs);
        set(gca, 'XTick',XT, 'XLim',[min(XT) max(XT)]);
        xlabel ('\Ree(k_{s}a)', 'FontWeight','Bold', 'FontSize',fs);
        ylabel('T_{2}^{Cs}', 'FontWeight','Bold', 'FontSize',fs);
        legend('\Ree(T_{2}^{Cs})','\Imm(T_{2}^{Cs})','\Ree(T_{2s}^{Cs})','\Imm(T_{2s}^{Cs})');
        
        %T2SC
        figure('NumberTitle','off', 'Name','T_2N^CS');
        plot (real(c3), real (T2SC./j2s./c1.^3), 'Color', '[0 1 1]', 'LineStyle','-.','LineWidth',lw);hold on;
        
        plot (real(c3), imag (T2SC./j2s./c1.^3), 'Color', 'b', 'LineStyle','-.','LineWidth',lw);
        
        plot (real(c3), real (T2SCs./j2s./c1.^3), 'Color', 'k', 'LineStyle','-','LineWidth',lw);
        
        plot (real(c3), imag (T2SCs./j2s./c1.^3), 'Color', '[1 0 1]', 'LineStyle','-','LineWidth',lw);
        set(gca, 'XScale','log', 'YScale','log', 'FontSize',fs);
        set(gca, 'XTick',XT, 'XLim',[min(XT) max(XT)]);
        xlabel ('\Ree(k_{s}a)', 'FontWeight','Bold', 'FontSize',fs);
        ylabel('T_{2}^{sc}', 'FontWeight','Bold', 'FontSize',fs);
        legend('\Ree(T_{2}^{sc})','\Imm(T_{2}^{sc})','\Ree(T_{2s}^{sc})','\Imm(T_{2s}^{sc})');
        
        
        %T2SS
        figure('NumberTitle','off', 'Name','T_2N^SS');
        plot (real(c3), real (T2SS.*h2s./j2s), 'Color', '[0 1 1]', 'LineStyle','-.','LineWidth',lw);hold on;
        
        plot (real(c3), imag (T2SS.*h2s./j2s), 'Color', 'b', 'LineStyle','-.','LineWidth',lw);
        
        plot (real(c3), real (T2SSs.*h2s./j2s), 'Color', 'k', 'LineStyle','-','LineWidth',lw);
        
        plot (real(c3), imag (T2SSs.*h2s./j2s), 'Color', '[1 0 1]', 'LineStyle','-','LineWidth',lw);
        set(gca, 'XScale','log', 'YScale','log', 'FontSize',fs);
        set(gca, 'XTick',XT, 'XLim',[min(XT) max(XT)]);
        xlabel ('\Ree(k_{s}a)', 'FontWeight','Bold', 'FontSize',fs);
        ylabel('T_{2}^{ss}', 'FontWeight','Bold', 'FontSize',fs);
        legend('\Ree(T_{2}^{ss})','\Imm(T_{2}^{ss})','\Ree(T_{2s}^{ss})','\Imm(T_{2s}^{ss})');
        
    
    
        
    case {2}
        %modulas/phase on two figures
        
        
        %T0CC
        
        figure('NumberTitle','off', 'Name',' abs T_0^CC');
        plot (real(c3), real (T0CC), 'Color', '[0 1 1]', 'LineStyle','-','LineWidth',lw);hold on;
        plot (real(c3), real (T0CCs), 'Color', 'k', 'LineStyle','-','LineWidth',lw);
        set(gca, 'XScale','log', 'YScale','log', 'FontSize',fs);
        set(gca, 'XTick',XT, 'XLim',[min(XT) max(XT)]);
        xlabel ('\Ree(k_{s}a)', 'FontWeight','Bold', 'FontSize',fs);
        ylabel(' T_{0}^{CC}', 'FontWeight','Bold', 'FontSize',fs);
        legend('(T_{0}^{CC})','(T_{0s}^{CC})');
        
        %arg
        figure('NumberTitle','off', 'Name',' arg T_0^CC');
        plot (real(c3), imag (T0CC), 'Color', '[0 1 1]', 'LineStyle','-','LineWidth',lw);hold on;
        plot (real(c3), imag (T0CCs), 'Color', 'k', 'LineStyle','-','LineWidth',lw);
        set(gca, 'XScale','log', 'YScale','log', 'FontSize',fs);
        set(gca, 'XTick',XT, 'XLim',[min(XT) max(XT)]);
        xlabel ('\Ree(k_{s}a)', 'FontWeight','Bold', 'FontSize',fs);
        ylabel(' T_{0}^{CC}', 'FontWeight','Bold', 'FontSize',fs);
        legend('(Arg(T_{0}^{CC})','Arg(T_{0s}^{CC})');
        
         %T1CC
        
        figure('NumberTitle','off', 'Name',' abs T_1^CC');
        plot (real(c3), real (T1CC./c1.^3), 'Color', '[0 1 1]', 'LineStyle','-','LineWidth',lw);hold on;
        plot (real(c3), real (T1CCs./c1.^3), 'Color', 'k', 'LineStyle','-','LineWidth',lw);
        set(gca, 'XScale','log', 'YScale','log', 'FontSize',fs);
        set(gca, 'XTick',XT, 'XLim',[min(XT) max(XT)]);
        xlabel ('\Ree(k_{s}a)', 'FontWeight','Bold', 'FontSize',fs);
        ylabel(' T_{1}^{CC}', 'FontWeight','Bold', 'FontSize',fs);
        legend('(T_{1}^{CC})','(T_{1s}^{CC})');
        
        %arg
        figure('NumberTitle','off', 'Name',' arg T_1^CC');
        plot (real(c3), imag (T1CC./c1.^3), 'Color', '[0 1 1]', 'LineStyle','-','LineWidth',lw);hold on;
        plot (real(c3), imag (T1CCs./c1.^3), 'Color', 'k', 'LineStyle','-','LineWidth',lw);
        set(gca, 'XScale','log', 'YScale','log', 'FontSize',fs);
        set(gca, 'XTick',XT, 'XLim',[min(XT) max(XT)]);
        xlabel ('\Ree(k_{s}a)', 'FontWeight','Bold', 'FontSize',fs);
        ylabel(' T_{1}^{CC}', 'FontWeight','Bold', 'FontSize',fs);
        legend('(Arg(T_{1}^{CC})','Arg(T_{1s}^{CC})');
        
        
         %T1CS
        
        figure('NumberTitle','off', 'Name',' abs T_1^CS');
        plot (real(c3), real (T1CS./h1s./c1), 'Color', '[0 1 1]', 'LineStyle','-','LineWidth',lw);hold on;
        plot (real(c3), real (T1CSs./h1s./c1), 'Color', 'k', 'LineStyle','-','LineWidth',lw);
        set(gca, 'XScale','log', 'YScale','log', 'FontSize',fs);
        set(gca, 'XTick',XT, 'XLim',[min(XT) max(XT)]);
        xlabel ('\Ree(k_{s}a)', 'FontWeight','Bold', 'FontSize',fs);
        ylabel(' T_{1}^{CS}', 'FontWeight','Bold', 'FontSize',fs);
        legend('(T_{1}^{CS})','(T_{1s}^{CS})');
        
        %arg
        figure('NumberTitle','off', 'Name',' arg T_1^CS');
        plot (real(c3), imag (T1CS./h1s./c1), 'Color', '[0 1 1]', 'LineStyle','-','LineWidth',lw);hold on;
        plot (real(c3), imag (T1CSs./h1s./c1), 'Color', 'k', 'LineStyle','-','LineWidth',lw);
        set(gca, 'XScale','log', 'YScale','log', 'FontSize',fs);
        set(gca, 'XTick',XT, 'XLim',[min(XT) max(XT)]);
        xlabel ('\Ree(k_{s}a)', 'FontWeight','Bold', 'FontSize',fs);
        ylabel(' T_{1}^{CS}', 'FontWeight','Bold', 'FontSize',fs);
        legend('(Arg(T_{1}^{CS})','Arg(T_{1s}^{CS})');
        
         %T1SC
        
        figure('NumberTitle','off', 'Name',' abs T_1^SC');
        plot (real(c3), real (T1SC./h1s./c1.^2), 'Color', '[0 1 1]', 'LineStyle','-','LineWidth',lw);hold on;
        plot (real(c3), real (T1SCs ./h1s./c1.^5), 'Color', 'k', 'LineStyle','-','LineWidth',lw);
        set(gca, 'XScale','log', 'YScale','log', 'FontSize',fs);
        set(gca, 'XTick',XT, 'XLim',[min(XT) max(XT)]);
        xlabel ('\Ree(k_{s}a)', 'FontWeight','Bold', 'FontSize',fs);
        ylabel(' T_{1}^{SC}', 'FontWeight','Bold', 'FontSize',fs);
        legend('(T_{1}^{SC})','(T_{1s}^{SC})');
        
        %arg
        figure('NumberTitle','off', 'Name',' arg T_1^SC');
        plot (real(c3), imag (T1SC ./h1s./c1.^5), 'Color', '[0 1 1]', 'LineStyle','-','LineWidth',lw);hold on;
        plot (real(c3), imag (T1SCs ./h1s./c1.^5), 'Color', 'k', 'LineStyle','-','LineWidth',lw);
        set(gca, 'XScale','log', 'YScale','log', 'FontSize',fs);
        set(gca, 'XTick',XT, 'XLim',[min(XT) max(XT)]);
        xlabel ('\Ree(k_{s}a)', 'FontWeight','Bold', 'FontSize',fs);
        ylabel(' T_{1}^{SC}', 'FontWeight','Bold', 'FontSize',fs);
        legend('(Arg(T_{1}^{SC})','Arg(T_{1s}^{SC})');
        
        %T1SS
        
        figure('NumberTitle','off', 'Name',' abs T_1^SS');
        plot (real(c3), real (T1SS./j1s./h1s), 'Color', '[0 1 1]', 'LineStyle','-','LineWidth',lw);hold on;
        plot (real(c3), real (T1SSs./j1s./h1s), 'Color', 'k', 'LineStyle','-','LineWidth',lw);
        set(gca, 'XScale','log', 'YScale','log', 'FontSize',fs);
        set(gca, 'XTick',XT, 'XLim',[min(XT) max(XT)]);
        xlabel ('\Ree(k_{s}a)', 'FontWeight','Bold', 'FontSize',fs);
        ylabel(' T_{1}^{SS}', 'FontWeight','Bold', 'FontSize',fs);
        legend('(T_{1}^{SS})','(T_{1s}^{SS})');
        
        %arg
        figure('NumberTitle','off', 'Name',' arg T_1^SS');
        plot (real(c3), imag (T1SS./j1s./h1s), 'Color', '[0 1 1]', 'LineStyle','-','LineWidth',lw);hold on;
        plot (real(c3), imag (T1SSs./j1s./h1s), 'Color', 'k', 'LineStyle','-','LineWidth',lw);
        set(gca, 'XScale','log', 'YScale','log', 'FontSize',fs);
        set(gca, 'XTick',XT, 'XLim',[min(XT) max(XT)]);
        xlabel ('\Ree(k_{s}a)', 'FontWeight','Bold', 'FontSize',fs);
        ylabel(' T_{1}^{SS}', 'FontWeight','Bold', 'FontSize',fs);
        legend('(Arg(T_{1}^{SS})','Arg(T_{1s}^{SS})');
        
        
        
        
        
        
        %T2CC
        
        figure('NumberTitle','off', 'Name',' abs T_2N^CC');
        plot (real(c3), real (T2CC./c1.^5), 'Color', '[0 1 1]', 'LineStyle','-','LineWidth',lw);hold on;
        plot (real(c3), real (T2CCs./c1.^5), 'Color', 'k', 'LineStyle','-','LineWidth',lw);
        set(gca, 'XScale','log', 'YScale','log', 'FontSize',fs);
        set(gca, 'XTick',XT, 'XLim',[min(XT) max(XT)]);
        xlabel ('\Ree(k_{s}a)', 'FontWeight','Bold', 'FontSize',fs);
        ylabel(' T_{2N}^{CC}', 'FontWeight','Bold', 'FontSize',fs);
        legend('(T_{2N}^{CC})','(T_{2Ns}^{CC})');
        
        %arg
        figure('NumberTitle','off', 'Name',' arg T_2N^CC');
        plot (real(c3), imag (T2CC./c1.^5), 'Color', '[0 1 1]', 'LineStyle','-','LineWidth',lw);hold on;
        plot (real(c3), imag (T2CCs./c1.^5), 'Color', 'k', 'LineStyle','-','LineWidth',lw);
        set(gca, 'XScale','log', 'YScale','log', 'FontSize',fs);
        set(gca, 'XTick',XT, 'XLim',[min(XT) max(XT)]);
        xlabel ('\Ree(k_{s}a)', 'FontWeight','Bold', 'FontSize',fs);
        ylabel(' T_{2N}^{CC}', 'FontWeight','Bold', 'FontSize',fs);
        legend('(Arg(T_{2N}^{CC})','Arg(T_{2Ns}^{CC})');
        
        
        %T2CS
        
        figure('NumberTitle','off', 'Name',' abs T_2N^CS');
        plot (real(c3), real (T2CS.*h2s./c1.^2), 'Color', '[0 1 1]', 'LineStyle','-','LineWidth',lw);hold on;
        plot (real(c3), real (T2CSs.*h2s./c1.^2), 'Color', 'k', 'LineStyle','-','LineWidth',lw);
        set(gca, 'XScale','log', 'YScale','log', 'FontSize',fs);
        set(gca, 'XTick',XT, 'XLim',[min(XT) max(XT)]);
        xlabel ('\Ree(k_{s}a)', 'FontWeight','Bold', 'FontSize',fs);
        ylabel(' T_{2N}^{CS}', 'FontWeight','Bold', 'FontSize',fs);
        legend('(T_{2N}^{CS})','(T_{2Ns}^{CS})');
        
        %arg
        figure('NumberTitle','off', 'Name',' arg T_2N^CS');
        plot (real(c3), imag (T2CS.*h2s./c1.^2), 'Color', '[0 1 1]', 'LineStyle','-','LineWidth',lw);hold on;
        plot (real(c3), imag (T2CSs.*h2s./c1.^2), 'Color', 'k', 'LineStyle','-','LineWidth',lw);
        set(gca, 'XScale','log', 'YScale','log', 'FontSize',fs);
        set(gca, 'XTick',XT, 'XLim',[min(XT) max(XT)]);
        xlabel ('\Ree(k_{s}a)', 'FontWeight','Bold', 'FontSize',fs);
        ylabel(' T_{2N}^{CS}', 'FontWeight','Bold', 'FontSize',fs);
        legend('(Arg(T_{2N}^{CS})','Arg(T_{2Ns}^{CS})');
        
        %T2SC
        
        figure('NumberTitle','off', 'Name',' abs T_2N^SC');
        plot (real(c3), real (T2SC./j2s./c1.^3), 'Color', '[0 1 1]', 'LineStyle','-','LineWidth',lw);hold on;
        plot (real(c3), real (T2SCs./j2s./c1.^3), 'Color', 'k', 'LineStyle','-','LineWidth',lw);
        set(gca, 'XScale','log', 'YScale','log', 'FontSize',fs);
        set(gca, 'XTick',XT, 'XLim',[min(XT) max(XT)]);
        xlabel ('\Ree(k_{s}a)', 'FontWeight','Bold', 'FontSize',fs);
        ylabel(' T_{2N}^{SC}', 'FontWeight','Bold', 'FontSize',fs);
        legend('(T_{2N}^{SC})','(T_{2Ns}^{SC})');
        
        %arg
        figure('NumberTitle','off', 'Name',' arg T_2N^SC');
        plot (real(c3), imag (T2SC./j2s./c1.^3), 'Color', '[0 1 1]', 'LineStyle','-','LineWidth',lw);hold on;
        plot (real(c3), imag (T2SCs./j2s./c1.^3), 'Color', 'k', 'LineStyle','-','LineWidth',lw);
        set(gca, 'XScale','log', 'YScale','log', 'FontSize',fs);
        set(gca, 'XTick',XT, 'XLim',[min(XT) max(XT)]);
        xlabel ('\Ree(k_{s}a)', 'FontWeight','Bold', 'FontSize',fs);
        ylabel(' T_{2N}^{SC}', 'FontWeight','Bold', 'FontSize',fs);
        legend('(Arg(T_{2N}^{SC})','Arg(T_{2Ns}^{SC})');
        
        %T2SS
        
        figure('NumberTitle','off', 'Name',' abs T_2N^SS');
        plot (real(c3), real (T2SS.*h2s./j2s), 'Color', '[0 1 1]', 'LineStyle','-','LineWidth',lw);hold on;
        plot (real(c3), real (T2SSs.*h2s./j2s), 'Color', 'k', 'LineStyle','-','LineWidth',lw);
        set(gca, 'XScale','log', 'YScale','log', 'FontSize',fs);
        set(gca, 'XTick',XT, 'XLim',[min(XT) max(XT)]);
        xlabel ('\Ree(k_{s}a)', 'FontWeight','Bold', 'FontSize',fs);
        ylabel(' T_{2N}^{SS}', 'FontWeight','Bold', 'FontSize',fs);
        legend('(T_{2N}^{SS})','(T_{2Ns}^{SS})');
        
        %arg
        figure('NumberTitle','off', 'Name',' arg T_2N^SS');
        plot (real(c3), imag (T2SS.*h2s./j2s), 'Color', '[0 1 1]', 'LineStyle','-','LineWidth',lw);hold on;
        plot (real(c3), imag (T2SSs.*h2s./j2s), 'Color', 'k', 'LineStyle','-','LineWidth',lw);
        set(gca, 'XScale','log', 'YScale','log', 'FontSize',fs);
        set(gca, 'XTick',XT, 'XLim',[min(XT) max(XT)]);
        xlabel ('\Ree(k_{s}a)', 'FontWeight','Bold', 'FontSize',fs);
        ylabel(' T_{2N}^{SC}', 'FontWeight','Bold', 'FontSize',fs);
        legend('(Arg(T_{2N}^{SS})','Arg(T_{2Ns}^{SS})');
end
end
        
        
        
        
        
        
        


function [T] = CoefficientThermalNeglected (NN,media)
T = [];
a = media.a;
k1 = media.k1;
k3 = media.k3;
K1 = media.K1;
K3 = media.K3;
rho_i = media.rho_i;
rho_o = media.rho_o;
rho_cap = rho_i/rho_o;
c1 = k1*a;
c3 = k3*a;
C1 = K1*a;
C3 = K3*a;
for m = 1:2 %m =1 for compressional incident and m = 2 for shear incidend
    
    T1 = []; %for compressional
    T3 = []; %shear
    for n = 0:NN %n is current order and NN is maximum order
        
        %spherical bessel and hankel functions 
     [jn1, jn_prime_1, hn1, hn_prime_1] =  SpherBess (n,c1);
     [jn3, jn_prime_3, hn3, hn_prime_3] =  SpherBess (n,c3);
     [Jn1, Jn_prime_1, Hn1, Hn_prime_1] =  SpherBess (n,C1);
     [Jn3, Jn_prime_3, Hn3, Hn_prime_3] =  SpherBess (n,C3);
     %Matrix components set up
     
     M11 = c1*hn_prime_1;
     M13 = n*(n+1)*hn3;
     M14 = -(C1*jn_prime_1);
     M16 = -(n*(n+1)*jn3);
     
     M21 = hn1;
     M23 = (c3*hn_prime_3)+hn3;
     M24 = -(Jn1);
     M26 = -c3*Jn_prime_3-Jn3;
     
     M31=-(K3/k3)^2*((c3^2-2*n*(n+1))*hn1+4*c1*hn_prime_1);
     M33=2*(K3/k3)^2*n*(n+1)*(c3*hn_prime_3-hn3);
     M34=(rho_i/rho_o)*((C3^2-2*n*(n+1))*Jn1+4*C1*Jn_prime_1);
     M36=-(2*rho_i/rho_o)*n*(n+1)*(C3*Jn_prime_3-Jn3);
     
     M41=2*(K3/k3)^2*(c1*hn_prime_1-hn1);
     M43=(K3/k3)^2*((2*n*(n+1)-2-c3^2)*hn3-2*c3*hn_prime_3);
     M44=-(2*rho_i/rho_o)*(C1*Jn_prime_1-Jn1);
     M46=-(rho_i/rho_o)*((2*n*(n+1)-2-C3^2)*Jn3-2*C3*Jn_prime_3);
     
     D11=-c1*jn_prime_1;
     D13=-n*(n+1)*jn3;
     
     D21=-jn1;
     D23=-c3*jn_prime_3-jn3;
     
     D31=(K3/k3)^2*((c3^2-2*n*(n+1))*jn1+4*c1*jn_prime_1);
     D33=-2*(K3/k3)^2*n*(n+1)*(c3*jn_prime_3-jn3);
     
     D41=-2*(K3/k3)^2*(c1*jn_prime_1-jn1);
     D43=-(K3/k3)^2*((2*n*(n+1)-2-c3^2)*jn3-2*c3*jn_prime_3);
     
     D=[D11 D13; D21 D23; D31 D33; D41 D43];
     M=[M11 M13 M14 M16; M21 M23 M24 M26; M31 M33 M34 M36;M41 M43 M44 M46];
     
     m1=Order_grand(abs(mean(M(:,1))));
     m2=Order_grand(abs(mean(M(:,2))));
     m3=Order_grand(abs(mean(M(:,3))));
     m4=Order_grand(abs(mean(M(:,4))));
     m5=Order_grand(abs(mean(D(:,1))));
     m6=Order_grand(abs(mean(D(:,2))));
     
     M(:,1)=M(:,1)*10^(-m1);
     M(:,2)=M(:,2)*10^(-m2);
     M(:,3)=M(:,3)*10^(-m3);
     M(:,4)=M(:,4)*10^(-m4);
     D(:,1)=D(:,1)*10^(-m5);
     D(:,2)=D(:,2)*10^(-m6);
     
     format long e
     
     if n == 0
        M=[M11*10^(-m1) M14 *10^(-m4); M31*10^(-m1) M34*10^(-m4)];
        D=[D11*10^(-m5) ; D31*10^(-m5)];
        X=linsolve(M, D);
         
         switch m
             case 1
                 TC=X(1)*10^(m5-m1); % T0CC
                 TS=0; % T0CS
             case 2
                 TC=0; % T0SC
                 TS=0; % T0SS
         end
     else 
         X=linsolve(M, D(:,m));
         switch m
             case 1 
                 TC=X(1)*10^(m5-m1); %TCC
                 TS=X(2)*10^(m5-m2); %TCS
             case 2
                 TC=X(1)*10^(m6-m1); %TCC
                 TS=X(2)*10^(m6-m2); %TCC
         end
     end %end if n==0
     
     T1=[T1;TC];
     T3=[T3;TS];
     
     
     
     
     
        
    end %end n loop over n
    T=[T,T1,T3];
end %end m loop
    
end  %end function

function  [m]=Order_grand(x)
y=abs(x);
m=0;
if y<1

    while y<10^(m-1)
        m=m-1;
    end
    if abs(y-10^m)>abs(y-10^(m-1))
       m=m-1;
    end
else 
    while y>10^(m+1)
        m=m+1;
    end
    if abs(y-10^m)>abs(y-10^(m+1))
       m=m+1;

    end
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
