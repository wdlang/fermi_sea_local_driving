% tight binding model; single level evolution; weak pumping
% exact results compared with 1st order perburbation theory
% 2014. April 5 in dresden
clear all; close all; clc; tic; format long; myfont = 16;
myfont2 = 14;

T=3;
step = 30;
omega = 2*pi/T;
Tmax =790;
L = 512;                             %  length of the arms
Nspin = 2*L;                         %  number of spins
klist = -L : (L-1);
phiklist = 2*pi*klist/Nspin;
eklist = -2*cos(phiklist);

Dlist = [0.25, 0.8, 1.5];

% kilist  = [ 39, 40, 41,42,43 ];
%  kilist = [0, 1, 2, 3, 4];
% kilist = [10, 11, 12, 13, 14];
% kilist = [20, 21, 22, 23, 24];
% kilist = [30, 31, 32, 33, 34];
%  kilist = [50, 51, 52, 53, 54];
kilist = [ 61, 62, 63, 64, 65];
% kilist = 70+[0, 1, 2, 3, 4];

plist = zeros(length(Dlist), length(kilist),Tmax+1);
plist2 = zeros(length(Dlist), length(kilist),Tmax+1);

for sssss1 = 1:length(Dlist)
    sssss1
    Delta = Dlist(sssss1);
    
    Eilist = -2*cos(2*pi*kilist/Nspin)';
    Eflist = Eilist + omega;
    rate = 2*pi*4*Delta^2*omega./(pi^3*Nspin*sqrt(4- Eflist.^2));
      
    Tc = Nspin/sqrt(4-Eflist(1)^2);
    
    A0 = zeros(Nspin, Nspin);
    for s=1:Nspin-1
        A0(s,s+1)= -1;
        A0(s+1,s)= -1;
    end
    A0(Nspin,1)= -1;
    A0(1,Nspin)= -1;
    
    U = eye(Nspin, Nspin);
    for s = 1:step
        A1 = A0;
        A1(L+1, L+1) =  Delta*(cos(omega*(T/step)*s)- cos(omega*(T/step)*(s-1)))/(omega*(T/step));
        [V1, D1] = eig(A1);
        U1 = V1*diag(exp(-i*(T/step)*diag(D1)))*V1';
        U = U1*U;
    end
    
    vec = zeros(Nspin, length(kilist));
    for s11 = 1: length(kilist)
        ki = kilist(s11);
        for s = 1: Nspin
            vec(s, s11) = exp(i*2*pi*ki*s/Nspin)/sqrt(Nspin);
        end
    end
    
    vec0 = vec';
    
    A00 =  sparse(A0);
    for s = 0:Tmax
        plist (sssss1,:,s+1) = abs(diag(vec0*vec)).^2;
        
        vec = U*vec;
    end
    
    ek2 = sort(eklist);
    
    for sss = 1: length(kilist)
        
        Ef = Eflist (sss);
        index1 = find(Ef>ek2, 1,'last');
        index2 = find(Ef<ek2, 1,'first');
        
        E1 = ek2(index1);
        E2 = ek2(index2);
        delta = E2 - E1;
        
        eklist2 = E1 + delta*(-Nspin/2:(Nspin/2-1));
        
        for s=0:Tmax
            t = T*s;
            
            plist2(sssss1,sss,s+1) = 1- 2*4*((Delta/2/Nspin)^2)*sum((sin((eklist2-Ef).*(t/2)).^2)./((eklist2-Ef).^2));
        end
    end
end


h1=figure;
ha = tight_subplot(1,3,[.05 .05],[.7 .01],[.1 .01])

axes(ha(1))
plot(0:Tmax, reshape(plist(1,:,:),length(kilist),Tmax+1))
hold on
plot( 0:Tmax, reshape(plist2(1,:,:),length(kilist),Tmax+1), ':')
% plot((Tc/T)*ones(1,100), 0.01:0.01:1,':',(2*Tc/T)*ones(1,100), 0.01:0.01:1,':')
% ylim([floor(min(min(plist))*200)/200 1])
xlabel('$n$','fontsize',myfont,'Interpreter','Latex');
ylabel('$p$','fontsize',myfont,'Interpreter','Latex');
xlim([0 Tmax])
ylim([0.95 1])
xlims = xlim;
ylims = ylim;
a = 0.02; b = 0.1;
str = ['(a) $\Delta=$ ',num2str(Dlist(1))];
text((1-a)*xlims(1)+a*xlims(2), b*ylims(2)+(1-b)*ylims(1),str,'fontsize',myfont2,'Interpreter','Latex')
set(gca,'ytick',[0.95 1] )
set(gca,'yticklabel',{'0.95', '1'})

axes(ha(2))
plot(0:Tmax, reshape(plist(2,:,:),length(kilist),Tmax+1))
hold on
plot( 0:Tmax, reshape(plist2(2,:,:),length(kilist),Tmax+1), ':')
% plot((Tc/T)*ones(1,100), 0.01:0.01:1,':',(2*Tc/T)*ones(1,100), 0.01:0.01:1,':')
% ylim([floor(min(min(plist))*200)/200 1])
xlabel('$n$','fontsize',myfont,'Interpreter','Latex');
xlim([0 Tmax])
ylim([0.6 1])
xlims = xlim;
ylims = ylim;
str = ['(b) $\Delta=$ ',num2str(Dlist(2))];
text((1-a)*xlims(1)+a*xlims(2), b*ylims(2)+(1-b)*ylims(1),str,'fontsize',myfont2,'Interpreter','Latex')

axes(ha(3))
plot(0:Tmax, reshape(plist(3,:,:),length(kilist),Tmax+1))
hold on
plot( 0:Tmax, reshape(plist2(3,:,:),length(kilist),Tmax+1), ':')
% plot((Tc/T)*ones(1,100), 0.01:0.01:1,':',(2*Tc/T)*ones(1,100), 0.01:0.01:1,':')
% ylim([floor(min(min(plist))*200)/200 1])
xlabel('$n$','fontsize',myfont,'Interpreter','Latex');
ylim([0.2 1])
xlim([0 Tmax])
xlims = xlim;
ylims = ylim;
str = ['(c) $\Delta=$ ',num2str(Dlist(3))];
text((1-a)*xlims(1)+a*xlims(2), b*ylims(2)+(1-b)*ylims(1),str,'fontsize',myfont2,'Interpreter','Latex')

print(h1,'-depsc','perfect2.eps')

