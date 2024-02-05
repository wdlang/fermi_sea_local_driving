% tight binding model; defect mode to continuum
% exact results compared with 1st order perburbation theory
% 2014. April 7 in dresden
clear all; close all; clc; tic; format long; myfont = 18;

defect = 1;
% Tmax =780;
% L = 250;                               %  length of the arms
% Nspin = 2*L+1;                         %  number of spins
% 
% TTlist = [2.5, 2.504, 2.508];
% Dlist  = [0.02, 0.1];

Tmax =380;
L = 100;                               %  length of the arms
Nspin = 2*L+1;                         %  number of spins

TTlist = [2.5,2.51,2.52];
Dlist  = [0.05, 0.15];

plist = zeros(length(Dlist),length(TTlist),Tmax+1);
plist2 = zeros(length(Dlist),length(TTlist),Tmax+1);
plist3 = zeros(length(Dlist),length(TTlist),Tmax+1);

for sssss1 = 1: length(Dlist)
    for sssss2 = 1: length(TTlist)
        
        Delta =Dlist(sssss1);
        T=TTlist(sssss2);
        step = 30;
        omega = 2*pi/T;

        vec = zeros(Nspin, 1);
        lambda = (sqrt(defect^2+4)-defect)/2;
        for s = 1: Nspin
            vec(s) = lambda^(abs(s-L-1));
        end
        vec = vec./norm(vec);
        vec0 = vec';
        
        Ei = -sqrt(4+defect^2);
        Ef = Ei + omega;
        
        A0 = zeros(Nspin, Nspin);
        for s=1:Nspin-1
            A0(s,s+1)= -1;
            A0(s+1,s)= -1;
        end
        A0(L+1,L+1)= - defect;
        A0(Nspin,1)= -1;
        A0(1,Nspin)= -1;
        
        U = eye(Nspin, Nspin);
        for s = 1:step
            A1 = A0;
            A1(L+1, L+1) = -defect+ Delta*(cos(omega*(T/step)*s)- cos(omega*(T/step)*(s-1)))/(omega*(T/step));
            [V1, D1] = eig(A1);
            U1 = V1*diag(exp(-i*(T/step)*diag(D1)))*V1';
            U = U1*U;
        end
        
        for s = 0:Tmax
            plist (sssss1, sssss2,s+1) = abs(diag(vec0*vec)).^2;
            vec = U*vec;
        end
        
        [V, D ] = eig(A0);
        ddd = diag(D);
        [ddd, order] = sort(ddd);
        V = V(:, order);
        
        Veven = zeros(Nspin, Nspin);
        dddeven = zeros(1, Nspin);
        dimeven = 0;
        for s = 2: Nspin
            vec1 = V(:,s);
            if vec1'*flipud(vec1) > 0
                dimeven = dimeven+1;
                Veven(:,dimeven) = vec1;
                dddeven(dimeven) = ddd(s);
            end
        end
        coupling = abs(V(L+1, 1)*Veven(L+1, 1:dimeven));
        dddeven = dddeven(1:dimeven);
        
        index1 = find(Ef>dddeven, 1,'last');
        index2 = find(Ef<dddeven, 1,'first');
        
        E1 = dddeven(index1);
        E2 = dddeven(index2);
        gap = E2 - E1;
        if abs(E1-Ef)< abs(E2- Ef)
            couple = coupling (index1);
        else
            couple = coupling (index2);
        end
        
        eklist2 = E1 + gap*(-Nspin:(Nspin-1));
        
        for s=0:Tmax
            t = T*s;
            
            plist2(sssss1, sssss2, s+1) = 1- 4*((Delta/2)*couple)^2*sum( ((sin((eklist2-Ef).*(t/2)).^2)./((eklist2-Ef).^2)));
            plist3(sssss1, sssss2, s+1) = 1- 4*sum( (((Delta/2)*coupling).^2).* (((sin((dddeven-Ef).*(t/2)).^2)./((dddeven-Ef).^2))));
            
        end
        
    end
end

Tc = 2*pi/gap/T

hh = figure;
axes('Position',[0.20 0.59 0.6 0.36])
hold on
plot(0:Tmax, reshape(plist(1,:,:), length(TTlist), Tmax+1),'linewidth',1.0)
plot(0:Tmax, reshape(plist2(1,:,:), length(TTlist), Tmax+1),'--','linewidth',1.0)
ylim([0.8 1])
xlim([0 Tmax])
set(gca,'XTicklabel',[])
set(gca,'fontsize',myfont)
ylabel('$p$','fontsize',myfont,'Interpreter','Latex')
box on
xlims = xlim;
ylims = ylim;
a = 0.02; b = 0.1;
str = ['(a) $\Delta=$ ',num2str(Dlist(1))];
text((1-a)*xlims(1)+a*xlims(2), b*ylims(2)+(1-b)*ylims(1),str,'fontsize',myfont,'Interpreter','Latex')
set(gca,'ytick',[0.8 0.9  1] )
set(gca,'yticklabel',{'0.8','0.9', '1'})

axes('Position',[0.20 0.19 0.6 0.36])
hold on
plot(0:Tmax, reshape(plist(2,:,:), length(TTlist), Tmax+1),'linewidth',1.0)
plot(0:Tmax, reshape(plist2(2,:,:), length(TTlist), Tmax+1),'--','linewidth',1.0)
ylim([0 1])
xlim([0 Tmax])
set(gca,'fontsize',myfont)
ylabel('$p$','fontsize',myfont,'Interpreter','Latex')
xlabel('$\omega t /2\pi $','fontsize',myfont,'Interpreter','Latex')
box on
xlims = xlim;
ylims = ylim;
str = ['(b) $\Delta=$ ',num2str(Dlist(2))];
text((1-a)*xlims(1)+a*xlims(2), b*ylims(2)+(1-b)*ylims(1),str,'fontsize',myfont,'Interpreter','Latex')

print(hh,'-depsc','local.eps')
