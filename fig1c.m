% fig1c 
clear all; close all; clc; tic; myfont = 22;


Nlist = -1000:1000;

Tlist = 0.01:0.01:pi*4;

alist = [0.1, 0.2,0.3,0.5];

plist = zeros(length(alist), length(Tlist));

for s1 = 1: length(alist)
    a = alist(s1);
    for s2 =1 : length(Tlist)
        T = Tlist(s2);
        
        array = (Nlist+a)*T;
        
        plist(s1, s2) = (T^2)*sum( (sin(array)./array).^2   );
    end
end

h1= figure;
axes('Position',[0.15 0.2 0.8 0.6])
plot(Tlist./pi, plist(1,:)./pi^2, Tlist./pi, plist(2,:)./pi^2,':',Tlist./pi,...
    plist(3,:)./pi^2,'--',Tlist./pi, plist(4,:)./pi^2,'-.','linewidth',2);

set(gca,'fontsize',myfont);

ylim([0 5.5])
text(3.6,5.1,'(c)','fontsize',26,'Interpreter','latex')


xlabel('$T/\pi$','fontsize',myfont,'Interpreter','Latex')
ylabel('$W_\alpha(T)/\pi^2$','fontsize',myfont,'Interpreter','Latex')
str1 = strcat('$\alpha =',num2str(alist(1)),'$');
str2 = strcat('$\alpha=',num2str(alist(2)),'$');
str3 = strcat('$\alpha=',num2str(alist(3)),'$');
str4 = strcat('$\alpha=',num2str(alist(4)),'$');
hleg = legend(str1,str2, str3, str4);
set(hleg,'location','Northwest','box','off','Interpreter','Latex')

print(h1,'-depsc','Wa.eps')
