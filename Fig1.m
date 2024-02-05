% to prepare the fourier transform of sinx/x 
clear all; close all; clc; tic; myfont = 22;
mywidth = 2;

a1 = 0.3;

b1 = 0.25;

c1 = 0.25;

h1 = figure;

x1 = -a1:0.01:(1+a1);
y1 = zeros(1,length(x1));

x2 = x1+ 1 + 2*a1+ b1;
y2 = zeros(1,length(x2));

plot(x1, y1 ,'linewidth',mywidth)

hold on

plot(x2, y2,'linewidth',mywidth)

y3 = 0:0.01:1.7;
x3 = 0.5+ zeros(1,length(y3));

y4 = 0:0.01:1.7;
x4 = 0.5+zeros(1,length(y4)) +1 + 2*a1+ b1;

plot(x3, y3,'linewidth',mywidth)

plot(x4, y4,'linewidth',mywidth)

x5 = 0:0.01:1;
y5 = ones(1, length(x5));
plot(x5, y5, 'linewidth',mywidth)


y6 = 0:0.01:1;
x6 = zeros(1, length(y5));
plot(x6, y6, ':','linewidth',mywidth)

y7 = 0:0.01:1;
x7 = 1+zeros(1, length(y5));
plot(x7, y7, ':','linewidth',mywidth)

x8 = (1+2*a1+b1):0.01:(1.5+2*a1+b1);
y8 = 0.02*(0:1:50);
plot(x8, y8 ,'linewidth',mywidth)

x9 = 2*(1.5+2*a1+b1)- x8;
y9 = 0.02*(0:1:50);
plot(x9, y9 ,'linewidth',mywidth)

arrow([0.5 0],[1+a1 0],'width',2,'TipAngle',24,'BaseAngle',30,'FaceColor','b','EdgeColor','b')
arrow([2+2*a1+b1 0 ],[ 2+3*a1+b1 0],'width',1.5,'TipAngle',24,'BaseAngle',30,'FaceColor','b','EdgeColor','b')


arrow([0.5 0],[0.5 1.7],'width',1.5,'TipAngle',24,'BaseAngle',30,'FaceColor','b','EdgeColor','b')
arrow([1.5+2*a1+b1 0 ],[ 1.5+2*a1+b1 1.7],'width',1.5,'TipAngle',24,'BaseAngle',30,'FaceColor','b','EdgeColor','b')


text(1+0.05,1,'$ \pi $','fontsize',25,'Interpreter','latex')

text(1.5+2*a1+b1+0.05,1,'$ \pi $','fontsize',25,'Interpreter','latex')


text(-a1,1.6,'(a)','fontsize',22,'Interpreter','latex')
text(1+a1+b1,1.6,'(b)','fontsize',22,'Interpreter','latex')


text(0.59,1.6,'$F_1$','fontsize',25,'Interpreter','latex')
text(1.5+2*a1+b1+0.09,1.6,'$F_2$','fontsize',25,'Interpreter','latex')


text(1+0.75*a1,-0.75*c1,'$q$','fontsize',25,'Interpreter','latex')
text(2+2*a1+b1+0.75*a1,-0.75*c1,'$q$','fontsize',25,'Interpreter','latex')


text(0,-0.75*c1,'$-1$','HorizontalAlignment','center','fontsize',25,'Interpreter','latex')
text(1,-0.75*c1,'$+1$','HorizontalAlignment','center','fontsize',25,'Interpreter','latex')
text(0.5,-0.75*c1,'$0$','HorizontalAlignment','center','fontsize',25,'Interpreter','latex')

text(1+2*a1+b1,-0.75*c1,'$-2$','HorizontalAlignment','center','fontsize',25,'Interpreter','latex')
text(2+2*a1+b1,-0.75*c1,'$+2$','HorizontalAlignment','center','fontsize',25,'Interpreter','latex')
text(1.5+2*a1+b1,-0.75*c1,'$0$','HorizontalAlignment','center','fontsize',25,'Interpreter','latex')

text(2+2*a1+b1+0.75*a1,-0.75*c1,'$q$','fontsize',25,'Interpreter','latex')

axis equal 


axis off

print(h1,'-depsc','fourier.eps')
