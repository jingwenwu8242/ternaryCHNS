%轴对称图形的画法
clear; clc; clf;
close all;
ss1=sprintf('./phi1.m');
ss2=sprintf('./phi3.m');
% ss=sprintf('data/data/phi.m');
ssu=sprintf('./u.m');
ssw=sprintf('./w.m');
CC1=load(ss1);
CC2=load(ss2);
CC3=load(ssu);
CC4=load(ssw);
nxy=1;
ny=size(CC1,2);
nx=ny/nxy;
nn=size(CC1,1)/nx;
n=nx;
N=size(CC1,1)/nx;
xleft=0;
xright=1;
yleft=xleft;
yright=ny/nx*xright;
x=linspace(xleft,xright*2,nx*2);
y=linspace(yleft,yright,ny);
[xx,yy]=meshgrid(x,y);


aa=1:10:101;
for ij=1:size(aa,2)
    kk=aa(ij);
    figure(kk)
    AAA1=CC1(1+(kk-1)*n:kk*n,:);
    AAA2=CC2(1+(kk-1)*n:kk*n,:);
    AAA3=CC3(1+(kk-1)*n:kk*n,:);
    AAA4=CC4(1+(kk-1)*n:kk*n,:);
    BBB1=cat(2,0*AAA1',AAA1')';
    BBB2=cat(2,0*AAA2',AAA2')';
    BBB3=cat(2,0*AAA3',AAA3')';
    BBB4=cat(2,0*AAA4',AAA4')';
    for i=1:nx
        BBB1(i,:)=AAA1(nx+1-i,:);
        BBB2(i,:)=AAA2(nx+1-i,:);
        BBB3(i,:)=AAA3(nx+1-i,:);
        BBB4(i,:)=AAA4(nx+1-i,:);
    end
  
    hh=contourf(yy',xx',BBB1,[0.5 0.5],'facecolor','None','edgecolor','r','linewidth',2,'linestyle','-');hold on;
    contourf(yy',xx',BBB2,[0.5 0.5],'facecolor','None','edgecolor','k','linewidth',2,'linestyle','-');hold on;


    XX=hh(1,2:end);
    number=size(XX,2);
    XX(number+1)=XX(1);
    YY=hh(2,2:end);
    XX(number+1)=XX(1);
    YY(number+1)=YY(1);
    masspolygon(ij)=0;
    for i =1:number

        masspolygon(ij)= masspolygon(ij)+0.5*(XX(i)*YY(i+1)-YY(i)*XX(i+1));

    end
    colormap gray
    axis image
    axis([0 1 0.0 2.0 ])
    view(-90,90)

    set(gca,'fontsize',20)
    ss=sprintf('figa_%d.png',ij);
    print('-dpng',ss)


end

ss=sprintf('./mass1.m');
massphi1=load(ss);
a= massphi1(1);
fig2=figure(137);
hold on;
number=size(massphi1,1);
t = linspace(0,1,number);
massphi1=massphi1/a;
plot(t,massphi1(1:end),'ro-','linewidth',2);hold on;
box on;

ss1=sprintf('mass.png');
print(fig2,'-dpng',ss1)

fig2=figure(187);
b=masspolygon(1);
masspolygon=masspolygon/b;
number=size(masspolygon,2);
t = linspace(0,1,number);
plot(t,masspolygon(1:end),'ro-','linewidth',2);hold on;
box on;
ss1=sprintf('masspolygon.png');
print(fig2,'-dpng',ss1);