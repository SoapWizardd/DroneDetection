%program finite_difference_calculation_2D_with_screen_v1_scaled.m
%first edition 29112017
%addition of absorbing boundary layers (PML) Lei 19-01-2018
%screen 14012018
clear all
%close all

ax=100;
ay=100;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%USER DEFINITIONS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c0=340;
rho0=1.3;
factor=4;
xmax=10;
ymax=10;
nx=ax*factor;
ny=ay*factor;

sources=[5,5];

source_width_x=0.1;
source_width_y=0.1;
ywallbottom=0;

%ywalltop=1.5;
ywalltop=0.01;

xwall=5;
imovie=1;
ipauze=50*5/2;
nt=2000;

%yymin=-0.05;
%yymax=0.08;

yymin=-0.15;
yymax=0.15;

schaal_v=[0 1.5];
schaal_p=[-0.1 0.1];




%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%



vx=zeros(nx,ny);
vy=zeros(nx,ny);

x=linspace(0,xmax,nx)';
y=linspace(0,ymax,ny)';
xmat=x*ones(1,ny);
ymat=ones(nx,1)*y';

p=zeros(nx,ny);
%%%Ingang t=0%%%

for i=1:size(sources,1),
    p=p+exp(-((xmat-sources(i,1))/source_width_x).^2-((ymat-sources(i,2))/source_width_y).^2);
end

%p=exp(-((xmat-xsource)/source_width_x).^2-((ymat-ysource)/source_width_x).^2);


dx=x(2)-x(1);
dy=y(2)-y(1);

c=ones(nx,ny)*c0;

rho=ones(nx,ny)*rho0;

%rho(:,150:200)=0.1;

rhox=0.5*p./(c.^2);
rhoy=0.5*p./(c.^2);


dt=0.2*sqrt(dx^2+dy^2)/max(max(c));



%T=(1:nt)';
image_counter=0;
AA=clock;AAA=['y',int2str(AA(1)),'m',int2str(AA(2)),'d',int2str(AA(3)),'_',int2str(AA(4)),'h',int2str(AA(5))];




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%staggered grid



for ii=1:nt;
 voldx=vx;
 voldy=vy;
 
 p0=zeros(nx,ny);
 %%%Tijdsafhankelijke ingang%%%
 
 %for i=1:size(sources,1),
  %   p0=p0+abs(z(ii+1000000,1)+z(ii+1000000,2))*exp(-((xmat-sources(i,1))/0.1).^2-((ymat-sources(i,2))/0.1).^2);
 %end
 
 
 pold=p+p0;
 
 %pold=p;
 
 vx=voldx-1./rho.*(circshift(p,[1,0])-circshift(p,[0,0]))/(dx)*dt;
 vy=voldy-1./rho.*(circshift(p,[0 1])-circshift(p,[0 0]))/(dx)*dt;
 p=pold-rho.*c.^2.*((circshift(vx,[0 0])-circshift(vx,[-1 0]))/(dx)+(circshift(vy,[0 0])-circshift(vy,[0 -1]))/(dy))*dt;
 vx(1,:)=0;
 vx(end,:)=0;
 vy(:,1)=0;
 vy(:,end)=0;
 
 
 for i=1:size(mics,1),
     mic(ii,i)=p((nx/xmax)*mics(i,1),(ny/ymax)*mics(i,2));
 end
 
  tijd=ii*dt;
 
 if mod(ii,ipauze)==1 || ii == nt
   

  figure(3);
   imagesc(p,schaal_p);
   title({['Time (s): ',num2str(tijd)]},'FontSize',40);
   xlabel('x (m)','FontSize',30)
   ylabel('y (m)','FontSize',30)
   
   zoom on
   xlim([0 nx]);
   ylim([0 ny]);
   set(gca,'Ydir','normal');
   hold on
   %h3=plot([IXW(1),IXW(1)],[IYbottom(1),IYtop(1)],'r');
   %set(h3,'LineWidth',3);
   %hold on
   plot((nx/xmax)*sources(:,2),(ny/ymax)*sources(:,1),'ko','MarkerSize', 20)
   
   hold off
   axis equal
   axis on
   %set(cga,'XTickLabel',{0,1,2,3,4,5,6,7,8,9,10});
   xticks([linspace(0,nx,11)])
   xticklabels({'0','1','2','3','4','5','6','7','8','9','10'})
   yticks([linspace(0,nx,11)])
   yticklabels({'0','1','2','3','4','5','6','7','8','9','10'})
   ax = gca;
   ax.FontSize = 30; 
   a = colorbar;
   ylabel(a,"p' (Pa)",'FontSize',30,'Rotation',270);
   a.Label.Position(1) = 6;
   %pause(0.3);
   %pause
   
   
   
   if imovie==1,
     hh=figure(3);
     %set(hh,'position',[100 50 1400 700]);
     image_counter=image_counter+1;
     fprintf('storing image %i\n',image_counter);
     G=getframe(hh);
     %FF(image_counter)=G;
     [imind,cm] = rgb2ind(frame2im(G),256);
     if image_counter==1, 
       imwrite(imind,cm,['simulation-diffraction-screen_',AAA,'.gif'],'gif', 'Loopcount',inf); 
       else 
       imwrite(imind,cm,['simulation-diffraction-screen_',AAA,'.gif'],'gif','WriteMode','append'); 
       end;%end if image_counter==1
     pause(0.5);
     end;%end if imovie==1,
   
     
   end;%end if mod
 end;%end for ii

