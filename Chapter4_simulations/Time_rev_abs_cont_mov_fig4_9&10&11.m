%program finite_difference_calculation_2D_with_screen_v1_scaled.m
%first edition 29112017
%addition of absorbing boundary layers (PML) Lei 19-01-2018
%screen 14012018
clear all
%close all

ax=100;
ay=100;

[z,Fs] = audioread('drone_felix.mp3');
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

sources=[6.5,4];
afstand=2;

source_width_x=0.1;
source_width_y=0.1;
ywallbottom=0;

%ywalltop=1.5;
ywalltop=0.01;

xwall=5;
imovie=2;
ipauze=200;
nt=1000;
T0=10;

%yymin=-0.05;
%yymax=0.08;

yymin=-0.15;
yymax=0.15;

schaal_v=[0 1.5];
schaal_p=[-0.1 0.1];


%mics=[ 9.0,6.5,
 %      .5,8.5,
  %     4.5,1.5,
   %    7.5,2.5,
    %   2.5,5.0,
     %  1.0,1.0,
      % 6.0,4.0];

%%%Microfoons op cirkel
aantal=10;
radius = 4;
xCenter = 5;
yCenter = 5;
angles = linspace(0, 2*pi, aantal+1);
mics=zeros(aantal,2);
for i=1:aantal
    mics(i,1)=round(radius * cos(angles(i)) + xCenter);
    mics(i,2)=round(radius * sin(angles(i)) + yCenter);
end      
      
mic=zeros(nt,size(mics,1));

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

%for i=1:size(sources,1),
 %   p=p+exp(-((xmat-sources(i,1))/source_width_x).^2-((ymat-sources(i,2))/source_width_y).^2);
%end

%p=exp(-((xmat-xsource)/source_width_x).^2-((ymat-ysource)/source_width_x).^2);


dx=x(2)-x(1);
dy=y(2)-y(1);

c=ones(nx,ny)*c0;

rho=ones(nx,ny)*rho0;

%rho(:,150:200)=0.1;

rhox=0.5*p./(c.^2);
rhoy=0.5*p./(c.^2);


dt=0.4360*sqrt(dx^2+dy^2)/max(max(c));


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% PML
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pml_x=ones(nx,1);
pml_y=ones(1,ny);

pml_size=20;  %%in number of pixels
decay_factor=4;

pml_index=1:pml_size;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
pml_x_up=decay_factor*(max(c(:))/dx)*((pml_index-pml_size-1)./(0-pml_size)).^4;
pml_x_down=decay_factor*(max(c(:))/dx)*(pml_index./pml_size).^4;

pml_x_up=exp(-pml_x_up*dt/2);
pml_x_down=exp(-pml_x_down*dt/2);


pml_y_left=decay_factor*(max(c(:))/dy)*((pml_index-pml_size-1)./(0-pml_size)).^4;
pml_y_right=decay_factor*(max(c(:))/dy)*(pml_index./pml_size).^4;

pml_y_left=exp(-pml_y_left*dt/2);
pml_y_right=exp(-pml_y_right*dt/2);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
pml_x(1:pml_size,1)=pml_x_up;
pml_x(end-pml_size+1:end,1)=pml_x_down;
pml_y(1,1:pml_size)=pml_y_left;
pml_y(1,end-pml_size+1:end)=pml_y_right;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
decay_matrix_x=pml_x*ones(1,ny);
decay_matrix_y=ones(nx,1)*pml_y;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%staggered grid



T=(1:nt)';
image_counter=0;
AA=clock;AAA=['y',int2str(AA(1)),'m',int2str(AA(2)),'d',int2str(AA(3)),'_',int2str(AA(4)),'h',int2str(AA(5))];


noice=randn(nt,1);
for ii=1:nt
 voldx=vx;
 voldy=vy;
 
 
 p0=zeros(nx,ny);
 %%%Tijdsafhankelijke ingang%%%
 
 %if mod(ii,T0)==1 
    if ii==1
        stappen=0;
    else
        stappen=afstand*ii/nt;
    end
    
    for i=1:size(sources,1),
        p0=p0+z(ii+100000,1)*exp(-((xmat-sources(i,2)-stappen)/source_width_x).^2-((ymat-sources(i,1))/source_width_y).^2);
        %p0=p0+noice(ii)*exp(-((xmat-sources(i,2)-stappen)/source_width_x).^2-((ymat-sources(i,1))/source_width_y).^2);
    end
 %end
 
 
 %pold=p+p0;
 %pold=p;
 
 rhoxold=rhox+0.5*p0./(c.^2);
 rhoyold=rhoy+0.5*p0./(c.^2);
 
 %rhoxold=rhox;
 %rhoyold=rhoy;


 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 vx=decay_matrix_x.*(decay_matrix_x.*voldx-1./rho.*(circshift(p,[1,0])-circshift(p,[0,0]))/(dx)*dt);
 vy=decay_matrix_y.*(decay_matrix_y.*voldy-1./rho.*(circshift(p,[0 1])-circshift(p,[0 0]))/(dy)*dt);

 vx(1,:)=0;
 vx(end,:)=0;
 vy(:,1)=0;
 vy(:,end)=0;
 IYbottom=find(y>=ywallbottom);
 IYtop=find(y>=ywalltop);
 IXW=find(x>=xwall);
 vy(IYbottom(1):IYtop(1),IXW(1))=0;
 
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 rhox=decay_matrix_x.*(decay_matrix_x.*rhoxold-rho.*(circshift(vx,[0,0])-circshift(vx,[-1,0]))/(dx)*dt);
 rhoy=decay_matrix_y.*(decay_matrix_y.*rhoyold-rho.*(circshift(vy,[0,0])-circshift(vy,[0,-1]))/(dy)*dt);
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 
 %%%%%%%%
 p=c.^2.*(rhox+rhoy);
 %%%%%%%%
 
 
 for i=1:size(mics,1),
     mic(ii,i)=p((nx/xmax)*mics(i,2),(ny/ymax)*mics(i,1));
 end
 
 tijd=ii*dt;
  
 if mod(ii,ipauze)==1
   

   figure(4);
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
   plot((nx/xmax)*sources(:,1),(ny/ymax)*(sources(:,2)+(afstand*ii/nt)),'ko','MarkerSize', 20)
   hold on
   plot((nx/xmax)*mics(:,1),(ny/ymax)*mics(:,2),'rx','MarkerSize', 20)
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
   
   
   if imovie==2,
     hh=figure(4);
     %set(hh,'position',[100 50 1400 700]);
     image_counter=image_counter+1;
     fprintf('storing image %i\n',image_counter);
     G=getframe(hh);
     %FF(image_counter)=G;
     [imind,cm] = rgb2ind(frame2im(G),256);
     if image_counter==1, 
       imwrite(imind,cm,['simulation-diffraction-screen1_',AAA,'.gif'],'gif', 'Loopcount',inf); 
       else 
       imwrite(imind,cm,['simulation-diffraction-screen1_',AAA,'.gif'],'gif','WriteMode','append'); 
       end;%end if image_counter==1
     %pause(0.5);
     end;%end if imovie==1,
  
  end;%end if mod
 end;%end for ii

 
 
 
 
 
 
 
 
 
 %program finite_difference_calculation_2D_with_screen_v1_scaled.m
%first edition 29112017
%addition of absorbing boundary layers (PML) Lei 19-01-2018
%screen 14012018
%clear all
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
xsource=1;
ysource=4.5;
source_width_x=0.1;
source_width_y=0.1;
ywallbottom=0;

%ywalltop=1.5;
ywalltop=0.01;

xwall=5;
imovie=2;
ipauze=100;
nt=1000;

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
%p=exp(-((xmat-xsource)/source_width_x).^2-((ymat-ysource)/source_width_x).^2);


dx=x(2)-x(1);
dy=y(2)-y(1);

c=ones(nx,ny)*c0;

rho=ones(nx,ny)*rho0;

%rho(:,150:200)=0.1;

rhox=0.5*p./(c.^2);
rhoy=0.5*p./(c.^2);


dt=0.4360*sqrt(dx^2+dy^2)/max(max(c));


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% PML
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pml_x=ones(nx,1);
pml_y=ones(1,ny);

pml_size=20;  %%in number of pixels
decay_factor=4;

pml_index=1:pml_size;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
pml_x_up=decay_factor*(max(c(:))/dx)*((pml_index-pml_size-1)./(0-pml_size)).^4;
pml_x_down=decay_factor*(max(c(:))/dx)*(pml_index./pml_size).^4;

pml_x_up=exp(-pml_x_up*dt/2);
pml_x_down=exp(-pml_x_down*dt/2);


pml_y_left=decay_factor*(max(c(:))/dy)*((pml_index-pml_size-1)./(0-pml_size)).^4;
pml_y_right=decay_factor*(max(c(:))/dy)*(pml_index./pml_size).^4;

pml_y_left=exp(-pml_y_left*dt/2);
pml_y_right=exp(-pml_y_right*dt/2);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
pml_x(1:pml_size,1)=pml_x_up;
pml_x(end-pml_size+1:end,1)=pml_x_down;
pml_y(1,1:pml_size)=pml_y_left;
pml_y(1,end-pml_size+1:end)=pml_y_right;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
decay_matrix_x=pml_x*ones(1,ny);
decay_matrix_y=ones(nx,1)*pml_y;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%staggered grid



T=(1:nt)';
image_counter=0;
AA=clock;AAA=['y',int2str(AA(1)),'m',int2str(AA(2)),'d',int2str(AA(3)),'_',int2str(AA(4)),'h',int2str(AA(5))];

image_counter1=0;
AA1=clock;AAA1=['y',int2str(AA(1)),'m',int2str(AA(2)),'d',int2str(AA(3)),'_',int2str(AA(4)),'h',int2str(AA(5))];


sum2=zeros(nx,ny);


for ii=1:nt
 voldx=vx;
 voldy=vy;
 
 
 p0=zeros(nx,ny);
    
    %for i=1:size(mics,1),
     %   p0(round(mics(i,1)*nx/xmax),round(mics(i,2)*ny/ymax))=mic(nt-ii+1,i);
    %end
    
    for i=1:size(mics,1),
        p0=p0+mic(nt-ii+1,i)*exp(-((xmat-mics(i,2))/source_width_x).^2-((ymat-mics(i,1))/source_width_y).^2);
    end
    rhoxold=rhox+0.5*p0./(c.^2);
    rhoyold=rhoy+0.5*p0./(c.^2);
    %pold=p+p0;
 %pold=p;
  
 %rhoxold=0.5*pold./(c.^2);
 %rhoyold=0.5*pold./(c.^2);

 
 %rhoxold=rhox;
 %rhoyold=rhoy;


 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 vx=decay_matrix_x.*(decay_matrix_x.*voldx-1./rho.*(circshift(p,[1,0])-circshift(p,[0,0]))/(dx)*dt);
 vy=decay_matrix_y.*(decay_matrix_y.*voldy-1./rho.*(circshift(p,[0 1])-circshift(p,[0 0]))/(dy)*dt);

 vx(1,:)=0;
 vx(end,:)=0;
 vy(:,1)=0;
 vy(:,end)=0;
 IYbottom=find(y>=ywallbottom);
 IYtop=find(y>=ywalltop);
 IXW=find(x>=xwall);
 vy(IYbottom(1):IYtop(1),IXW(1))=0;
 
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 rhox=decay_matrix_x.*(decay_matrix_x.*rhoxold-rho.*(circshift(vx,[0,0])-circshift(vx,[-1,0]))/(dx)*dt);
 rhoy=decay_matrix_y.*(decay_matrix_y.*rhoyold-rho.*(circshift(vy,[0,0])-circshift(vy,[0,-1]))/(dy)*dt);
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 
 %%%%%%%%
 p=c.^2.*(rhox+rhoy);
 %%%%%%%%
 
 sum2 = sum2 + p.*p;
 
 tijd=ii*dt;
  
 if mod(ii,ipauze)==1 || ii == nt
   
   
     
   maximum = max(max(p));
   [x,y]=find(p==maximum);
   
   minimum = min(min(p));
   [a,b]=find(p==minimum);

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
   plot((nx/xmax)*sources(:,1),(ny/ymax)*(sources(:,2)+afstand-(afstand*ii/nt)),'ko','MarkerSize', 20)
   hold on
   plot((nx/xmax)*mics(:,1),(ny/ymax)*mics(:,2),'rx','MarkerSize', 20)
   %if ii == nt
    hold on
    plot(y,x,'kx','MarkerSize', 30)
    hold on
    plot(b,a,'wx','MarkerSize', 30)
   %end
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
   
    if imovie==2,
     hh=figure(3);
     %set(hh,'position',[100 50 1400 700]);
     image_counter=image_counter+1;
     fprintf('storing image %i\n',image_counter);
     G=getframe(hh);
     %FF(image_counter)=G;
     [imind,cm] = rgb2ind(frame2im(G),256);
     if image_counter==1, 
       imwrite(imind,cm,['p_',AAA,'.gif'],'gif', 'Loopcount',inf); 
       else 
       imwrite(imind,cm,['p_',AAA,'.gif'],'gif','WriteMode','append'); 
       end;%end if image_counter==1
     %pause(0.5);
     end;%end if imovie==1,
   
  
 %maximum = max(max(sum2));
 %[x,y]=find(sum2==maximum);
 regionalMaxima = imregionalmax(sum2);
 [rows, columns] = find(regionalMaxima);
   
  figure(5);
   imagesc(sum2);
   title({['SoS for times ','[',num2str(tijd-ipauze*dt), ',',num2str(tijd),'](s)']},'FontSize',40);
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
   plot((nx/xmax)*sources(:,1),(ny/ymax)*(sources(:,2)+afstand-(afstand*ii/nt)),'ko','MarkerSize', 20)
   hold on
   plot((nx/xmax)*mics(:,1),(ny/ymax)*mics(:,2),'rx','MarkerSize', 20)
   %if ii == nt
    %hold on
    %plot(y,x,'kx','MarkerSize', 30)
   %end
   maxim=zeros(size(rows,1));
   for r=1:size(rows,1)
       maxim(r)=sum2(rows(r),columns(r));
   end
   maxi=zeros(size(mics,1)+1,1);
   for s=1:size(maxi,1)
       maxi(s,1)=max(max(maxim));
       pos_max=find(maxim==maxi(s));
       maxim(pos_max)=0;
       [x,y]=find(sum2==maxi(s));
   %if ii == nt
    hold on
    plot(y,x,'kx','MarkerSize', 30)
   %end
   end
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
   a.Label.Position(2) = 0;
   %pause(0.3);
   %pause
   
    if imovie==2,
     hh1=figure(5);
     %set(hh,'position',[100 50 1400 700]);
     image_counter1=image_counter1+1;
     fprintf('storing image %i\n',image_counter1);
     G=getframe(hh1);
     %FF(image_counter)=G;
     [imind,cm] = rgb2ind(frame2im(G),256);
     if image_counter1==1, 
       imwrite(imind,cm,['sum_',AAA1,'.gif'],'gif', 'Loopcount',inf); 
       else 
       imwrite(imind,cm,['sum_',AAA1,'.gif'],'gif','WriteMode','append'); 
       end;%end if image_counter==1
     %pause(0.5);
     end;%end if imovie==1,
   
   sum2=zeros(nx,ny);
  end;%end if mod
 end;%end for ii

