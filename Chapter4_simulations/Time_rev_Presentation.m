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

%---FIGUUR---
sources=[ 6.5,4];
%------------

%---FIGUUR---
%sources=[1.5,1.5];
%------------


source_width_x=0.1;
source_width_y=0.1;
ywallbottom=0;

%ywalltop=1.5;
ywalltop=0.01;

xwall=5;
imovie=2;
ipauze=50;
nt=500;

%yymin=-0.05;
%yymax=0.08;

yymin=-0.15;
yymax=0.15;

schaal_v=[0 1.5];
schaal_p=[-0.1 0.1];

%%%Microfoons op cirkel
aantal=5;
radius = 4;
xCenter = 5;
yCenter = 5;
angles = linspace(0, 2*pi, aantal+1);
mics=zeros(aantal,2);
for i=1:aantal
    mics(i,1)=round(radius * cos(angles(i)) + xCenter);
    mics(i,2)=round(radius * sin(angles(i)) + yCenter);
end      

%mics=[1,5,
 %     3,8,
  %    3,2,
   %   7,8,
    %  7,2,
     % 9,5]
      

%---FIGUUR---
%aantal_mics=5; %Not to much, otherwise infinite loop
%mics=round(8*rand(1,2)+1);
%while size(mics,1) < aantal_mics
 %   random=round(8*rand(1,2)+1);
  %  tel=0;
   % for i=1:size(mics,1)
    %    if random ~= mics(i,:)
     %       tel=tel+1;
      %  end
    %end
   % if tel==size(mics,1)
  %      mics=[mics;random];
 %   end
%end
        

%aantal_mics=5;    
%mics=round(8*rand(aantal_mics,2)+1);
%------------

%---FIGUUR---
%mics=[ 5,1,
 %      5,2,
  %     5,3,
   %    5,4];
%------------
   
%---FIGUUR---
%mics=[1,1,
 %     1,1.5,
  %    1.5,1,
   %   1.5,1.5];
%------------

%---FIGUUR---
%mics=[1,1,
 %     1,9,
  %    9,1,
   %   9,9];
%------------
  
%---FIGUUR---
%mics=[5,3,
 %     3,5,
  %    7,5];
%------------

%mics=[5,3];

%mics=[5,3,
 %     3,5];
 
 
%mics=[1,2,
 %     1,3];
 
%mics=[5,2,
 %     5,6,
  %    9,8];
   
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

for i=1:size(sources,1),
    p=p+exp(-((xmat-sources(i,2))/source_width_x).^2-((ymat-sources(i,1))/source_width_y).^2);
end

%p=exp(-((xmat-xsource)/source_width_x).^2-((ymat-ysource)/source_width_x).^2);


dx=x(2)-x(1);
dy=y(2)-y(1);

c=ones(nx,ny)*c0;

rho=ones(nx,ny)*rho0;

%rho(:,150:200)=0.1;

rhox=0.5*p./(c.^2);
rhoy=0.5*p./(c.^2);


dt=0.5*sqrt(dx^2+dy^2)/max(max(c));


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



%T=(1:nt)';
image_counter=0;
AA=clock;AAA=['y',int2str(AA(1)),'m',int2str(AA(2)),'d',int2str(AA(3)),'_',int2str(AA(4)),'h',int2str(AA(5))];

noice=rand(nt,1);
for ii=1:nt
 voldx=vx;
 voldy=vy;
 
 
 p0=zeros(nx,ny);
 %%%Tijdsafhankelijke ingang%%%
 %if mod(ii,10)==1
 
 %for i=1:size(sources,1),
  %   p0=p0+noice(ii)*exp(-((xmat-sources(i,1))/source_width_x).^2-((ymat-sources(i,2))/source_width_y).^2);
 %end
 
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
 %IYbottom=find(y>=ywallbottom);
 %IYtop=find(y>=ywalltop);
 %IXW=find(x>=xwall);
 %vy(IYbottom(1):IYtop(1),IXW(1))=0;
 
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
  
 if mod(ii,ipauze)==1 || ii == nt
   

   figure(4);
   tiledlayout(1,2)
   nexttile
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
   plot((nx/xmax)*sources(:,1),(ny/ymax)*sources(:,2),'ko','MarkerSize', 20)
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
   nexttile
 for i=1:size(mics,1)
     rec=mic(:,i);
 plot(rec(1:ii),'DisplayName',['(',num2str(mics(i,1)),'m,',num2str(mics(i,2)),'m)'])
   title({'Recordings at microphones'},'FontSize',40);
   xlabel('t (ms)','FontSize',30)
   ylabel("p' (Pa)",'FontSize',30)
   xlim([0 nt]);
   ylim([-0.04 0.06]);
 h = gca; 
h.XTickMode = 'manual'; 
NumTicks = 5;
L = get(gca,'XLim');
set(gca,'XTick',linspace(L(1),L(2),NumTicks))
h.XTickLabel = round((h.XTick * dt)*1000)/10 ; 
h.FontSize = 30;
 hold on
 end
 legend
 
   
   
   if imovie==2,
     hh=figure(4);
     %set(hh,'position',[100 50 1400 700]);
     image_counter=image_counter+1;
     fprintf('storing image %i\n',image_counter);
     G=getframe(hh);
     %FF(image_counter)=G;
     [imind,cm] = rgb2ind(frame2im(G),256);
     if image_counter==1, 
       imwrite(imind,cm,['simultion-diffraction-screen_',AAA,'.gif'],'gif', 'Loopcount',inf); 
       else 
       imwrite(imind,cm,['simultion-diffraction-screen_',AAA,'.gif'],'gif','WriteMode','append'); 
       end;%end if image_counter==1
     pause(0.5);
     end;%end if imovie==1,
   
   
  end;%end if mod
 end;%end for ii

 
  figure(1)
 tiledlayout(1,2)
   nexttile
  for i=1:size(mics,1)
 plot((mic(:,i)),'DisplayName',['(',num2str(mics(i,1)),'m,',num2str(mics(i,2)),'m)'])
   title({'Recordings at microphones'},'FontSize',40);
   xlabel('t (ms)','FontSize',30)
   ylabel("p' (Pa)",'FontSize',30)
   xlim([0 nt]);
   ylim([-0.04 0.06]);
 h = gca; 
h.XTickMode = 'manual'; 
NumTicks = 5;
L = get(gca,'XLim');
set(gca,'XTick',linspace(L(1),L(2),NumTicks))
h.XTickLabel = round((h.XTick * dt)*1000)/10 ; 
h.FontSize = 30;
 hold on
 end
 legend
 
 nexttile
  for i=1:size(mics,1)
 plot(flip(mic(:,i)),'DisplayName',['(',num2str(mics(i,1)),'m,',num2str(mics(i,2)),'m)'])
   title({'Time reversed recordings'},'FontSize',40);
   xlabel('t (ms)','FontSize',30)
   ylabel("p' (Pa)",'FontSize',30)
   xlim([0 nt]);
   ylim([-0.04 0.06]);
 h = gca; 
h.XTickMode = 'manual'; 
NumTicks = 5;
L = get(gca,'XLim');
set(gca,'XTick',linspace(L(1),L(2),NumTicks))
h.XTickLabel = round((h.XTick * dt)*1000)/10 ; 
h.FontSize = 30;
 hold on
 end
 legend('Location','northwest')
 
 pause(1);
 
 
 
 
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
ipauze=50;
nt=500;

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


dt=0.5*sqrt(dx^2+dy^2)/max(max(c));


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
 
 tijd=ii*dt;
  
 if mod(ii,ipauze)==1 || ii == nt
  
 %maximum = max(max(p));
 %[x,y]=find(p==maximum);
 regionalMaxima = imregionalmax(p);
 [rows, columns] = find(regionalMaxima);
   
  figure(4);
   tiledlayout(1,2)
   nexttile
    for i=1:size(mics,1)
        rev=flip(mic(:,i));
 plot(rev(1:ii),'DisplayName',['(',num2str(mics(i,1)),'m,',num2str(mics(i,2)),'m)'])
   title({'Time reversed recordings'},'FontSize',40);
   xlabel('t (ms)','FontSize',30)
   ylabel("p' (Pa)",'FontSize',30)
   xlim([0 nt]);
   ylim([-0.04 0.06]);
 h = gca; 
h.XTickMode = 'manual'; 
NumTicks = 5;
L = get(gca,'XLim');
set(gca,'XTick',linspace(L(1),L(2),NumTicks))
h.XTickLabel = round((h.XTick * dt)*1000)/10 ; 
h.FontSize = 30;
 hold on
    end
 legend('Location','northwest')
   
   nexttile
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
   plot((nx/xmax)*sources(:,1),(ny/ymax)*sources(:,2),'ko','MarkerSize', 20)
   hold on
   plot((nx/xmax)*mics(:,1),(ny/ymax)*mics(:,2),'rx','MarkerSize', 20)
   
   %if ii == nt
    %hold on
    %plot(y,x,'kx','MarkerSize', 30)
   %end
   maxim=zeros(size(rows,1));
   for r=1:size(rows,1)
       maxim(r)=p(rows(r),columns(r));
   end
   maxi=zeros(size(sources,1),1);
   for s=1:size(sources,1)
       maxi(s,1)=max(max(maxim));
       pos_max=find(maxim==maxi(s));
       maxim(pos_max)=0;
       [x,y]=find(p==maxi(s));
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
       imwrite(imind,cm,['simulation-diffraction-screen_',AAA,'.gif'],'gif', 'Loopcount',inf); 
       else 
       imwrite(imind,cm,['simulation-diffraction-screen_',AAA,'.gif'],'gif','WriteMode','append'); 
       end;%end if image_counter==1
     pause(0.5);
     end;%end if imovie==1,
   
   
  end;%end if mod
 end;%end for ii
 
% for i=1:size(mics,1)
 %figure(i+10)
 %plot(flip(mic(:,i)))
  % title({['Recording at: (',num2str(mics(i,1)),'m,',num2str(mics(i,2)),'m)']},'FontSize',40);
   %xlabel('t (ms)','FontSize',30)
%   ylabel("p' (Pa)",'FontSize',30)
 %  xlim([0 nt]);
  % ylim([-0.04 0.06]);
 %h = gca; 
%h.XTickMode = 'manual'; 
%NumTicks = 5;
%L = get(gca,'XLim');
%set(gca,'XTick',linspace(L(1),L(2),NumTicks))
%h.XTickLabel = round((h.XTick * dt)*1000)/10 ; 
%h.FontSize = 30;
 %end