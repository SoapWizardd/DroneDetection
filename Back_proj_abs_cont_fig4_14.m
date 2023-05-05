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

%xsource=1;
%ysource=4.5;
sources=[ 6.5,4];

source_width_x=0.1;
source_width_y=0.1;
ywallbottom=0;

%ywalltop=1.5;
ywalltop=0.01;

xwall=5;
imovie=2;
ipauze=100;
nt=800;


T1=300;
stappen=0;

rij_sum=zeros(nx,ny,nt);
sumsum=zeros(nx,ny);


%yymin=-0.05;
%yymax=0.08;

yymin=-0.15;
yymax=0.15;

schaal_v=[0 1.5];
schaal_p=[-0.3 0.3];

%---FIGUUR---
%mics=[ 5,1,
 %      5,2,
  %     5,3,
   %    5,4];
%------------

%aantal_mics=5;
%mics=round(8*rand(aantal_mics,2)+1);

%mics=[ 9.0,6.5,
 %      .5,8.5,
  %     4.5,1.5,
   %    7.5,2.5,
    %   2.5,5.0,
     %  1.0,1.0,
      % 6.0,4.0];
      
      
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
%p=exp(-((xmat-xsource)/source_width_x).^2-((ymat-ysource)/source_width_x).^2);
%p=zeros(nx,ny);
%p_back=zeros(nx,ny);
%for i=1:size(sources,1),
 %   p=p+exp(-((xmat-sources(i,1))/0.1).^2-((ymat-sources(i,2))/0.1).^2);
%end

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



noice=randn(nt,size(sources,1));
for ii=1:nt
 voldx=vx;
 voldy=vy;
 
 p0=zeros(nx,ny);
 %%%Tijdsafhankelijke ingang%%%
 %if mod(ii,10)==1
 for i=1:size(sources,1)
     p0=p0+noice(ii,i)*exp(-((xmat-sources(i,2))/source_width_x).^2-((ymat-sources(i,1))/source_width_y).^2);
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
 
 
for m=1:size(mics,1)      
    mic(ii,m)=p(round((nx/xmax)*mics(m,2)),round((ny/ymax)*mics(m,1)));
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
     %pause(0.5);
     end;%end if imovie==1,
   


   
  end;%end if mod
 end;%end for ii
 
 
 
 sum=zeros(nx,ny,nt);
 for m=1:size(mics,1)
        for i=1:nx
            for j=1:ny
                d_m=(sqrt(((nx/xmax)*(mics(m,2))-i)^2 + (((ny/ymax)*mics(m,1))-j)^2));
                verschoven=circshift(mic(:,m),nt-round(d_m/(c0*(nx/xmax)*dt))).*sqrt(d_m);
                for ii=1:nt
                    sum(i,j,ii)=sum(i,j,ii)+verschoven(ii);
                end
            end
        end
 end
  
image_counter1=0;
AA1=clock;AAA1=['y',int2str(AA(1)),'m',int2str(AA(2)),'d',int2str(AA(3)),'_',int2str(AA(4)),'h',int2str(AA(5))];

 
 for ii=1:nt

    tijd=ii*dt;
      maximum = max(max(sum(:,:,ii)));
      [x,y]=find(sum(:,:,ii)==maximum);

    if mod(ii,ipauze)==1 || ii == nt
        figure(5);
   imagesc(sum(:,:,ii),[-10,10]);
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
    hold on
    plot(y,x,'kx','MarkerSize', 30)
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
   
    end;%end if mod
    
    
end
 
 
 gemiddelde=zeros(nx,ny);
 for i=1:nx
     for j=1:ny
         for ii=1:nt
            gemiddelde(i,j)=gemiddelde(i,j)+sum(i,j,ii)^2;
         end
     end
 end
 
  maximum = max(max(gemiddelde));
 [x,y]=find(gemiddelde==maximum);
 
 
 figure(1)
   imagesc(gemiddelde);
   title({'Sum of squares'},'FontSize',40);
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
    hold on
    plot(y,x,'kx','MarkerSize', 30)
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
   a.Label.Position(1) = 6.5;
   %pause(0.3);
   %pause
   
   
     hh1=figure(1);
     %set(hh,'position',[100 50 1400 700]);
     image_counter1=1;
     fprintf('storing image %i\n',image_counter1);
     G=getframe(hh1);
     %FF(image_counter)=G;
     [imind,cm] = rgb2ind(frame2im(G),256);
     if image_counter1==1, 
       imwrite(imind,cm,['su',AAA1,'.gif'],'gif', 'Loopcount',inf); 
       else 
       imwrite(imind,cm,['su',AAA1,'.gif'],'gif','WriteMode','append'); 
       end;%end if image_counter==1
     %pause(0.5);