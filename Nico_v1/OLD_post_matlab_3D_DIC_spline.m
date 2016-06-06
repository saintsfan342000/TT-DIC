%Nicolas Tardif
%01/07/2013 (dd,mm,yyyy)
%Post processing of the 3D DIC results when exported from a spline strain calculation using Nico_3D_if_spline 
clear
clf
%path for the files that the PERL program produces
PATH = 'C:\Students\Martin\Martin_Experiments\TT2-9_June27_a2p0\TT2-9_19-6_spline\post\'; %MUST CHANGE
%prefix of the name
prefix = 'TT2-9_19-6_spline'; %MUST CHANGE don't include underscore after linear
% # of the last picture
last = 467; %MUST CHANGE




%Here's the path for the perl file:
% perl "C:\Students\Martin\Martin_Experiments\DIC processing files\post_3D_DIC_spline.pl"

%Specimen dimensions in in                      %MUST CHANGE
Rm = 0.89340;                                    %MUST CHANGE
th = 0.0378;                                    %MUST CHANGE
Lg = 0.4;                                       %MUST CHANGE
S = 2 * pi * Rm * th;                           %MUST CHANGE
%Change mm->in
coef = 1/25.4;

%plotting parameters
%screen dimension
screen = 4/3;                                    %MUST CHANGE
%fontsize
fts = 18;
ftstics = round(fts*3/4);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%DIC paramaters
Facet_size = 19;%pix                            %MUST CHANGE                               %%VERY IMPORTANT MUST CHANGE
Step_size = 6;%pix                              %MUST CHANGE
%%%%%%%%%%%%%%!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%create the new folder where all the exported files are
folder = sprintf('export_FS_%d_SS_%d',Facet_size,Step_size);
mkdir(folder);
cd(folder);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%1-Determination of the 10 highest value of the cumulative plastic strain
%%%in the  the last picture.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%name of the file
name = sprintf('%s%s_%d.dat',PATH,prefix,last);
%open the file A = matrix (nx13) with
%A(:,1) index i [] 
%A(:,2) index j [] 
%A(:,3:5) undeformed coordinates X Y Z [mm] 
%A(:,6:8) displacement U V W [mm] 
%A(:,9) Equivalent plstic strain []
A = load(name);

%Coordinate of the center of the facets in the deformed configuration 
x=(A(:,3)+A(:,6))*coef/Rm;
y=(A(:,4)+A(:,7))*coef/Lg;
z=(A(:,5)+A(:,8))*coef/Lg;
%Calculation for each facet of the Logarithmic cumulative plastic strain
LEp=A(:,9);

%Creation of a grid in order to plot the data

    %Determination of the min distance between facet size
min_disp = 1000;
for i = [2 : round(size(x,1)/30) : size(x,1)-1]
    for j= [1:1:i-1,i+1:1:size(x,1)]
        min_disp = min ( min_disp , ( (x(i)-x(j))^2 + (y(i)-y(j) )^2 + ( z(i)-z(j) )^2 )^0.5 );
    end
end

 %creation of the grid [x,y]
dx=min_disp/5;
dy=min_disp/5;
%dz=min_disp/10;

x_edge=[min(x):dx:max(x)];
y_edge=[min(y):dy:max(y)];
%z_edge=[floor(min(z)):dz:ceil(max(z))];

[X,Y]=meshgrid(x_edge,y_edge);

agree = 0;

while agree == 0

%Find the 10 biggest values
[value,index_sort]=sort(-1*LEp);
index_crop = index_sort(1:10);


%projection of the cumulative plastic strain on the grid
LEP=griddata(x,y,LEp,X,Y);

close all;
%plot the data
 surf(X,Y,LEP,'EdgeColor','none')
 hold on
 for i= 1:1:10     
 plot3(x(index_crop(i)),y(index_crop(i)),LEp(index_crop(i)),'marker','o')
 end
xlabel('X coordinate / Rm []','Fontsize',fts)
ylabel('Y coordinate / Lg [] (axial)','Fontsize',fts)
zlabel('e^p','Fontsize',fts)
axis([min(x) max(x) min(y) max(y) min(LEp) max(LEp) min(LEp) max(LEp)])
colorbar('location','NorthOutside','Fontsize',fts) 
set (gca,'Xgrid','on','Ygrid','on','LineWidth',2)
set (gca,'Fontsize',ftstics)
set(gcf,'Units','Normalized','Outerposition',[0 0 1 1]);
hold off

figure
contourf(X,Y,LEP,10)
 hold on
 for i= 1:1:10     
 plot3(x(index_crop(i)),y(index_crop(i)),LEp(index_crop(i)),'marker','o','MarkerEdgeColor','y','MarkerFaceColor','y')
 end
xlabel('X coordinate / Rm []','Fontsize',fts)
ylabel('Y coordinate / Lg [] (axial)','Fontsize',fts)
set (gca,'Fontsize',ftstics)
set(gcf,'Units','Normalized','Outerposition',[0 0 1 1]); 
hold off

agree = input('do you agree with the selected points 0->no / 1->yes');

if agree == 0

fprintf('click on 2 points to define the crop zone');    
    
crop = ginput(2);
Xmin_crop = min(crop(:,1)); 
Xmax_crop = max(crop(:,1));
Ymin_crop = min(crop(:,2));
Ymax_crop = max(crop(:,2));

cropX = find(x < Xmin_crop | x > Xmax_crop);
cropY = find(y < Ymin_crop | y > Ymax_crop);
cropi = union(cropX,cropY);
x1 = x(cropi);
y1 = y(cropi);
LEp1 = LEp(cropi);
clear x y LEp;
x = x1;
y = y1;
LEp = LEp1;
clear x1 y1 LEp1;
end


end

close all;
%plot the data
 surf(X,Y,LEP,'EdgeColor','none')
 hold on
 for i= 1:1:10     
 plot3(x(index_crop(i)),y(index_crop(i)),LEp(index_crop(i)),'marker','o')
 end
xlabel('X coordinate / Rm []','Fontsize',fts)
ylabel('Y coordinate / Lg [] (axial)','Fontsize',fts)
zlabel('e^p','Fontsize',fts)
axis([min(x) max(x) min(y) max(y) min(LEp) max(LEp) min(LEp) max(LEp)])
colorbar('location','NorthOutside','Fontsize',fts) 
set (gca,'Xgrid','on','Ygrid','on','LineWidth',2)
set (gca,'Fontsize',ftstics)
set(gcf,'Units','Normalized','Outerposition',[0 0 1 1]);
hold off
print(gcf,'-dpng','last_picture.png');



 %Index of the facets which corresponds to the 10 highest cumulative
 %plastic strain in the last picture.
 index(:,1) = A(index_crop(:),1);
 index(:,2) = A(index_crop(:),2);
 
 %Determination of the initial step size and facet size compared to the thickness of the specimen
 %Coordinate of the center of the facets in the undeformed configuration 
X=(A(:,3))*coef/th;
Y=(A(:,4))*coef/th;
Z=(A(:,5))*coef/th;
 %Determination of the min distance between facet
min_disp = 1000;
for i = [2 : round(size(X,1)/30) : size(X,1)-1]
    for j= [1:1:i-1,i+1:1:size(X,1)]
        min_disp = min ( min_disp , ( (X(i)-X(j))^2 + (Y(i)-Y(j) )^2 + ( Z(i)-Z(j) )^2 )^0.5 );
    end
end
Step_size_th = min_disp;
Facet_size_th = Facet_size / Step_size * min_disp;
 
 
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%2-For the ten facets calculate all the technical strain and the 
%%%logarithmic plastic strain as a fonction of the time  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%time and Force data
name_force = sprintf('%stime_force.dat',PATH);
Force = load(name_force);

%Find the value of the index of the maximum force
[max_force,ind_max_force] = max(Force(:,3));
nb_profile = 10;
incr_Force = round((last-ind_max_force) / (nb_profile-1));
data_Force_point = [ind_max_force:incr_Force:ind_max_force+(nb_profile-2)*incr_Force,last]';
profile=1;
profile2=1;

%last pictures
last_pictures = [last-4:1:last]';

export_max=[];
export_mean=[];
export_stdv=[];


for k = 1 : last
    clear A;
    clear LEp;
    %name of the file
    name = sprintf('%s%s_%d.dat',PATH,prefix,k);
    %open the file A = matrix (nx13) with
    %A(:,1) index i [] 
    %A(:,2) index j [] 
    %A(:,3:5) undeformed coordinates X Y Z [mm] 
    %A(:,6:8) displacement U V W [mm] 
    %A(:,9) e^p []
    A = load(name);
    %Find the facets
    top10 = [];
    for i = 1 : size(index,1)
    top10 = [top10;find( A(:,1)==index(i,1)  & A(:,2)==index(i,2) )];
    end
    for i = 1 : size(top10)
    LEp(i) = A(top10(i),9);
    end
    export_max = [export_max;[Force(k+1,2), Force(k+1,3)/S/1000,LEp(1)]];
    export_mean = [export_mean;[Force(k+1,2), Force(k+1,3)/S/1000,nanmean(LEp)]];
    export_stdv = [export_stdv;[Force(k+1,2), Force(k+1,3)/S/1000,std(LEp)]];
    %Profile of the cumulative plastic strain along y for data_force_points
    %set
    if any(data_Force_point==k) == true
      index_ep = find(A(:,1)==index(1,1));
        for l = 1 : size(index_ep,1)
      %calculate Y coordinate
      pro_LEp(l,1,profile) = (A(index_ep(l),4)+A(index_ep(l),7))*coef/th;     
      %Logarithmic Cumulative Plastic strain
      pro_LEp(l,2,profile) = A(index_ep(l),9) ;      
        end
      profile=profile+1;
    end  
    %Profile of the cumulative plastic strain along y for last_pictures set
      if any(last_pictures==k) == true
      index_ep2 = find(A(:,1)==index(1,1));
        for l = 1 : size(index_ep2,1)
      %calculate Y coordinate
      pro_LEp2(l,1,profile2) = (A(index_ep2(l),4)+A(index_ep2(l),7))*coef/th;     
      %Logarithmic Cumulative Plastic strain
      pro_LEp2(l,2,profile2) =  A(index_ep2(l),9) ;       
        end
      profile2=profile2+1;
    end   
end        



%Plot time-cumulative plastic strain
figure
errorbar(export_mean(:,1),export_mean(:,3),1.5*export_stdv(:,3));
hold on
plot(export_max(:,1),export_max(:,3),'ro')
legend('mean value over 10','max value','Location','NorthWest')
axis([0 max(export_max(:,1)) 0 1.1*max(export_max(:,3))])
xlabel('time [s]','Fontsize',fts)
ylabel('cumulative plastic strain e^p []','Fontsize',fts)
set (gca,'Xgrid','on','Ygrid','on','LineWidth',2)
set (gca,'Fontsize',ftstics)
set(gcf,'Units','Normalized','Outerposition',[0 0 1/screen 1]);
hold off
print(gcf,'-dpng','time_ep.png');

%plot y,e^p for data_force_points set
figure
set(gcf,'Units','Normalized','Outerposition',[0 0 1 1]);

marker_type = ['+','o','*','.','x']';
marker_color = ['r','g','b','c','m','y','k']';
marker_style = [];
n = size(marker_type,1);
m = size(marker_color,1);
m0=1;
n0=1;
i=1;
while i <= nb_profile
    marker_style = [marker_style ; marker_color(m0) marker_type(n0)]; 
    i=i+1;
    m0 = m0+1;
    while m0 > m
        m0=1;
        n0=n0+1;
    end
end
marker_style = cellstr(marker_style);

subplot(1,2,1);
plot(Force(:,2),Force(:,3)/S/1000,'LineWidth',2);
hold on
for i=1:nb_profile
plot(Force(data_Force_point(i),2),Force(data_Force_point(i),3)/S/1000,char(marker_style(i)),'MarkerSize',10);
end
axis([0 max(Force(:,2)) 0 1.1*max(Force(:,3))/S/1000])
xlabel('time [s]','Fontsize',fts)
ylabel('\sigma_Y [ksi]','Fontsize',fts)
set (gca,'Xgrid','on','Ygrid','on','LineWidth',2)
set (gca,'Fontsize',ftstics)
hold off

subplot(1,2,2);
plot(pro_LEp(:,2,1),pro_LEp(:,1,1),char(marker_style(1)));
hold on
for i = 2 : nb_profile
plot(pro_LEp(:,2,i),pro_LEp(:,1,i),char(marker_style(i)));
end
xlabel('cumulative plastic strain e^p []','Fontsize',fts)
ylabel('axial coordinate Y/th []','Fontsize',fts)
set (gca,'Xgrid','on','Ygrid','on','LineWidth',2)
set (gca,'Fontsize',ftstics)
hold off
print(gcf,'-dpng','ep_profile')

%plot y,e^p for last pictures set
figure
set(gcf,'Units','Normalized','Outerposition',[0 0 1 1]);


marker_style = [];
m0=1;
n0=1;
i=1;
while i <= size(last_pictures,1)
    marker_style = [marker_style ; '-' marker_color(m0) marker_type(n0)]; 
    i=i+1;
    m0 = m0+1;
    while m0 > m
        m0=1;
        n0=n0+1;
    end
end
marker_style = cellstr(marker_style)

plot(pro_LEp2(:,2,1),pro_LEp2(:,1,1),char(marker_style(1)));
hold on
for i = 2 : size(last_pictures,1)
plot(pro_LEp2(:,2,i),pro_LEp2(:,1,i),char(marker_style(i)));
end
xlabel('cumulative plastic strain e^p []','Fontsize',fts)
ylabel('axial coordinate Y/th []','Fontsize',fts)
title('Max value of ep in Five last pictures','Fontsize',fts)
legend('1','2','3','4','last','Location','NorthWest')
set (gca,'Xgrid','on','Ygrid','on','LineWidth',2)
set (gca,'Fontsize',ftstics)
hold off
print(gcf,'-dpng','ep2_profile')
 

save('max.dat','export_max','-ASCII');
save('mean.dat','export_mean','-ASCII');
save('stdv.dat','export_stdv','-ASCII');

for i=1:nb_profile
  name = sprintf('export_Y_LEp_time_%6.3f_stress_%6.3f.dat',Force(data_Force_point(i),2),Force(data_Force_point(i),3)/S/1000);
  out = pro_LEp(:,:,i);
  save(name,'out','-ASCII');
end

fileID = fopen('facet_step.txt','w');
fprintf(fileID,'Facet size/th = %1.4f \n Step size/th = %1.4f',Facet_size_th,Step_size_th);
fclose(fileID);
    


% Must update path in lines 41, 143, 172