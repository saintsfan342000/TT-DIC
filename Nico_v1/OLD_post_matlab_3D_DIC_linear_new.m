function OLD_post_matlab_3D_DIC_linear_new();

%Nicolas Tardif
%01/07/2013 (dd,mm,yyyy)
%Post processing of the 3D DIC results when exported from a linear strain calculation using Nico_3D_linear3_fv 
clear
clf
%path for the files that the PERL program produces
    PATH = 'F:\Martin_Experiments\TT2-17\postOLDMETH\'; %MUST CHANGE
%prefix of the name
    prefix = 'TT2-17_OLD'; %MUST CHANGE don't include underscore after linear
%DIC paramaters
    Facet_size = 19;%pix                            %MUST CHANGE                               %%VERY IMPORTANT MUST CHANGE
    Step_size = 6;%pix                              %MUST CHANGE
% # of the last picture
    last = 609; %MUST CHANGE
    TT2=17;
    alpha=1;
    calctype='Linear';
%plotting parameters
    %screen dimension
        screen = 4/3;                                    %MUST CHANGE
    %fontsize
        fts = 18;
        ftstics = round(fts*3/4);

%Specimen dimensions in in                      %MUST CHANGE
    Rm = 0.8938;                                    %MUST CHANGE
    th = 0.0384;                                    %MUST CHANGE
    Lg = 0.4;                                       %MUST CHANGE
    S = 2 * pi * Rm * th;                           %MUST CHANGE
%Change mm->in
    coef = 1/25.4;

%create the new folder where all the exported files are
folder = sprintf('TT2-%d_FS_%d_SS_%d_OLDCODE',TT2,Facet_size,Step_size);
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
%A(:,6:9) displacement U V W ur=w [mm] 
%A(:,10:13) Surface gradient transformation F(0,0) F(0,1) F(1,0) F(1,1) []
A = load(name);

%Coordinate of the center of the facets in the deformed configuration 
x=(A(:,3)+A(:,6))*coef/Rm;
y=(A(:,4)+A(:,7))*coef/Lg;
z=(A(:,5)+A(:,8))*coef/Rm;

%Calculation for each facet of the Logarithmic cumulative plastic strain
for i = [1:1:size(x,1)]
    %transformation gradient F=RU
F = [[A(i,10),A(i,11)] ;[A(i,12),A(i,13)]];
    %stretching tensor???
U = transpose(F)*F;
    %principal stretch
diagU = eig(U);
    %logarithmic strain in the principal coordinate system
LE_calc = 0.5 * log(diagU);
    %Logarithmic Cumulative Plastic strain
LEp(i) = ( 2./3. * LE_calc(1)^2 + LE_calc(2)^2 + (-LE_calc(1)-LE_calc(2))^2 )^0.5 ;
end

%{
%%Find the max LEp along each "vertical section", that is, along
%%each group of rows with identical i indices (column 2 of A)
maxi=max(A(:,1));
for j=0:maxi;
    k=1;
    for i=1:length(A(:,1)); 
        if A(i,1)==j;               %%Note - not all indices are necessarily represented.  In otherwise, there might be some "gaps" in the counting.
            B{j+1}(k)=LEp(i);
            k=k+1;
        else;
        end;
    end;
end;    

%Remove the empties
emptyCells = cellfun('isempty', B);
B(emptyCells) = [];

for i=1:length(B);
    [emax(i),loc(i)]=max(B{i});
end;

HistData=sprintf('TT2-%d_HistogramData_%d-%d',TT2,Facet_size,Step_size);
save(HistData,'A','LEp','B','emax','loc');
figure;hist(emax,25)
mytitlestring=sprintf('Histogram of Max e^p in Each Column of Comp.Pts-TT2-%d-\\alpha=%.1f-FS%d SS%d-%s',TT2,alpha,Facet_size,Step_size,calctype);
xlabel('Cum. Plastic Strain e^p ','Fontsize',14)
ylabel('Number','Fontsize',14)
title(mytitlestring,'FontSize',10)
set(gca,'XLim',[.8*min(emax) 1.1*max(emax)]);
histtitle=sprintf('TT2-%d_Histogram_%d-%d',TT2,Facet_size,Step_size);
print(gcf,'-dpng', histtitle)
%close
%close
%}
%{ Martin's histogram work

%Alternative method where I first identify rows
%{
%Find horizontal sections by identifying common j indice
maxj=max(A(:,2));
for j=0:maxj;
    k=0;
    for i=1:length(A(:,1)); 
        if A(i,2)==j;
            k=k+1;
            B{j+1}(k)=LEp(i);
        else;
        end;
    end;
end;    

%identify length of each row, and find the min length
for i=1:length(B)
    B_elem_length(i)=length(B{i});
end;
B_min_length=min(B_elem_length);

%truncate all rows to the minimum
for i=1:length(B);
    C{i}=B{i}(1:B_min_length);
end;

%Turn rows into columns
for i=1:length(C);
    for j=1:B_min_length;
        D{j}(i)=C{i}(j);
    end;
end;

%Find the max in each column
for i=1:length(D);
    [emax(i),loc(i)]=max(D{i});
end;        
        
figure
hist(emax)
  %}

%Various histogram and pdf fit functions
%{
To get normal dist mu and sigma:
[mu,sigma]=normfit(emax);  %returns mu and sigma for normal distribution
normaldist=normpdf(emax,mu,sigma); % generates pdf value for each emax value
figure;plot(emax,normaldist)  %plots the normal fit


figure;normplot(emax)      %plots the normal fit (but mu, sig not known)
figure;plot(emax,normaldist,'o')

%Histfit fits the histogram w/ a
distribution, but the characteristic parameters of the distribution are
unknown
figure;histfit(emax,20,'normal')  
figure;histfit(emax,20,'nakagami')
figure;histfit(emax,20,'lognormal')
figure;histfit(emax,20,'gev')    <--------Nice
%}


%%{    
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

%Here is the plot generation and code for the "Do you agree with this"?
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
axis([min(x) max(x) min(y) max(y) min(LEp) max(LEp) min(LEp) max(LEp)])
colorbar('location','EastOutside','Fontsize',fts) 
set (gca,'Fontsize',ftstics)
set(gcf,'Units','Normalized','Outerposition',[0 0 1 1]); 
hold off

agree = input('do you agree with the selected points 0->no / 1->yes');

if agree == 0
end;
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

%Now the final ep_surf plot is actually made
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
%mytitlestring=sprintf('TT2-%d - \\alpha=%.1f - FS%d SS%d - %s (Deviation Check OFF)',TT2,alpha,Facet_size,Step_size,calctype)
mytitlestring=sprintf('TT2-%d - \\alpha=%.1f - FS%d SS%d - %s',TT2,alpha,Facet_size,Step_size,calctype);
title(mytitlestring,'FontSize',14)
axis([min(x) max(x) min(y) max(y) min(LEp) max(LEp) min(LEp) max(LEp)])
colorbar('location','EastOutside','Fontsize',fts) 
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
    clear NEx;
    clear NEy;
    clear NExy;
    clear gamma;
    clear LEp;
    clear LE_calc;
    %name of the file
    name = sprintf('%s%s_%d.dat',PATH,prefix,k);
    %open the file A = matrix (nx13) with
    %A(:,1) index i [] 
    %A(:,2) index j [] 
    %A(:,3:5) undeformed coordinates X Y Z [mm] 
    %A(:,6:9) displacement U V W ur=w [mm] 
    %A(:,10:13) Surface gradient transformation F(0,0) F(0,1) F(1,0) F(1,1) []
    A = load(name);
    %Find the facets
    top10 = [];
    for i = 1 : size(index,1)
    top10 = [top10;find( A(:,1)==index(i,1)  & A(:,2)==index(i,2) )];
    end
    for i = 1 : size(top10)
    %transformation gradient F=RU
    F = [[A(top10(i),10),A(top10(i),11)] ;[A(top10(i),12),A(top10(i),13)]];
    %stretching tensor???
    U = transpose(F)*F;
    %principal stretch
    diagU = eig(U);
    %logarithmic strain in the principal coordinate system
    LE_calc = 0.5 * log(diagU);
    %Logarithmic Cumulative Plastic strain
    LEp(i) = ( 2./3. * LE_calc(1)^2 + LE_calc(2)^2 + (-LE_calc(1)-LE_calc(2))^2 )^0.5 ;
    %Rotation tensor
    R = F * U^(-0.5);
    %Technical strain calculation in the stretching coordinate system
    NE_calc = U^0.5 - eye([2 2]);
    %Technical strain in the x',y' coordinate system with z'=normal to the
    %surface, x' parallel to the (X,Z) plane (Y axis of the cylinder in the undeformed stage), (x',y',z') orthonormal
    %coordinate system
    NE_calc_rot = R * NE_calc * transpose(R);
    NEx(i) = NE_calc_rot(1,1);
    NEy(i) = NE_calc_rot(2,2);
    NExy(i) = NE_calc_rot(1,2);
    gamma(i) = atan(NExy(i)/(1+NEx(i))) + atan(NExy(i)/(1+NEy(i)));
    end
    export_max = [export_max;[Force(k+1,2), Force(k+1,3)/S/1000,NEx(1),NEy(1),abs(NExy(1)),abs(gamma(1)),LEp(1)]];
    export_mean = [export_mean;[Force(k+1,2), Force(k+1,3)/S/1000,nanmean(NEx),nanmean(NEy),abs(nanmean(NExy)),abs(nanmean(gamma)),nanmean(LEp)]];
    export_stdv = [export_stdv;[Force(k+1,2), Force(k+1,3)/S/1000,std(NEx),std(NEy),std(NExy),std(gamma),std(LEp)]];
    %Profile of the cumulative plastic strain along y for data_force_points
    %set
    if any(data_Force_point==k) == true
      index_ep = find(A(:,1)==index(1,1));
        for l = 1 : size(index_ep,1)
      %calculate Y coordinate
      pro_LEp(l,1,profile) = (A(index_ep(l),4)+A(index_ep(l),7))*coef/th;     
      %transformation gradient F=RU
      F = [[A(index_ep(l),10),A(index_ep(l),11)] ;[A(index_ep(l),12),A(index_ep(l),13)]];
      %stretching tensor???
      U = transpose(F)*F;
      %principal stretch
      diagU = eig(U);
      %logarithmic strain in the principal coordinate system
      LE_calc = 0.5 * log(diagU);
      %Logarithmic Cumulative Plastic strain
      pro_LEp(l,2,profile) = ( 2./3. * LE_calc(1)^2 + LE_calc(2)^2 + (-LE_calc(1)-LE_calc(2))^2 )^0.5 ;      
        end
      profile=profile+1;
    end  
    %Profile of the cumulative plastic strain along y for last_pictures set
      if any(last_pictures==k) == true
      index_ep2 = find(A(:,1)==index(1,1));
        for l = 1 : size(index_ep2,1)
      %calculate Y coordinate
      pro_LEp2(l,1,profile2) = (A(index_ep2(l),4)+A(index_ep2(l),7))*coef/th;     
      %transformation gradient F=RU
      F = [[A(index_ep2(l),10),A(index_ep2(l),11)] ;[A(index_ep2(l),12),A(index_ep2(l),13)]];
      %stretching tensor???
      U = transpose(F)*F;
      %principal stretch
      diagU = eig(U);
      %logarithmic strain in the principal coordinate system
      LE_calc = 0.5 * log(diagU);
      %Logarithmic Cumulative Plastic strain
      pro_LEp2(l,2,profile2) = ( 2./3. * LE_calc(1)^2 + LE_calc(2)^2 + (-LE_calc(1)-LE_calc(2))^2 )^0.5 ;      
        end
      profile2=profile2+1;
    end   
end        



%Plot time-cumulative plastic strain
figure
errorbar(export_mean(:,1),export_mean(:,7),1.5*export_stdv(:,7));
hold on
plot(export_max(:,1),export_max(:,7),'ro')
%mytitlestring=sprintf('TT2-%d - \\alpha=%.1f - FS%d SS%d - %s (Deviation Check OFF)',TT2,alpha,Facet_size,Step_size,calctype)
mytitlestring=sprintf('TT2-%d - \\alpha=%.1f - FS%d SS%d - %s',TT2,alpha,Facet_size,Step_size,calctype);
title(mytitlestring,'FontSize',12)
legend('mean value over 10','max value','Location','NorthWest')
axis([0 max(export_max(:,1)) 0 1.1*max(export_max(:,7))])
xlabel('time [s]','Fontsize',fts)
ylabel('cumulative plastic strain e^p []','Fontsize',fts)
set (gca,'Xgrid','on','Ygrid','on','LineWidth',2)
set (gca,'Fontsize',ftstics)
set(gcf,'Units','Normalized','Outerposition',[0 0 1/screen 1]);
hold off
print(gcf,'-dpng','time_ep.png');

%plot gamma,epsilon y, 
figure
plot(export_mean(:,6),export_mean(:,4),'bx');
hold on
plot(export_max(:,6),export_max(:,4),'ro')
%mytitlestring=sprintf('TT2-%d - \\alpha=%.1f - FS%d SS%d - %s (Deviation Check OFF)',TT2,alpha,Facet_size,Step_size,calctype)
mytitlestring=sprintf('TT2-%d - \\alpha=%.1f - FS%d SS%d - %s',TT2,alpha,Facet_size,Step_size,calctype);
title(mytitlestring,'FontSize',12)
legend('mean value over 10','max value','Location','SouthWest','Fontsize',2)
xlabel('\gamma [rad]','Fontsize',fts)
ylabel('\epsilon _y []','Fontsize',fts)
set (gca,'Xgrid','on','Ygrid','on','LineWidth',2)
set (gca,'Fontsize',ftstics)
set(gcf,'Units','Normalized','Outerposition',[0 0 1/screen 1]);
hold off
print(gcf,'-dpng','gamma_epsy.png');

%plot gamma,epsilon y, highlighting last 10
figure
plot(export_mean(:,6),export_mean(:,4),'bx');
hold on
plot(export_max(:,6),export_max(:,4),'ro')
%fillcolors=['m','c','r','g','b','k','m','c','r','k'];
fillcolors={[205/255 0 0] [255 165 0]/255 [255 215 0]/255 [205 205 0]/300 [0 .5 0] [0 1 1]*.8 [0 0 1] [153 50 204]/255 [110 123 139]/255 [0 0 0]};
for i=0:9;
    plot(export_mean((last-i),6),export_mean((last-i),4),'Marker','p','MarkerEdgeColor',fillcolors{i+1},'MarkerFaceColor',fillcolors{i+1})
end;
for i=0:9;
    plot(export_max((last-i),6),export_max((last-i),4),'Marker','o','MarkerEdgeColor',fillcolors{i+1},'MarkerFaceColor',fillcolors{i+1})
end;
%mytitlestring=sprintf('TT2-%d - \\alpha=%.1f - FS%d SS%d - %s (Deviation Check OFF)',TT2,alpha,Facet_size,Step_size,calctype)
mytitlestring=sprintf('TT2-%d - \\alpha=%.1f - FS%d SS%d - %s',TT2,alpha,Facet_size,Step_size,calctype);
title(mytitlestring,'FontSize',12)
legend('mean value over 10','max value','Final Stage mean','2nd-to-final mean','3rd to final mean','4th to final mean','5th to final mean','6th to final mean','7th to final mean','8th to final mean','9th to final mean','10th to final mean',...
    'Final stage max','2nd to final max','3rd to final max','4th to final max','5th to final max','6th to final max','7th to final max','8th to final max','9th to final max','10th to final max','Location','BestOutside','FontWeight','light')
xlabel('\gamma [rad]','Fontsize',fts)
ylabel('\epsilon _y []','Fontsize',fts)
set (gca,'Xgrid','on','Ygrid','on','LineWidth',2)
set (gca,'Fontsize',ftstics)
xmin=min(min(export_max(length(export_max)-9:length(export_max),6),export_mean(length(export_mean)-9:length(export_mean),6)))*.5;
ymin=min(min(export_max(length(export_max)-9:length(export_max),4),export_mean(length(export_mean)-9:length(export_mean),4)))*.5;
xmax=max(max(export_max(length(export_max)-9:length(export_max),6),export_mean(length(export_mean)-9:length(export_mean),6)))*1.1;
ymax=max(max(export_max(length(export_max)-9:length(export_max),4),export_mean(length(export_mean)-9:length(export_mean),4)))*1.1;
axis([xmin xmax ymin ymax])
set(gcf,'Units','Normalized','Outerposition',[0 0 1/screen 1]);
hold off
print(gcf,'-dpng','gamma_epsy_last10.png');

%plot y,e^p for data_force_points set
figure
set(gcf,'Units','Normalized','Outerposition',[0 0 1 1]);
    %above makes the figure full screen
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
%mytitlestring=sprintf('TT2-%d - \\alpha=%.1f - FS%d SS%d - %s (Deviation Check OFF)',TT2,alpha,Facet_size,Step_size,calctype)
mytitlestring=sprintf('TT2-%d - \\alpha=%.1f - FS%d SS%d - %s',TT2,alpha,Facet_size,Step_size,calctype);
title(mytitlestring,'FontSize',12)
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
marker_style = cellstr(marker_style);

plot(pro_LEp2(:,2,1),pro_LEp2(:,1,1),char(marker_style(1)));
%mytitlestring=sprintf('TT2-%d - \\alpha=%.1f - FS%d SS%d - %s (Deviation Check OFF)',TT2,alpha,Facet_size,Step_size,calctype)
mytitlestring=sprintf('TT2-%d - \\alpha=%.1f - FS%d SS%d - %s',TT2,alpha,Facet_size,Step_size,calctype);
title(mytitlestring,'FontSize',12)
hold on
for i = 2 : size(last_pictures,1)
plot(pro_LEp2(:,2,i),pro_LEp2(:,1,i),char(marker_style(i)));
end
xlabel('cumulative plastic strain e^p []','Fontsize',fts)
ylabel('axial coordinate Y/th []','Fontsize',fts)
%mytitlestring=sprintf('Max value of e^p in last five pictures - TT2-%d - \\alpha=%.1f - FS%d SS%d - %s (Deviation Check OFF)',TT2,alpha,Facet_size,Step_size,calctype)
mytitlestring=sprintf('Max value of e^p in last five pictures - TT2-%d - \\alpha=%.1f - FS%d SS%d - %s',TT2,alpha,Facet_size,Step_size,calctype);
title(mytitlestring,'FontSize',10)
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

cd('C:\Users\mfs279\Documents\MATLAB')
%}