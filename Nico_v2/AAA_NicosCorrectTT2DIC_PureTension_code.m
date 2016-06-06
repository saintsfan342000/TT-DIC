% USE THIS FOR PURE TENSION

%Nicolas Tardif
%01/07/2013 (dd,mm,yyyy)
%Post processing of the 3D DIC results when exported from a linear strain calculation using Nico_3D_linear3_fv 
%
clear
close all;

%path for the files that the PERL program produces
    PATH = 'E:\Martin_Experiments\TT2-20_Oct3_a0p5\SBS Computation\perl_post'; %MUST CHANGE
%prefix of the name
    prefix = TT2_20_DC15_FS19_SS6'; %MUST CHANGE don't include underscore after linear
%DIC paramaters
    Facet_size = 19;%pix                            %MUST CHANGE                               %%VERY IMPORTANT MUST CHANGE
    Step_size = 6;%pix                              %MUST CHANGE
% # of the last picture
    last = 580; %MUST CHANGE
    TT2=20;
    alpha=0.5;
    calctype='Linear';
%plotting parameters
    %screen dimension
        screen = 4/3;                                    %MUST CHANGE
    %fontsize
        fts = 18;
        ftstics = round(fts*3/4);
%Specimen dimensions in in                      %MUST CHANGE
    Rm = 0.8940;                                    %MUST CHANGE
    th = 0.0379;                                    %MUST CHANGE
    Lg = 0.4;                                       %MUST CHANGE
        Measured_output = 'Axial';
        Stress_normalized = 2 * pi  * Rm * th;
        cal = 5; %1V = 5 kips
        displayed = '\sigma_Y [ksi]' 
%Change mm->in
    coef = 1/25.4;
%Size of the side of the Area of calculation of Strain : Scott (grid method)
    Size_Scott = 1/16; %[in]
%create the new folder where all the exported files are
folder = sprintf('export_FS_%d_SS_%d',Facet_size,Step_size);
mkdir(folder);
folder2 = sprintf('export_FS_%d_SS_%d\profile',Facet_size,Step_size);
mkdir(folder2);

%Watch result
watch = [1:50:last,last];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%1-Determination of the area to scan.%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%name of the file
name = sprintf('%s%s_%d.dat',PATH,prefix,last);
%open the file A = matrix (nx13) with
%A(:,1) index i [] corresponding to the x axis
%A(:,2) index j [] corresponding to the y axis
%A(:,3:5) undeformed coordinates X Y Z [mm] 
%A(:,6:9) displacement U V W ur=w [mm] 
%A(:,10:13) Surface gradient transformation F(0,0) F(0,1) F(1,0) F(1,1) []
A = load(name);

%Coordinate of the center of the facets in the undeformed configuration 
x=(A(:,3) * coef); % [in]
y=(A(:,4)) * coef; % [in]
z=(A(:,5)) *coef; % [in]

%Calculation for each facet of the Logarithmic cumulative plastic strain
for i = [1:1:size(x,1)]                    
    F = [[A(i,10),A(i,11)] ;[A(i,12),A(i,13)]];
    %stretching tensor???
    U = transpose(F)*F;
    %principal stretch
    [RotU,diagU] = eig(U); 
    %logarithmic strain in the principal coordinate system
    LE_calc = 0.5 * [[log(diagU(1,1)),0];[0,log(diagU(2,2))]];
    %Logarithmic Cumulative Plastic strain
    LEp(i) = ( 2./3. * ( LE_calc(1,1)^2 + LE_calc(2,2)^2 + (-LE_calc(1,1)-LE_calc(2,2))^2 ))^0.5 ;
end

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

%projection of the cumulative plastic strain on the grid
LEP=griddata(x,y,LEp,X,Y);

close all;
%plot the data
contourf(X,Y,LEP,20)
 hold on
shading flat
xlabel('undeformed x [in]','Fontsize',fts)
ylabel('undeformed y [in] (axial)','Fontsize',fts)
axis([min(x) max(x) min(y) max(y) min(LEp) max(LEp) min(LEp) max(LEp)])
colorbar('location','EastOutside','Fontsize',fts) 
set (gca,'Fontsize',ftstics)
set(gcf,'Units','Normalized','Outerposition',[0 0 1 1]); 
hold off

disp('crop the zone you want to study (ie highest strain zone)');
crop = ginput(2);
Xmin_crop = min(crop(:,1)) ; 
Xmax_crop = max(crop(:,1)) ;
Ymin_crop = min(crop(:,2)) ;
Ymax_crop = max(crop(:,2)) ;

box_x = [Xmin_crop , Xmax_crop , Xmax_crop , Xmin_crop , Xmin_crop]';
box_y = [Ymin_crop , Ymin_crop , Ymax_crop , Ymax_crop , Ymin_crop]';
box_z = [max(LEp) , max(LEp) , max(LEp) , max(LEp) , max(LEp)]';

close all;
%plot the data
contourf(X,Y,LEP,20)
 hold on
plot3(box_x,box_y,box_z,'y-o','MarkerEdgeColor','y','MarkerFaceColor','y') 
xlabel('undeformed x [in]','Fontsize',fts)
ylabel('undeformed y [in] (axial)','Fontsize',fts)
zlabel('e^p','Fontsize',fts)
mytitlestring=sprintf('TT2-%d - \\alpha=%.2f - FS%d SS%d - %s',TT2,alpha,Facet_size,Step_size,calctype);
shading flat
title(mytitlestring,'FontSize',14)
axis([min(x) max(x) min(y) max(y) min(LEp) max(LEp) min(LEp) max(LEp)])
colorbar('location','EastOutside','Fontsize',fts) 
set (gca,'Fontsize',ftstics)
set(gcf,'Units','Normalized','Outerposition',[0 0 1 1]); 
hold off
save_path = sprintf('%s/last_picture',folder);
print(gcf,'-dpng',save_path); %save the data



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



%Definition of the size of the averaging zone for grid method
size_av = Size_Scott - th * Facet_size_th;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%2-Calculate strain values in the scanned area.%%%%%%%%%%%%%%%%
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

%initialization
export_max=[];
export_mean=[];
export_stdv=[];
export_Scott=[];

for k = 1 : last
    compt=0; %count the number of scanned columns
    compt2=0; %count the number of point in the Scott method
    clear A;
    clear LEp_kept;
    clear NEx_kept;
    clear NEy_kept;
    clear NExy_kept;
    clear gamma_kept;
    %name of the file
    name = sprintf('%s%s_%d.dat',PATH,prefix,k);
    %open the file A = matrix (nx13) with
    %A(:,1) index i [] corresponding to the x axis
    %A(:,2) index j [] corresponding to the y axis (axial)
    %A(:,3:5) undeformed coordinates X Y Z [mm] 
    %A(:,6:9) displacement U V W ur=w [mm] 
    %A(:,10:13) Surface gradient transformation F(0,0) F(0,1) F(1,0) F(1,1) []
    A = load(name);
    %Find the facets corresponding to the scanned zone
    scan = find( A(:,3)*coef >= Xmin_crop & A(:,3)*coef <= Xmax_crop & A(:,4)*coef >= Ymin_crop & A(:,4)*coef <= Ymax_crop);  
    %Find the index i that are used to calculate the strain ie without
    %index j missing
        for l = unique(A(scan,1))'
            column = find( A(scan,1) == l );
            index_j = A(scan(column),2);
            if (size(index_j,1) == max(index_j) - min(index_j) +1)
                clear LEp;
                clear NEx;
                clear NEy;
                clear NExy;
                clear gamma;
                compt = compt+1;
                for m = 1 : size(index_j,1)
                   %transformation gradient F=RU
                    F = [[A(scan(column(m)),10),A(scan(column(m)),11)] ;[A(scan(column(m)),12),A(scan(column(m)),13)]]; 
                    %stretching tensor???
                    U = transpose(F)*F;
                    %principal stretch
                    [RotU,diagU] = eig(U); 
                    %logarithmic strain in the principal coordinate system
                    LE_calc = 0.5 * [[log(diagU(1,1)),0];[0,log(diagU(2,2))]];
                    %Logarithmic Cumulative Plastic strain
                    LEp(m) = ( 2./3. * ( LE_calc(1,1)^2 + LE_calc(2,2)^2 + (-LE_calc(1,1)-LE_calc(2,2))^2 ))^0.5 ;
                    %Rotation tensor
                    %R = F * U^(-0.5);
                    %Technical strain calculation in the stretching coordinate system
                    NE_calc = (U^0.5 - eye([2 2]));
                    %Technical strain in the x',y' material coordinate
                    %system
                    x_undef(m) = A(scan(column(m)),3)*coef;
                    y_undef(m) = A(scan(column(m)),4)*coef;
                    NEx(m) = NE_calc(1,1);
                    NEy(m) = NE_calc(2,2);
                    NExy(m) = NE_calc(1,2);
                    gamma(m) = atan(NExy(m)/(1+NEx(m))) + atan(NExy(m)/(1+NEy(m)));
                end 
            %index for which we find the maximum equivalent strain
            [LEp_kept(compt) , index_m] = max(LEp);
            NEx_kept(compt) = NEx(index_m);
            NEy_kept(compt) = NEy(index_m);
            gamma_kept(compt) = gamma(index_m);
            x_undef_kept(compt) = x_undef(index_m);
            y_undef_kept(compt) = y_undef(index_m);
            end
        end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%        
%Calculation of the max equivalent strain and so on using a filtering method based on the strain ratio           
    ratio = NEx_kept ./ NEy_kept;
    mean_ratio = nanmean(ratio);
    stdv_ratio = nanstd(ratio);
    %Criteria to remove points that have a bad strain ratio
    compt_kept2 = find(ratio >= mean_ratio - 0.5 * stdv_ratio & ratio <= mean_ratio + 0.5 * stdv_ratio);
    
    [max_LEP, index_max] = max(LEp_kept(compt_kept2));
     export_max = [export_max;[size(compt_kept2,2), Force(k+1,2), cal * Force(k+1,3)/Stress_normalized,NEx_kept(compt_kept2(index_max)),NEy_kept(compt_kept2(index_max)),abs(gamma_kept(compt_kept2(index_max))),LEp_kept(compt_kept2(index_max))]];
     export_mean = [export_mean;[size(compt_kept2,2), Force(k+1,2), cal * Force(k+1,3)/Stress_normalized,nanmean(NEx_kept(compt_kept2)),nanmean(NEy_kept(compt_kept2)),abs(nanmean(gamma_kept(compt_kept2))),nanmean(LEp_kept(compt_kept2))]];
     export_stdv = [export_stdv;[size(compt_kept2,2), Force(k+1,2), cal * Force(k+1,3)/Stress_normalized,nanstd(NEx_kept(compt_kept2)),nanstd(NEy_kept(compt_kept2)),abs(nanstd(gamma_kept(compt_kept2))),nanstd(LEp_kept(compt_kept2))]];
%Profile of ep for different data_force_points :
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%     
    if any(data_Force_point==k) == true
      index_MAX = find( A(:,3) * coef == x_undef_kept(compt_kept2(index_max)) & A(:,4) * coef == y_undef_kept(compt_kept2(index_max)));
      index_ep = find( A(:,1) == A(index_MAX,1))
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
      pro_LEp(l,2,profile) = ( 2./3. * (LE_calc(1)^2 + LE_calc(2)^2 + (-LE_calc(1)-LE_calc(2))^2 ))^0.5 ;      
        end
      profile=profile+1;
    end  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%     
%Calculation of the max equivalent strain and so on using an average value
%over an area equivalent of the one of Scott (grid method)
    for i = size(x_undef_kept)
       clear av_index;
       clear F;
       clear U;
       clear RotU;
       clear diagU;
       clear LE_calc;
       clear LEp;
       clear NE_calc;
       clear x_undef;
       clear y_undef;
       clear NEx;
       clear NEy;
       clear NExy;
       clear gamma;
       compt2 = compt2 + 1;
       av_index = find( A(:,3)*coef <= x_undef_kept(i) + size_av/2 &  A(:,3)*coef >= x_undef_kept(i) - size_av/2  & A(:,4)*coef <= y_undef_kept(i) + size_av/2 &  A(:,4)*coef >= y_undef_kept(i) - size_av/2); 
                for m = 1 : size(av_index,1)
                   %transformation gradient F=RU
                    F = [[A(av_index(m),10),A(av_index(m),11)] ;[A(av_index(m),12),A(av_index(m),13)]]; 
                    %stretching tensor???
                    U = transpose(F)*F;
                    %principal stretch
                    [RotU,diagU] = eig(U); 
                    %logarithmic strain in the principal coordinate system
                    LE_calc = 0.5 * [[log(diagU(1,1)),0];[0,log(diagU(2,2))]];
                    %Logarithmic Cumulative Plastic strain
                    LEp(m) = ( 2./3. * ( LE_calc(1,1)^2 + LE_calc(2,2)^2 + (-LE_calc(1,1)-LE_calc(2,2))^2 ))^0.5 ;
                    %Rotation tensor
                    %R = F * U^(-0.5);
                    %Technical strain calculation in the stretching coordinate system
                    NE_calc = (U^0.5 - eye([2 2]));
                    %Technical strain in the x',y' material coordinate
                    %system
                    x_undef(m) = A(av_index(m),3)*coef;
                    y_undef(m) = A(av_index(m),4)*coef;
                    NEx(m) = NE_calc(1,1);
                    NEy(m) = NE_calc(2,2);
                    NExy(m) = NE_calc(1,2);
                    gamma(m) = atan(NExy(m)/(1+NEx(m))) + atan(NExy(m)/(1+NEy(m)));
                end
            size2_Scott(compt2) = Facet_size_th * th + mean( [max(x_undef) - min(x_undef) , max(y_undef) - min(y_undef)] );
            NEx_Scott(compt2) = nanmean(NEx); 
            NEy_Scott(compt2) = nanmean(NEy);
            gamma_Scott(compt2) = nanmean(gamma);
            LEp_Scott(compt2) = nanmean(LEp);
    end   
    export_Scott = [export_Scott;[nanmean(size2_Scott), Force(k+1,2), cal * Force(k+1,3) / Stress_normalized,nanmean(NEx_Scott),nanmean(NEy_Scott),abs(nanmean(gamma_Scott)),nanmean(LEp_Scott)]];
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%     
     %Part to be removed if speed is needed
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        if isempty(find(watch==k)) == 0
        clear x_graph;
        clear y_graph;
        clear F_graph;
        clear U_graph;
        clear RotU_graph;
        clear diagU_graph;
        clear LE_calc_graph;
        clear LEp_graph;
        clear x_edge;
        clear y_edge;
        clear X_graph;
        clear Y_graph;
        clear LEP_graph;
        %graph
        %Plot time-cumulative plastic strain
        %Coordinate of the center of the facets in the undeformed configuration 
        x_graph=(A(:,3) * coef); % [in]
        y_graph=(A(:,4)) * coef; % [in]
        %Calculation for each facet of the Logarithmic cumulative plastic strain
        for i = [1:1:size(x_graph,1)]                    
        F_graph = [[A(i,10),A(i,11)] ;[A(i,12),A(i,13)]];
        %stretching tensor???
        U_graph = transpose(F_graph)*F_graph;
        %principal stretch
        [RotU_graph,diagU_graph] = eig(U_graph); 
        %logarithmic strain in the principal coordinate system
        LE_calc_graph = 0.5 * [[log(diagU_graph(1,1)),0];[0,log(diagU_graph(2,2))]];
        %Logarithmic Cumulative Plastic strain
        LEp_graph(i) = ( 2./3. * ( LE_calc_graph(1,1)^2 + LE_calc_graph(2,2)^2 + (-LE_calc_graph(1,1)-LE_calc_graph(2,2))^2 ))^0.5 ;
        end
        x_edge=[min(x_graph):dx:max(x_graph)];
        y_edge=[min(y_graph):dy:max(y_graph)];
        [X_graph,Y_graph]=meshgrid(x_edge,y_edge);
        LEP_graph=griddata(x_graph,y_graph,LEp_graph',X_graph,Y_graph);
        clf();
        %plot the data
        contourf(X_graph,Y_graph,LEP_graph,20)
        hold on
        plot3(x_undef_kept(compt_kept2),y_undef_kept(compt_kept2),LEp_kept(compt_kept2),'o','MarkerEdgeColor','w','MarkerFaceColor','w') 
        plot3(box_x,box_y,box_z,'w-') ;
        shading flat
        xlabel('undeformed x [in]','Fontsize',fts)
        ylabel('undeformed y [in] (axial)','Fontsize',fts)
        zlabel('e^p','Fontsize',fts)
        mytitlestring=sprintf('stage %d over %d',k,last);
        title(mytitlestring,'FontSize',14)
        axis([min(x) max(x) min(y) max(y) min(LEp) max(LEp) min(LEp) max(LEp)])
        colorbar('location','EastOutside','Fontsize',fts) 
        set (gca,'Fontsize',ftstics)
        set(gcf,'Units','Normalized','Outerposition',[0 0 1 1]); 
        hold off
        pause(2);
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end    
   

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%3-PLOT DATA%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%plot gamma,epsilon y, 
figure
plot(export_mean(:,4),export_mean(:,5),'bx');
hold on
plot(export_max(:,4),export_max(:,5),'ro');
plot(export_Scott(:,4),export_Scott(:,5),'k+');
%mytitlestring=sprintf('TT2-%d - \\alpha=%.1f - FS%d SS%d - %s (Deviation Check OFF)',TT2,alpha,Facet_size,Step_size,calctype)
mytitlestring=sprintf('Strain in the rotating material coordinate system \n TT2-%d - \\alpha=%s - FS%d SS%d - %s \n Averaging size for like-Scott results = %f in',TT2,alpha,Facet_size,Step_size,calctype,size2_Scott(compt2));
title(mytitlestring,'FontSize',12)
legend('mean value','max value', 'averaged value (like-Scott)' , 'Location','SouthWest','Fontsize',2)
xlabel('\epsilon_x [rad]','Fontsize',fts)
ylabel('\epsilon_y []','Fontsize',fts)
set (gca,'Xgrid','on','Ygrid','on','LineWidth',2)
set (gca,'Fontsize',ftstics)
set(gcf,'Units','Normalized','Outerposition',[0 0 1/screen 1]);
hold off 
save_path = sprintf('%s/gamma_epsy',folder);
print(gcf,'-dpng',save_path);



     
 %Plot time-cumulative plastic strain
figure
errorbar(export_mean(:,2),export_mean(:,7),0.5*export_stdv(:,7));
hold on
plot(export_max(:,2),export_max(:,7),'ro')
plot(export_Scott(:,2),export_Scott(:,7),'k+');
%mytitlestring=sprintf('TT2-%d - \\alpha=%.1f - FS%d SS%d - %s (Deviation Check OFF)',TT2,alpha,Facet_size,Step_size,calctype)
mytitlestring=sprintf('TT2-%d - \\alpha=%s - FS%d SS%d - %s \n Averaging size for like-Scott results = %f in',TT2,alpha,Facet_size,Step_size,calctype,size2_Scott(compt2))
title(mytitlestring,'FontSize',12)
legend('mean value','max value', 'averaged value (like-Scott)' ,'Location','NorthWest')
axis([0 max(export_max(:,2)) 0 1.1*max(export_max(:,7))])
xlabel('time [s]','Fontsize',fts)
ylabel('cumulative plastic strain e^p []','Fontsize',fts)
set (gca,'Xgrid','on','Ygrid','on','LineWidth',2)
set (gca,'Fontsize',ftstics)
set(gcf,'Units','Normalized','Outerposition',[0 0 1/screen 1]);
hold off  
save_path = sprintf('%s/time_ep.png',folder);
print(gcf,'-dpng',save_path);  

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
plot(Force(:,2),Force(:,3)*cal / Stress_normalized,'LineWidth',2);
%mytitlestring=sprintf('TT2-%d - \\alpha=%.1f - FS%d SS%d - %s (Deviation Check OFF)',TT2,alpha,Facet_size,Step_size,calctype)
mytitlestring=sprintf('TT2-%d - \\alpha=%s - FS%d SS%d - %s',TT2,alpha,Facet_size,Step_size,calctype);
title(mytitlestring,'FontSize',12)
hold on
for i=1:nb_profile
plot(Force(data_Force_point(i),2),Force(data_Force_point(i),3) * cal / Stress_normalized,char(marker_style(i)),'MarkerSize',10);
end
axis([0 max(Force(:,2)) 0 1.1*max(Force(:,3)) * cal / Stress_normalized])
xlabel('time [s]','Fontsize',fts)
ylabel(displayed,'Fontsize',fts)
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
save_path = sprintf('%s/ep_profile',folder);
print(gcf,'-dpng',save_path)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%4-EXPORT DATA%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


save_path = sprintf('%s/max.dat',folder);
output = fopen(save_path,'w');
fprintf(output,'%s nb of points [] Time [s] %s Nominal Stress [ksi] Technical NEx [] NEy [] Gamma [rad] True cumulative strain [] \n','%',Measured_output);
fprintf(output, '%i %f %f %f %f %f %f\n', export_max');
fclose(output);

save_path = sprintf('%s/mean.dat',folder);
output = fopen(save_path,'w');
fprintf(output,'%s nb of points [] Time [s] %s Nominal Stress [ksi] Technical NEx [] NEy [] Gamma [rad] True cumulative strain [] \n','%',Measured_output);
fprintf(output, '%i %f %f %f %f %f %f\n', export_mean');
fclose(output);

save_path = sprintf('%s/stdv.dat',folder);
output = fopen(save_path,'w');
fprintf(output,'%s nb of points [] Time [s] %s Nominal Stress [ksi] Technical NEx [] NEy [] Gamma [rad] True cumulative strain [] \n','%',Measured_output);
fprintf(output, '%i %f %f %f %f %f %f\n', export_stdv');
fclose(output);

save_path = sprintf('%s/like_Scott.dat',folder);
output = fopen(save_path,'w');
fprintf(output,'%s size of the averaging zone for like-Scott [in] Time [s] %s Nominal Stress [ksi] Technical NEx [] NEy [] Gamma [rad] True cumulative strain [] \n','%',Measured_output);
fprintf(output, '%i %f %f %f %f %f %f\n', export_Scott');
fclose(output);



for i=1:nb_profile
  clear out;
  name = sprintf('%s/export_Y_LEp_time_%6.3f_stress_%6.3f.dat',folder2,Force(data_Force_point(i),2),Force(data_Force_point(i),3)*cal/Stress_normalized);
  output = fopen(name,'w');
  fprintf(output,'%s Y/th [] True Cumulative plastic strain [] \n','%');
  out = pro_LEp(:,:,i);
  fprintf(output,'%f %f',out');
  save(name,'out','-ASCII');
end

  
%Writing of the facet size and step size compared to the thicknes of the
%specimen
save_path = sprintf('%s/facet_step.txt',folder);
fileID = fopen(save_path,'w');
fprintf(fileID,'Facet size/th = %1.4f \n Step size/th = %1.4f',Facet_size_th,Step_size_th);
fclose(fileID);   

 