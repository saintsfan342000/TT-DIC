function TT2_DIC_Analysis_contours;
clear;
close all;

% For Proposal : Stg 200, 

%%%%%% PRELIMINARY DATA to ENTER %%%%%% 
    last = 609;         %Last stage in which specimen is not fractures
    TT2=17;             %Expt number
    alpha=1.0;          %Alpha
    Rm = 0.8942;        % Mean radius
    thickness = 0.0382;  %Wall thickness
    curdir=pwd;
    addpath(sprintf('%s\\Matlab\\extras',curdir(1:2)));     %Adds export_fig
    PATH = 'E:\Martin_Experiments\AAA_TensionTorsion\TT2-17\AramisExport_MissingRemoved';
    prefix = 'TT2-17_DC15_FS19_SS6_';                %MUST CHANGE don't include underscore after linear
    savepath = 'E:\Martin_Experiments\AAA_TensionTorsion\Paper and Conference\Base-Case Figures\TT2-17-Contours';
    
    % LIMIT LOAD IS AT STAGE 436
    
    %STAGES = [200   320  444   468   492   516   540   564   588   609];
    STAGES = 330;

%fontsize
fts = 18;                       %Font size for figures
ftstics = round(fts*3/4);       %Font size for figure axes
        
%Change mm->in
MMtoIn = 1/25.4;
    

% Find facet spacing to setup meshgrid increment
    name = sprintf('%s\\%s%d.dat',PATH,prefix,last);                  %Name of the stage file to open
    A=load(name);
    x=A(:,3) .* MMtoIn;     % Output file is in mm.   Converts to inches.
    y=A(:,4) .* MMtoIn;
    z=A(:,5) .* MMtoIn;
    %Creation of a grid in order to plot the data
    %Determination of the min distance between any two facets in the stage
    min_disp = 1000;    %Initialize to erroneously large value
    for k = [2 : round(length(x)/30) : size(x,1)-1] %i will go from 2 to length(X)-1 in 30 equally spaced increments
        for q= [1:1:k-1,k+1:1:size(x,1)]    %j will cover all indices except i itself
            min_disp = min ( min_disp , ( (x(k)-x(q))^2 + (y(k)-y(q) )^2 + ( z(k)-z(q) )^2 )^0.5 ); %Recursive process
        end
    end
    dx=(min_disp / 5);  
    dy=(min_disp / 5);
    
% Plot force and points
stf = load(sprintf('%s\\time_force.dat',PATH));
plot(stf(:,1),stf(:,3));
hold on
plot(stf(STAGES+1,1),stf(STAGES+1,3),'ro','markerfacecolor','r')
%export_fig(sprintf('%s\\Stages.jpg',savepath))    
zz = 1;  
%figure
newcaxis = load('E:\Martin_Experiments\AAA_TensionTorsion\Paper and Conference\Base-Case Figures\TT2-17-Contours\CaxisAndTix.dat');
for k= STAGES;
    %for k = 320
     %if k ~= 609
            clear xpoint ypoint F_graph U_graph RotU_graph diagU_graph LE_calc_graph LEp_graph x_edge y_edge X_graph Y_graph LEP_graph box_z;
            name = sprintf('%s\\%s%d.dat',PATH,prefix,k);                  %Name of the stage file to open
            A = load(name);
            xpoint=(A(:,3) ) .* MMtoIn ; 
            ypoint=(A(:,4)) .* MMtoIn ; 
            for i = 1:length(xpoint)                                        %Calculation for each facet of the Logarithmic cumulative plastic strain
                F_graph = [[A(i,10),A(i,11)] ;[A(i,12),A(i,13)]];
                U_graph = transpose(F_graph)*F_graph;                       %stretching tensor???
                [~,diagU_graph] = eig(U_graph);                    %principal stretch
                LE_calc_graph = 0.5 * [[log(diagU_graph(1,1)),0];[0,log(diagU_graph(2,2))]];    %logarithmic strain in the principal coordinate system
                LEp_graph(i) = ( 2./3. * ( LE_calc_graph(1,1)^2 + LE_calc_graph(2,2)^2 + (-LE_calc_graph(1,1)-LE_calc_graph(2,2))^2 ))^0.5 ;    %Logarithmic Cumulative Plastic strain
            end
            x_edge=[min(xpoint):dx:max(xpoint)];                            %Determine extremes (dx was defined above as min_disp / 5;
            y_edge=[min(ypoint):dy:max(ypoint)];
            [X_graph,Y_graph]=meshgrid(x_edge,y_edge);
            LEP_graph=griddata(xpoint,ypoint,LEp_graph,X_graph,Y_graph);    %Interpolated strain contour
            %clf();
            %plot the data
            
            %subplot(5,2,zz)
            figure
            contourf(X_graph / Rm,Y_graph / thickness,LEP_graph,20)
            %plot3(xcoord(passed),ycoord(passed),LEp(passed),'o','MarkerEdgeColor','w','MarkerFaceColor','w')
            %box_z=[max(LEp_graph) , max(LEp_graph) , max(LEp_graph) , max(LEp_graph) , max(LEp_graph)]';
            %plot3(box_x,box_y,box_z,'w--') 
            shading flat
            %xlabel('Undeformed X (in)','Fontsize',fts)
            %ylabel('Undeformed Y (in) (axial)','Fontsize',fts)
            %zlabel('e^p','Fontsize',fts)
            %title(sprintf('e^p: Stage %d of %d',k,last),'FontSize',14)
            axis([min(x_edge)/Rm max(x_edge/Rm) -.06/thickness .06/thickness min(LEp_graph) max(LEp_graph) min(LEp_graph) max(LEp_graph)])
%      elseif k == 609;
%          openfig('Stg609_old.fig');
     end;
            cbr = colorbar('location','EastOutside','Fontsize',fts) 
            caxis(newcaxis(zz,[1 2]))
            set(cbr,'Ytick',newcaxis(zz,[3:end]));
            %caxis([.1 .8])
%              if k == 200
%                  caxis([0.05 .15])
%              elseif k == 436
%                  caxis([.12 .26])
%              elseif k == 480
%                  caxis([.13 .29])
%              elseif k == 609
%                  caxis([.2 .8])
%              end
            set (gca,'Fontsize',ftstics);                               %Axis label fontsize
            set(gcf,'Units','inches','Position',[3 3 8 8*7/9]);    %Figure size and position
            set(gcf, 'color', [1 1 1] );
            %set(gca, 'XTick', 'Color','k','YTick','color','k');
            
            %savefig(sprintf('%s\\Contour_Stg%d_Mark%d',savepath,k,zz));
            %print(gcf,'-dpng',sprintf('%s\\Contour_Stg%d_Mark%d',savepath,k,zz));
            %close
            zz = zz+1;
end
%{
figure
plot(STF(:,2),STF(:,3),'b','LineWidth',2)
hold on
for i=profStages'
    a=plot(STF(i,2),STF(i,3),'o','MarkerEdgeColor','r','MarkerFaceColor','r');
    i=i+1;
end;
b=plot(STF(profStages(1),2),STF(profStages(1),3),'o','MarkerEdgeColor','y','MarkerFaceColor','y');
hold off
axis([0 1.1*STF(end,2) 0 1.1*maxf])
xlabel('Time (sec)','Fontsize',fts) 
ylabel('Axial Force (kips)','Fontsize',fts)
mytitlestring=sprintf('TT2-%d - \\alpha=%.2f - FS%d SS%d - %s',TT2,alpha,Facet_size,Step_size,calctype);
title(mytitlestring,'Fontsize',fts)
L=legend([b a],{'Before Load Max', 'After Load Max'},'Location','Southeast');
set(L,'Fontsize',fts)
export_fig(sprintf('%s\\Time-Force.png',folder));
   %}