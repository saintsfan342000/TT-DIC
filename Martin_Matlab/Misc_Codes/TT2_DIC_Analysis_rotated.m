function TT2_DIC_Analysis_rotated;
clear;
close all;

% The same as Nico's code "post_matlab_3D_DIC_linear_v2", except I changed variable names to my liking, and redid a couple small things to make them consisten with the way I code

%%%%%% PRELIMINARY DATA to ENTER %%%%%% 
    last = 562;         %Last stage in which specimen is not fractures
    TT2=27;             %Expt number
    alpha=1.0;          %Alpha
%Specimen measurements
    Rm = 0.8938;        % Mean radius
    thickness = 0.0384;  %Wall thickness
    Lg = 0.4;         
    tube=16;
%DIC paramaters
    Facet_size = 19;    %pix                            
    Step_size = 6;      %pix                              
    calctype='Linear';

curdir=pwd;
addpath(sprintf('%s\\Matlab\\extras',curdir(1:2)));     %Adds export_fig

%path for the files that the PERL program produces
    PATH = sprintf('%s\\Martin_Experiments\\TT2-%d\\AramisExport_MissingRemoved\\',curdir(1:2),TT2); %MUST CHANGE
%prefix of the name
    prefix = sprintf('TT2-%d_DC%d_FS%d_SS%d',TT2,tube,Facet_size,Step_size);                %MUST CHANGE don't include underscore after linear
%create the new folder where all the exported files are
    folder = sprintf('TT2-%d_MatlabResults',TT2);
    mkdir(folder);
    folder2 = sprintf('%s\\profile\\',folder);
    mkdir(folder2);
        
%plotting parameters
    %screen dimension
        screen = 4/3;                               %MUST CHANGE
    %fontsize
        fts = 18;                       %Font size for figures
        ftstics = round(fts*3/4);       %Font size for figure axes
    if alpha >= 3.25
        Measured_output = 'Hoop';                   % Column labels for export files
        stsfactor = 2 * pi  * Rm^2 * thickness;     %Converts force/torque to sts
        cal = 2;                                    %1V = 2 kips.in
        ststype = '\tau (ksi)' ;                    %For DUAL PLOT axis label
    else
        Measured_output = 'Axial';
        stsfactor = 2 * pi  * Rm * thickness;       %Converts force/torque to sts
        cal = 5;                                    %1V = 5 kips
        ststype = '\sigma_Y (ksi)';                 %For DUAL PLOT axis label
    end   
%Change mm->in
    MMtoIn = 1/25.4;
%Size of the side of the Area of calculation of Strain : Scott (grid method)
    Size_Scott = 1/16; %[in]
%Watch result
watch = [1:50:last,last];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%   1-Determination of the area to scan.    %%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%   Make the plot in which you'll draw a box %%%%%%%%%%%%%%%%%%%% 

    name = sprintf('%s%s_%d.dat',PATH,prefix,last);             %name of the file
    A = load(name);          % Open the file A = matrix (nx13) with
    x=A(:,3) .* MMtoIn;     % Output file is in mm.   Converts to inches.
    y=A(:,4) .* MMtoIn;
    z=A(:,5) .* MMtoIn;

    %Calculation for each facet of the Logarithmic cumulative plastic strain
    for i = 1:length(x); %Note length(x) is simply # of points in this last stage output file
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
    %Determination of the min distance between any two facets in the stage
    min_disp = 1000;    %Initialize to erroneously large value
    for k = [2 : round(length(x)/30) : size(x,1)-1] %i will go from 2 to length(X)-1 in 30 equally spaced increments
        for q= [1:1:k-1,k+1:1:size(x,1)]    %j will cover all indices except i itself
            min_disp = min ( min_disp , ( (x(k)-x(q))^2 + (y(k)-y(q) )^2 + ( z(k)-z(q) )^2 )^0.5 ); %Recursive process
        end
    end
    dx=min_disp/5;  
    dy=min_disp/5;
    %dz=min_disp/10;
    x_edge=[min(x):dx:max(x)];
    y_edge=[min(y):dy:max(y)];
    %z_edge=[floor(min(z)):dz:ceil(max(z))];
    [Xgrd,Ygrd]=meshgrid(x_edge,y_edge);
    LEP=griddata(x,y,LEp,Xgrd,Ygrd);  %Interpolate the strain 

    close all;  %plot the data
    contourf(Xgrd,Ygrd,LEP,20)        %20 = number of contour levels
    hold on
    shading flat                %Gets rid of isolines
    xlabel('Undeformed X (in)','Fontsize',fts)
    ylabel('Undeformed Y (in) (axial)','Fontsize',fts)
    axis([min(x) max(x) min(y) max(y) min(LEp) max(LEp) min(LEp) max(LEp)])
    colorbar('location','EastOutside','Fontsize',fts) 
    set (gca,'Fontsize',ftstics)
    set(gcf,'Units','Normalized','Outerposition',[0 0 1 1]); 
    hold off

    disp('Draw a box around the zone you want to study (ie highest strain zone)');
    crop = ginput(2);
    Xmin = min(crop(:,1));  Xmax = max(crop(:,1));  Ymin = min(crop(:,2)); Ymax = max(crop(:,2)); 

%     Xmin=-0.371807014419573;   %For comparing to Nico's original code
%     Xmax=0.327031253296962;
%     Ymin=0.019659302207930;
%     Ymax=0.060041973162229;

      Xmin=  -0.158893144157042;
      Xmax=  0;
      Ymin=  0.023808206758029;
      Ymax=  0.055616474975456;

    box_x = [Xmin , Xmax , Xmax , Xmin , Xmin]';
    box_y = [Ymin , Ymin , Ymax , Ymax , Ymin]';
    box_z = [max(LEp) , max(LEp) , max(LEp) , max(LEp) , max(LEp)]';

    close all;
    contourf(Xgrd,Ygrd,LEP,20)
    hold on
    plot3(box_x,box_y,box_z,'y-o','MarkerEdgeColor','y','MarkerFaceColor','y') 
    xlabel('Undeformed X (in)','Fontsize',fts)
    ylabel('Undeformed Y (in) (Axial)','Fontsize',fts)
    zlabel('e^p','Fontsize',fts)
    mytitlestring=sprintf('TT2-%d - \\alpha=%.2f - FS%d SS%d - %s',TT2,alpha,Facet_size,Step_size,calctype);
    shading flat
    title(mytitlestring,'FontSize',14)
    axis([min(x) max(x) min(y) max(y) min(LEp) max(LEp) min(LEp) max(LEp)])
    colorbar('location','EastOutside','Fontsize',fts) 
    set (gca,'Fontsize',ftstics)    %Axis label fontsize
    set(gcf,'Units','Normalized','Outerposition',[0 0 1 1]); %Plot size and position
    hold off
    %print(gcf,'-dpng',sprintf('%s/last_picture',folder)); %save the figure
    export_fig(sprintf('%s\\last_picture.png',folder));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%2-Calculate strain values in the scanned area.%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Step_size_th = (min_disp/thickness);
Facet_size_th = Facet_size / Step_size * (min_disp/thickness);

%Definition of the size of the averaging zone for grid method
size_av = Size_Scott - thickness * Facet_size_th;

%Stage Time Force data...Create the 10 points we'll highlight in Dr K's classic figure
STF = load(sprintf('%stime_force.dat',PATH));                    %Read in and load stage time force data.
STF((last+2):length(STF),:)=[];                                  %In case I exported  more un-cleaned Aramis files than the last stage PRIOR to failure
                                                                 %Remember - last is the last STAGE NUMBER.  So length(STF) = last+1.  So we want to cut last+2 and beyond.
STF(1,:)=[];                                                     %Get rid of the first line (stage zero) so that indexing is consistent
STF(:,3)=STF(:,3)*cal/stsfactor;                                 % Convert voltage to stress
[maxf,locf] = max(STF(:,3));                                     %Need the index of the max force
prof_num = 10;                                                   %Number of profile points we'll highlight in our figures
prof_inc = round((last-locf) / (prof_num-1));                    %Determine the equally-spaced stage increment between each of these points
profStages = [locf:prof_inc:locf+(prof_num-2)*prof_inc,last]';   %profStages is the stage numbers of our 10 pts.  Column vector for printf purposes 

%Initialize a few things before looping and calculating every stage
export_max=zeros(last,7);
export_mean=zeros(last,7);
export_stdv=zeros(last,7);
export_Scott=zeros(last,7);
prof_count=1;
profLEp{prof_num}=[];

%Cycle through the stages
for k = 1 : last
    col_count=0;                                                  %Count the number of point columns that we'll calculate 
    scot_count=0;                                                     %Count the number of points in the Scott method
    clear A LEp NEx NEy NExy gamma;
    name = sprintf('%s%s_%d.dat',PATH,prefix,k);                  %Name of the stage file to open
    A = load(name);
    A(:,3:5)=A(:,3:5) * MMtoIn;      % Output file is in mm.   Converts to inches for all further calculations
    A(:,6:9)=A(:,6:9) * MMtoIn;
    %Need to find the unbroken vertical colums of points in the scan zone, unbroken meaning there aren't any missing points in each column of computation points
    %RowsInAInTheBox will give us the matlab row index in A in which there's a computation point that lies inside the box we drew.  It returns a COLUMN VECTOR of varying length.
    RowsInAInTheBox = find( A(:,3) >= Xmin & A(:,3) <= Xmax & A(:,4) >= Ymin & A(:,4) <= Ymax);
    %NICO: "Find the index i that are used to calculate the strain ie without index j missing"
    AramisXIndices=A(RowsInAInTheBox,1); %AramisXIndices=A(RowsInAInTheBox,1) will give us the aramis x-axis indices of those points that lie in the box.  COLUMN VECTOR since RowsInAInTheBox is a column
    %Loop through each unique aramis x-axis index that lies in the box
    for uniqueX = unique(AramisXIndices)';                  %<-- transpose because this must be a row vector for the for loop to incrememnt correctly!
        uniqeXloc = find( AramisXIndices == uniqueX );      %uniqueXloc is location in AramisXIndices, which is also the row in A, in which the particular value of uniqueX lies.  Usually a column vector.
        AramisYIndices = A(RowsInAInTheBox(uniqeXloc),2);   %AramisYIndices is a vector of all the different Aramis y-indexes (column 2 of A) that are paired with the particular uniqueX. Same dim as uniqueXloc.
        %This checks whether the y-indexes in AramisYIndices are contiguous (that is, if no index number is skipped), so if this "column" of particular Aramis x indices is unbroken then we'll calculate its points and "scan 
        if (length(AramisYIndices) == max(AramisYIndices) - min(AramisYIndices) +1)
            clear colLEp colNEx colNEy colNExy colgamma NE_calc_rot NE_calc R U diagU LE_calc;
            col_count = col_count+1;
            %Initialize / pre-allocate length of the vectors we'll create in the for loop below
                colLEp=zeros(1,length(AramisYIndices)); colNEx=zeros(1,length(AramisYIndices)); colNEy=zeros(1,length(AramisYIndices)); 
                colNExy=zeros(1,length(AramisYIndices)); colgamma=zeros(1,length(AramisYIndices)); 
            % Calculate strain for every point in this column
            for m = 1 : length(AramisYIndices);
                F = [[A(RowsInAInTheBox(uniqeXloc(m)),10),A(RowsInAInTheBox(uniqeXloc(m)),11)] ;[A(RowsInAInTheBox(uniqeXloc(m)),12),A(RowsInAInTheBox(uniqeXloc(m)),13)]]; %transformation gradient F=RU
                U = transpose(F)*F; %stretching tensor???
                [~,diagU] = eig(U); %principal stretch
                LE_calc = 0.5 * [[log(diagU(1,1)),0];[0,log(diagU(2,2))]];  %logarithmic strain in the principal coordinate system
                colLEp(m) = ( 2./3. * ( LE_calc(1,1)^2 + LE_calc(2,2)^2 + (-LE_calc(1,1)-LE_calc(2,2))^2 ))^0.5;   %Logarithmic Cumulative Plastic strain
                R = F * U^(-0.5);  %Rotation tensor
                NE_calc_rot= (U^0.5 - eye([2 2]));     
                NE_calc=R*NE_calc_rot*R.';
                %Technical strain in the x',y' material coordinate system
                colNEx(m) = NE_calc(1,1);       %% col prefix implies column, as in NEx for each point in this column currently being calculated
                colNEy(m) = NE_calc(2,2);
                colNExy(m) = NE_calc(1,2);
                colgamma(m) = atan(colNExy(m)/(1+colNEx(m))) + atan(colNExy(m)/(1+colNEy(m)));
                col_x(m) = A(RowsInAInTheBox(uniqeXloc(m)),3);    
                col_y(m) = A(RowsInAInTheBox(uniqeXloc(m)),4);
            end 
        [LEp(col_count) , locLEp] = max(colLEp);  %Max LEp in the current colum (the one that passed the if test) of computation points and its location in AramisYIndices
        NEx(col_count) = colNEx(locLEp);          %Store the other technical strain data at the point with max LEp
        NEy(col_count) = colNEy(locLEp);          %Notice maxLEp and things below it are being kept as ROW VECTORS, indexing with col_count
        gamma(col_count) = colgamma(locLEp);      %So all of the strain data and coordinates of the max LEp point in each column will be saved
        xcoord(col_count) = col_x(locLEp);        %After all iterations of uniqueX each of these will be a row vector with the same length (though this won't necessarily be the same length as, say, AramisXIndices since we're only scanning columns whose j-indices are unbroken
        ycoord(col_count) = col_y(locLEp);
        end                                         %End of "if column of X-index is unbroken"
    end                                             %End of  "for uniqueX=unique(AramisXIndices)'"
 
    
    %%%% Keep strain data for only those points satisfying criteria
    %%%% Criteria is to omit points whose NEy/gamma ratio doesn't fall within half of a StdDev of the mean of this ratio
    ratio = NEy ./ gamma;                           %Find the ratio of the NEy/gamma each scanned column
    ratioAvg = nanmean(ratio);                      %Avg ratio
    ratioSDEV = nanstd(ratio);                      %StdDev of ratio
    passed = find(ratio >= ratioAvg - 0.5 * ratioSDEV & ratio <= ratioAvg + 0.5 * ratioSDEV);    %Locations in LEp, NEx, etc vectors of the points that passed.  A row vector.
    LEp_pass=LEp(passed); NEx_pass=NEx(passed); NEy_pass=NEy(passed); gamma_pass=gamma(passed); xcoord_pass=xcoord(passed); ycoord_pass=ycoord(passed); %Keep just those that passed.
    [~, index_max] = max(LEp_pass);                 %Find where the max LEp_pass was and export.  These will end up being last x 7 matrices
    export_max(k,:) = [length(passed), STF(k,2), STF(k,3),NEx_pass(index_max),NEy_pass(index_max),abs(gamma_pass(index_max)),LEp_pass(index_max)];
    export_mean(k,:) = [length(passed), STF(k,2),STF(k,3),nanmean(NEx_pass),nanmean(NEy_pass),abs(nanmean(gamma_pass)),nanmean(LEp_pass)];
    export_stdv(k,:) = [length(passed), STF(k,2),STF(k,3),nanstd(NEx_pass),nanstd(NEy_pass),abs(nanstd(gamma_pass)),nanstd(LEp_pass)];

%%%% Profile of ep for different data_force_points....Still in the loop for "k=1:last" %%%%
    if any(profStages==k) == 1       %Checks to see if the current stage being analyzed (k) is one of our profStages
        profMaxlocA = find( A(:,3)  == xcoord_pass(index_max) & A(:,4)  == ycoord_pass(index_max)); %What row in A is the maximum strain point in this stage located?
        profColum = find( A(:,1) == A(profMaxlocA,1));  %Find all the other points that share this point's Aramis X index (column vector)
        for j = 1:length(profColum) 
            %Create a cell array, each array being a profCount stage. Each array has j rows, j=# points in the Aramis X-index column
            %COLUMN 1 - Calculate Y coordinate layer :  Undef y-coord + V displacement normalized by wall thickness
            profLEp{prof_count}(j,1) = (A(profColum(j),4)+A(profColum(j),7))/thickness;     
            F = [[A(profColum(j),10),A(profColum(j),11)] ;[A(profColum(j),12),A(profColum(j),13)]]; %transformation gradient F=RU
            U = transpose(F)*F;             %stretching tensor???
            diagU = eig(U);                 %principal stretch
            LE_calc = 0.5 * log(diagU);     %logarithmic strain in the principal coordinate system
            %COLUMN 2 - Log Cum Plastic strain Column 2
            profLEp{prof_count}(j,2) = ( 2./3. * (  LE_calc(1)^2 + LE_calc(2)^2 + (-LE_calc(1)-LE_calc(2))^2) )^0.5 ;      
        end
        prof_count=prof_count+1;        %Incrememnt prof_count
    end  

%%%% Calculation of the max equivalent strain and so on using an average value
%over an area equivalent of the one of Scott (grid method)
% Appears to operate over all previously scanned columns and points, not just those which passed the ratio test (but this is acceptable as there's no guarantee that the Scot-grid box around those points which passed will also have passed)
    for i = 1:length(xcoord); %<-----My results will differ from Nico's code b/c his has an error
        clear RowsInAInScotGrid F U RotU diagU LE_calc scotLEp NE_calc_rot NE_calc scot_x scot_y scotNEx scotNEy scotNExy scotgamma;
        scot_count = scot_count + 1;     %Increment
        %xcoord(i), for example, will give us the xcoord of max LEp point in the ith scanned column
        %So here we're finding all rows in A that contain computation points whose coords lie inside this inside a Scott-grid-sized box that's centered around x/y-coord(i)
        %RowsinAInScotGrid is a COLUMN VECTOR
        RowsInAInScotGrid = find( A(:,3)<=(xcoord(i)+size_av/2) & A(:,3)>=(xcoord(i)-size_av/2) & A(:,4)<=(ycoord(i)+size_av/2) &  A(:,4)>=(ycoord(i)-size_av/2)); 
            %Now calculate strain data for each of the points that are in this Scot-sized box
            for m = 1:length(RowsInAInScotGrid);
                F = [[A(RowsInAInScotGrid(m),10),A(RowsInAInScotGrid(m),11)] ;[A(RowsInAInScotGrid(m),12),A(RowsInAInScotGrid(m),13)]];     %transformation gradient F=RU
                U = transpose(F)*F;                         %stretching tensor???
                [RotU,diagU] = eig(U);                      %principal stretch
                LE_calc = 0.5 * [[log(diagU(1,1)),0];[0,log(diagU(2,2))]];  %logarithmic strain in the principal coordinate system
                scotLEp(m) = ( 2./3. * ( LE_calc(1,1)^2 + LE_calc(2,2)^2 + (-LE_calc(1,1)-LE_calc(2,2))^2 ))^0.5 ;   %Logarithmic Cumulative Plastic strain
                R = F * U^(-0.5);                          %Rotation tensor
                NE_calc_rot = (U^0.5 - eye([2 2]));             %Technical strain calculation in the stretching coordinate system
                NE_calc=R*NE_calc_rot*R.';
                scot_x(m) = A(RowsInAInScotGrid(m),3);      %Technical strain in the x',y' material coordinate system
                scot_y(m) = A(RowsInAInScotGrid(m),4);
                scotNEx(m) = NE_calc(1,1);
                scotNEy(m) = NE_calc(2,2);
                scotNExy(m) = NE_calc(1,2);
                scotgamma(m) = atan(scotNExy(m)/(1+scotNEx(m))) + atan(scotNExy(m)/(1+scotNEy(m)));
            end
            NEx_Scott(scot_count) = nanmean(scotNEx); %For this one DIC point, we average all the m Scott-grid points that were near it and store this average strain data
            NEy_Scott(scot_count) = nanmean(scotNEy);   %scot_count(k) will be equal to length(xcoord), but not the same as col_count. Especialy true in later stages when facets drop out, so fewer columns will be contiguous
            gamma_Scott(scot_count) = nanmean(scotgamma);
            LEp_Scott(scot_count) = nanmean(scotLEp);
            OurSize_Scott(scot_count) = Facet_size_th * thickness + mean( [max(scot_x) - min(scot_x) , max(scot_y) - min(scot_y)] );  %Calculate our Scott-grid size based on the points we're averaging over
    end   
    export_Scott(k,:) = [nanmean(OurSize_Scott), STF(k,2), STF(k,3),nanmean(NEx_Scott),nanmean(NEy_Scott),abs(nanmean(gamma_Scott)),nanmean(LEp_Scott)];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%     
     %This doesn't produce any data...it just makes the strain contour every 50 stages
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%{        
    if any(watch==k) == 1
            close all
            clear xpoint ypoint F_graph U_graph RotU_graph diagU_graph LE_calc_graph LEp_graph x_edge y_edge X_graph Y_graph LEP_graph box_z;
            xpoint=(A(:,3) ); % [in]                                %Coordinates of all points
            ypoint=(A(:,4)) ; % [in]
            for i = 1:length(xpoint)                                        %Calculation for each facet of the Logarithmic cumulative plastic strain
                F_graph = [[A(i,10),A(i,11)] ;[A(i,12),A(i,13)]];
                U_graph = transpose(F_graph)*F_graph;                       %stretching tensor???
                [RotU_graph,diagU_graph] = eig(U_graph);                    %principal stretch
                LE_calc_graph = 0.5 * [[log(diagU_graph(1,1)),0];[0,log(diagU_graph(2,2))]];    %logarithmic strain in the principal coordinate system
                LEp_graph(i) = ( 2./3. * ( LE_calc_graph(1,1)^2 + LE_calc_graph(2,2)^2 + (-LE_calc_graph(1,1)-LE_calc_graph(2,2))^2 ))^0.5 ;    %Logarithmic Cumulative Plastic strain
            end
            x_edge=[min(xpoint):dx:max(xpoint)];                            %Determine extremes (dx was defined above as min_disp / 5;
            y_edge=[min(ypoint):dy:max(ypoint)];
            [X_graph,Y_graph]=meshgrid(x_edge,y_edge);
            LEP_graph=griddata(xpoint,ypoint,LEp_graph,X_graph,Y_graph);    %Interpolated strain contour
            clf();
            %plot the data
            contourf(X_graph,Y_graph,LEP_graph,20)
            hold on
            plot3(xcoord(passed),ycoord(passed),LEp(passed),'o','MarkerEdgeColor','w','MarkerFaceColor','w')
            box_z=[max(LEp_graph) , max(LEp_graph) , max(LEp_graph) , max(LEp_graph) , max(LEp_graph)]';
            plot3(box_x,box_y,box_z,'w--') 
            shading flat
            xlabel('Undeformed X (in)','Fontsize',fts)
            ylabel('Undeformed Y (in) (axial)','Fontsize',fts)
            zlabel('e^p','Fontsize',fts)
            title(sprintf('e^p: Stage %d of %d',k,last),'FontSize',14)
            axis([min(x) max(x) min(y) max(y) min(LEp_graph) max(LEp_graph) min(LEp_graph) max(LEp_graph)])
            colorbar('location','EastOutside','Fontsize',fts) 
            set (gca,'Fontsize',ftstics);                               %Axis label fontsize
            set(gcf,'Units','Normalized','Outerposition',[0 0 1 1]);    %Figure size and position
            hold off
            export_fig(sprintf('%s\\StnStage%d.png',folder,k));
            %pause(2);
    end
%}
 end    %%%  ENDS THE BIG "for k=1:last" loop.  Onto the next stage
   

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%3-PLOT DATA%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Remember:  export_mean, export_max, etc. are those which passed the ratio test
%%export_max(k,:) = [length(passed), STF(k,2), cal * STF(k,3),NEx_pass(index_max),NEy_pass(index_max),abs(gamma_pass(index_max)),LEp_pass(index_max)];


% Epsilon-Y vs Gamma 
figure
%(:,6) is gamma_pass, (:,5) is NEy_pass
plot(export_mean(:,6),export_mean(:,5),'bx');       
hold on
plot(export_max(:,6),export_max(:,5),'ro');
plot(export_Scott(:,6),export_Scott(:,5),'k+');
mytitlestring=sprintf('Strain in  Rotated Local Coordinate System \n TT2-%d - \\alpha=%.2f - FS%d SS%d - %s \n Averaging Size for Scott-grid Results = %f in',TT2,alpha,Facet_size,Step_size,calctype,OurSize_Scott(scot_count));
title(mytitlestring,'FontSize',12);
legend('Mean','Max', 'Average Value (Scott Grid)' , 'Location','SouthWest','Fontsize',2)
xlabel('\gamma (rad)','Fontsize',fts)
ylabel('\epsilon _y ','Fontsize',fts)
set (gca,'Xgrid','on','Ygrid','on','LineWidth',2)
set (gca,'Fontsize',ftstics)
set(gcf,'Units','Normalized','Outerposition',[0 0 1/screen 1]);
hold off 
%print(gcf,'-dpng',sprintf('%s/gamma_epsy',folder));
export_fig(sprintf('%s\\gamma_epsy.png',folder));

% e^p vs time
figure
%(:,2) is time, (:,7) is LEp_pass
errorbar(export_mean(:,2),export_mean(:,7),0.5*export_stdv(:,7));
hold on
plot(export_max(:,2),export_max(:,7),'ro')
plot(export_Scott(:,2),export_Scott(:,7),'k+');
mytitlestring=sprintf('TT2-%d - \\alpha=%.2f - FS%d SS%d - %s \n Averaging size for like-Scott results = %f in',TT2,alpha,Facet_size,Step_size,calctype,OurSize_Scott(scot_count));
title(mytitlestring,'FontSize',12)
legend('Mean','Max', 'Average Value (Scott Grid)' ,'Location','NorthWest')
axis([0 max(export_max(:,2)) 0 1.1*max(export_max(:,7))])
xlabel('Time (sec)','Fontsize',fts)
ylabel('e^p','Fontsize',fts)
set (gca,'Xgrid','on','Ygrid','on','LineWidth',2)
set (gca,'Fontsize',ftstics)
set(gcf,'Units','Normalized','Outerposition',[0 0 1/screen 1]);
hold off  
%print(gcf,'-dpng',sprintf('%s/time_ep.png',folder));  
export_fig(sprintf('%s\\time_ep.png',folder));

% Prepare marker types and colors for profile figures
marker_type = ['+';'o';'*';'x'];
marker_color = ['r';'g';'b';'m';'k'];
marker_style(prof_num,:) = '  ';
m0=1; n0=1; i=1;  %Initialize
while i <= prof_num %prof_num is the number of profile stages we're highlighting (10)
    marker_style(i,:) = [marker_color(m0) marker_type(n0)]; 
    i=i+1;
    if m0<length(marker_color)
        m0=m0+1;
    else
        m0=1;
    end;
    if n0<length(marker_type)
        n0=n0+1;
    else
        n0=1;
    end;
end

% DUAL PLOT:  Stress vs time and Y vs e^p profile
figure
set(gcf,'Units','Normalized','Outerposition',[0 0 1 1]); %Makes the figure full screen
subplot(1,2,1);
% Stress vs time
plot([0;STF(:,2)],[0;STF(:,3)],'LineWidth',2);
mytitlestring=sprintf('TT2-%d - \\alpha=%.2f - FS%d SS%d - %s',TT2,alpha,Facet_size,Step_size,calctype);
title(mytitlestring,'FontSize',12)
hold on
for i=1:prof_num
    %%%NICO's code contains an error.  We must use STF(profStages(i)+1,#)
    %%%because STF contains the stage zero in spot one, but k only cycles
    %%%from 1 to last.  Not a *huge* deal, but what we've reported as sts
    %%%at final stage isn't actuall that.  It's the 2nd to final.  Easier
    %%%solution would be to chop STF(0,:) snd get rid of all these +1s
    plot(STF(profStages(i),2),STF(profStages(i),3),marker_style(i,:),'Linewidth',2,'MarkerSize',10);
end
axis([0 1.1*max(STF(:,2)) 0 1.1*max(STF(:,3))])
xlabel('Time (sec)','Fontsize',fts)
ylabel(ststype,'Fontsize',fts)
set (gca,'Xgrid','on','Ygrid','on','LineWidth',2)
set (gca,'Fontsize',ftstics)
hold off

subplot(1,2,2);
%Y coord vs e^p profiles
for i = 1 : prof_num
    hold on
    plot(profLEp{i}(:,2),profLEp{i}(:,1),marker_style(i,:),'Linewidth',2,'MarkerSize',10);
end
xlabel('e^p','Fontsize',fts)
ylabel('Axial Coordinate (normalized by thickness)','Fontsize',fts)
set (gca,'Xgrid','on','Ygrid','on','LineWidth',2)
set (gca,'Fontsize',ftstics)
hold off
%print(gcf,'-dpng',sprintf('%s/ep_profile',folder))
export_fig(sprintf('%s\\ep_profile.png',folder));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%4-EXPORT DATA%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

output = fopen(sprintf('%s/max.dat',folder),'w');       %MAX
fprintf(output,'%% nb of points [] Time [s] %s Nominal Stress [ksi] Technical NEx [] NEy [] Gamma [rad] True cumulative strain [] \n',Measured_output);
fprintf(output, '%d %f %f %f %f %f %f\n', export_max');
fclose(output);

output = fopen(sprintf('%s/mean.dat',folder),'w');      %MEAN
fprintf(output,'%% nb of points [] Time [s] %s Nominal Stress [ksi] Technical NEx [] NEy [] Gamma [rad] True cumulative strain [] \n',Measured_output);
fprintf(output, '%d %f %f %f %f %f %f\n', export_mean');
fclose(output);

output = fopen(sprintf('%s/stdv.dat',folder),'w');      %STANDARD DEVIATION
fprintf(output,'%% nb of points [] Time [s] %s Nominal Stress [ksi] Technical NEx [] NEy [] Gamma [rad] True cumulative strain [] \n',Measured_output);
fprintf(output, '%d %f %f %f %f %f %f\n', export_stdv');
fclose(output);

output = fopen(sprintf('%s/like_Scott.dat',folder),'w');    %SCOT-SIZED GRID
fprintf(output,'%% size of the averaging zone for like-Scott [in] Time [s] %s Nominal Stress [ksi] Technical NEx [] NEy [] Gamma [rad] True cumulative strain [] \n',Measured_output);
fprintf(output, '%d %f %f %f %f %f %f\n', export_Scott');
fclose(output);

for i=1:prof_num                                        %PROFILE DATA
  clear out;
  name = sprintf('%s/export_Y_LEp_time_%6.3f_stress_%6.3f.dat',folder2,STF(profStages(i),2),STF(profStages(i),3));
  output = fopen(name,'w');
  fprintf(output,'%s Y/th [] True Cumulative plastic strain [] \n','%');
  out = profLEp{i}(:,:);
  fprintf(output,'%f %f',out');
  save(name,'out','-ASCII');
end
  
%Record facet size and step size compared to the thicknes of the specimen
fileID = fopen(sprintf('%s/facet_step.txt',folder),'w');
fprintf(fileID,'Facet size/thickness = %1.4f \n Step size/thickness = %1.4f',Facet_size_th,Step_size_th);
fclose(fileID);   

 
% Key of the different columns in the files we're reading in
%A(:,1) index i [] corresponding to the x axis. Note Aramis i,j index seems to always start in bottom left corner of facet field, though this is not he origin of the coordinate system
%A(:,2) index j [] corresponding to the y axis
%A(:,3:5) undeformed coordinates X Y Z [mm] 
%A(:,6:9) displacement U V W ur=w [mm] 
%A(:,10:13) Surface gradient transformation F(0,0) F(0,1) F(1,0) F(1,1) []

fclose all;