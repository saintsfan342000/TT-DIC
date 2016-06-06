clear;
close all;


%Modified to where time_force should include (10stage (2)time (3)Axial Sts(3)Shear sts

%Run just the final stage to find out the i-j indices of the max point.  Then pull out that whole column (same i-indices) for every stage.

% The same as Nico's code "post_matlab_3D_DIC_linear_v2", except I changed variable names to my liking, and redid a couple small things to make them consisten with the way I code

%%%%%% PRELIMINARY DATA to ENTER %%%%%% 
    last = 609;         %Last stage in which specimen is not fractures
    TT2=17;             %Expt number
    alpha=1;          %Alpha
%Specimen measurements
    Rm = 0.8942;        % Mean radius
    thickness = 0.0382;  %Wall thickness
    Lg = 0.4;         
    tube=16;
%DIC paramaters
    Facet_size = 19;    %pix                            
    Step_size = 6;      %pix                              
    calctype='Linear';
%path for the files that the PERL program produces
    PATH = 'F:\Martin_Experiments\AAA_TensionTorsion\TT2-17\AramisExport_MissingRemoved_expanded\'; %MUST CHANGE
%prefix of the name
    %prefix = sprintf('TT2-%d_DC%d_FS%d_SS%d',TT2,tube,Facet_size,Step_size);                %MUST CHANGE don't include underscore after linear
    prefix = 'TT2-17_expanded_FS19SS6_';
%create the new folder where all the exported files are
%     folder = sprintf('TT2-%d_MatlabResults',TT2);
%     mkdir(folder);
%     folder2 = sprintf('%s\\profile\\',folder);
%     mkdir(folder2);

addpath(sprintf('%s\\Matlab\\extras',PATH(1:2)));
    
%plotting parameters
    %screen dimension
        screen = 4/3;                               %MUST CHANGE
    %fontsize
        fts = 18;                       %Font size for figures
        ftstics = round(fts*3/4);       %Font size for figure axes
    if alpha >= 3.25
        Measured_output = 'Shear';                   % Column labels for export files
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


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%   1-Determination of the area to scan.    %%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%   Make the plot in which you'll draw a box %%%%%%%%%%%%%%%%%%%% 

    name = sprintf('%s%s%d.dat',PATH,prefix,last);             %name of the file
    A = load(name);          % Open the file A = matrix (nx13) with
    x=A(:,3) .* MMtoIn;     % Output file is in mm.   Converts to inches.
    y=A(:,4) .* MMtoIn;
    z=A(:,5) .* MMtoIn;

    %Creation of a grid in order to plot the data
    %Determination of the min distance between any two facets in the stage
    min_disp = 1000;    %Initialize to erroneously large value
    for k = [1:round(length(x)/30)] %i will go from 2 to length(X)-1 in 30 equally spaced increments
        for q= [1:round(length(x)/30)]    %j will cover all indices except i itself
            dist= ( (x(k)-x(q))^2 + (y(k)-y(q) )^2 + ( z(k)-z(q) )^2 )^0.5 ; %Recursive process
            if dist==0
                dist = 11000;
            end
            min_disp=min(min_disp,dist);
        end
    end
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%2-Calculate strain values in the scanned area.%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Step_size_th = (min_disp/thickness);
Facet_size_th = Facet_size / Step_size * (min_disp/thickness);

%Definition of the size of the averaging zone for grid method
size_av = Size_Scott - thickness * Facet_size_th;



%Initialize a few things before looping and calculating every stage


% Key of the different columns in the files we're reading in
%A(:,1) index i [] corresponding to the x axis. Note Aramis i,j index seems to always start in bottom left corner of facet field, though this is not he origin of the coordinate system
%A(:,2) index j [] corresponding to the y axis
%A(:,3:5) undeformed coordinates X Y Z [mm] 
%A(:,6:9) displacement U V W ur=w [mm] 
%A(:,10:13) Surface gradient transformation F(0,0) F(0,1) F(1,0) F(1,1) []

for k = 1:609;
    scot_count = 0;
    k
    clear A
    name = sprintf('%s%s%d.dat',PATH,prefix,k);             %name of the file
    A = load(name);          % Open the file A = matrix (nx13) with
    A(:,3:9) = A(:,3:9) * MMtoIn;
    MaxPt = A(A(:,1) == 196 & A(:,2) == 22,:);
    for i = 1; %<-----My results will differ from Nico's code b/c his has an error
        clear RowsInAInScotGrid F U RotU diagU LE_calc scotLEp NE_calc scot_x scot_y scotNEx scotNEy scotNExy scotgamma;
        scot_count = scot_count + 1;     %Increment
        %xcoord(i), for example, will give us the xcoord of max LEp point in the ith scanned column
        %So here we're finding all rows in A that contain computation points whose coords lie inside this inside a Scott-grid-sized box that's centered around x/y-coord(i)
        %RowsinAInScotGrid is a COLUMN VECTOR
        RowsInAInScotGrid = find( A(:,3)<=(MaxPt(3)+size_av/2) & A(:,3)>=(MaxPt(3)-size_av/2) & A(:,4)<=(MaxPt(4)+size_av/2) &  A(:,4)>=(MaxPt(4)-size_av/2)); 
            %Now calculate strain data for each of the points that are in this Scot-sized box
            for m = 1:length(RowsInAInScotGrid);
                F = [[A(RowsInAInScotGrid(m),10),A(RowsInAInScotGrid(m),11)] ;[A(RowsInAInScotGrid(m),12),A(RowsInAInScotGrid(m),13)]];     %transformation gradient F=RU
                U = transpose(F)*F;                         %stretching tensor???
                [RotU,diagU] = eig(U);                      %principal stretch
                LE_calc = 0.5 * [[log(diagU(1,1)),0];[0,log(diagU(2,2))]];  %logarithmic strain in the principal coordinate system
                scotLEp(m) = ( 2./3. * ( LE_calc(1,1)^2 + LE_calc(2,2)^2 + (-LE_calc(1,1)-LE_calc(2,2))^2 ))^0.5 ;   %Logarithmic Cumulative Plastic strain
                %R = F * U^(-0.5);                          %Rotation tensor
                NE_calc = (U^0.5 - eye([2 2]));             %Technical strain calculation in the stretching coordinate system
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
    export_Scott(k,:) = [nanmean(OurSize_Scott),length(RowsInAInScotGrid),nanmean(NEx_Scott),nanmean(NEy_Scott),abs(nanmean(gamma_Scott)),nanmean(LEp_Scott)];
    
%}
 end    %%%  ENDS THE BIG "for k=1:last" loop.  Onto the next stage
   