function here = TT2_DIC_Analysis;
clear;
close all;

%Modified to where time_force should include (10stage (2)time (3)Axial Sts(3)Shear sts

%Run just the final stage to find out the i-j indices of the max point.  Then pull out that whole column (same i-indices) for every stage.

% The same as Nico's code "post_matlab_3D_DIC_linear_v2", except I changed variable names to my liking, and redid a couple small things to make them consisten with the way I code

%%%%%% PRELIMINARY DATA to ENTER %%%%%% 
    last = 913;         %Last stage in which specimen is not fractures
    TT2=12;             %Expt number
    alpha=1.5;          %Alpha
%Specimen measurements
    Rm = 0.8940;        % Mean radius
    thickness = 0.0384;  %Wall thickness
    Lg = 0.4;         
    tube=16;
%DIC paramaters
    Facet_size = 19;    %pix                            
    Step_size = 6;      %pix                              
    calctype='Linear';
%path for the files that the PERL program produces
    PATH = 'F:\Martin_Experiments\AAA_TensionTorsion\TT2-12\AramisExport_MissingRemoved\'; %MUST CHANGE
%prefix of the name
    %prefix = sprintf('TT2-%d_DC%d_FS%d_SS%d',TT2,tube,Facet_size,Step_size);                %MUST CHANGE don't include underscore after linear
    prefix = 'TT2_21_DC15_FS19_SS6_';

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
%Watch result
watch = [1:50:last,last];
Xmin = -10; Xmax = 10; Ymax = 10; Ymin = -10;
%Stage Time Force data...Create the 10 points we'll highlight in Dr K's classic figure
STF = load(sprintf('%stime_force_new.dat',PATH));                    %Read in and load stage time force data.
STF((last+2):length(STF),:)=[];                                  %In case I exported  more un-cleaned Aramis files than the last stage PRIOR to failure
                                                                 %Remember - last is the last STAGE NUMBER.  So length(STF) = last+1.  So we want to cut last+2 and beyond.
STF(1,:)=[];                                                     %Get rid of the first line (stage zero) so that indexing is consistent
%STF(:,3)=STF(:,3)*cal/stsfactor;                                 % Convert voltage to stress
[maxf,locf] = max(STF(:,3));                                     %Need the index of the max force
prof_num = 10;                                                   %Number of profile points we'll highlight in our figures
prof_inc = round((last-locf) / (prof_num-1));                    %Determine the equally-spaced stage increment between each of these points
profStages = [locf:prof_inc:locf+(prof_num-2)*prof_inc,last]';   %profStages is the stage numbers of our 10 pts.  Column vector for printf purposes 

%Initialize a few things before looping and calculating every stage
export_max=zeros(last,8);
export_mean=zeros(last,8);
export_stdv=zeros(last,8);
export_Scott=zeros(last,8);
prof_count=1;
profLEp{prof_num}=[];
stgs = [];
here = [];
%Cycle through the stages

for k = [last-20:last]

            col_count=0;                                                  %Count the number of point columns that we'll calculate 
    scot_count=0;                                                     %Count the number of points in the Scott method
    clear A LEp NEx NEy NExy gamma;
    name = sprintf('%s%s%d.dat',PATH,prefix,k);                  %Name of the stage file to open
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
            clear colLEp colNEx colNEy colNExy colgamma;
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
                %R = F * U^(-0.5);  %Rotation tensor
                %NE_calc = R.'*(U^0.5 - eye([2 2]))*R;
                NE_calc = (U^0.5 - eye([2 2]));     
                %Technical strain in the x',y' material coordinate system
                colNEx(m) = NE_calc(1,1);       %% col prefix implies column, as in NEx for each point in this column currently being calculated
                colNEy(m) = NE_calc(2,2);
                colNExy(m) = NE_calc(1,2);
                colgamma(m) = atan(colNExy(m)/(1+colNEx(m))) + atan(colNExy(m)/(1+colNEy(m)));
                col_x(m) = A(RowsInAInTheBox(uniqeXloc(m)),3);    
                col_y(m) = A(RowsInAInTheBox(uniqeXloc(m)),4);
                col_I(m) = A(RowsInAInTheBox(uniqeXloc(m)),1);    
                col_J(m) = A(RowsInAInTheBox(uniqeXloc(m)),2);
            end 
        [LEp(col_count) , locLEp] = max(colLEp);  %Max LEp in the current colum (the one that passed the if test) of computation points and its location in AramisYIndices
        NEx(col_count) = colNEx(locLEp);          %Store the other technical strain data at the point with max LEp
        NEy(col_count) = colNEy(locLEp);          %Notice maxLEp and things below it are being kept as ROW VECTORS, indexing with col_count
        gamma(col_count) = colgamma(locLEp);      %So all of the strain data and coordinates of the max LEp point in each column will be saved
        xcoord(col_count) = col_x(locLEp);        %After all iterations of uniqueX each of these will be a row vector with the same length (though this won't necessarily be the same length as, say, AramisXIndices since we're only scanning columns whose j-indices are unbroken
        ycoord(col_count) = col_y(locLEp);
        aramI(col_count) = col_I(locLEp);        %After all iterations of uniqueX each of these will be a row vector with the same length (though this won't necessarily be the same length as, say, AramisXIndices since we're only scanning columns whose j-indices are unbroken
        aramJ(col_count) = col_J(locLEp);
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
     LEp_pass=LEp(passed);
    aramI_pass=aramI(passed); 
    aramJ_pass=aramJ(passed); %Keep just those that passed.
    if any(17 == aramI_pass) && any(21 == aramJ_pass)
        stgs = [stgs;k];
    end
end
   
   for z = 1:length(stgs);
    col_count=0;                                                  %Count the number of point columns that we'll calculate 
    scot_count=0;                                                     %Count the number of points in the Scott method
    clear A LEp NEx NEy NExy gamma;
    name = sprintf('%s%s%d.dat',PATH,prefix,stgs(z));                  %Name of the stage file to open
    A = load(name);
    A(:,3:5)=A(:,3:5) * MMtoIn;      % Output file is in mm.   Converts to inches for all further calculations
    A(:,6:9)=A(:,6:9) * MMtoIn;
    locs = find( A(:,1)== 17 & A(:,2)==21);
    locs;
    if ~isempty(locs)
        A = A(locs,:);
        F = [A(10) A(11);A(12) A(13)];
        U = sqrtm(F'*F);
        [~,diagU] = eig(U); %principal stretch
        LE_calc = [[log(diagU(1,1)),0];[0,log(diagU(2,2))]];  %logarithmic strain in the principal coordinate system
        colLEp = ( 2./3. * ( LE_calc(1,1)^2 + LE_calc(2,2)^2 + (-LE_calc(1,1)-LE_calc(2,2))^2 ))^0.5;   %Logarithmic Cumulative Plastic strain
        here = [here; stgs(z) locs colLEp];
    end
end;
    
end
