function TT2_DIC_Analysis_differentStnDef;
clear;
close all;

%Specialty script for making epsilon-v-gamma plots
% This works a lot like Nico's original script, but it employs Dr.K's definition of epsilon and gamma
% And exports eps and gamma for (1) The mean of all points and (2)The max LEp point in each stage (not the same point throughout, necessarily)

key = xlsread('F:\Martin_Experiments\AAA_TensionTorsion\TT-Summary.xlsx');
key = key(:,[1 end]);

exp = [32];
for G = 1:length(exp)
close all
last = key(key(:,1)==exp(G),2);

%path for the files that the PERL program produces
    PATH = sprintf('F:\\Martin_Experiments\\AAA_TensionTorsion\\TT2-%d\\AramisExport_MissingRemoved\\',exp(G));
    savepath = sprintf('F:\\Martin_Experiments\\AAA_TensionTorsion\\TT2-%d\\TT2-%d_MatlabResults\\',exp(G),exp(G)); 
    

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%   1-Determination of the area to scan.    %%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%   Make the plot in which you'll draw a box %%%%%%%%%%%%%%%%%%%% 

if exist([savepath 'box_limts.dat'])==2
    fid = fopen([savepath 'box_limts.dat'],'r')
    d = cell2mat(textscan(fid,'%f %f %f %f','delimiter','headerlines',1))
    Xmin=d(1); Xmax=d(2); Ymin=d(3); Ymax =d(4);
else

    prefix = ls(sprintf('%s*%d.dat',PATH,last))
    name = sprintf('%s%s',PATH,prefix);             %name of the file
    A = load(name);          % Open the file A = matrix (nx13) with
    x=A(:,3) ;     % Output file is in mm.   Converts to inches.
    y=A(:,4) ;
    z=A(:,5) ;

    %Calculation for each facet of the Logarithmic cumulative plastic strain
    for i = 1:length(x); %Note length(x) is simply # of points in this last stage output file
        F = [[A(i,10),A(i,11)] ;[A(i,12),A(i,13)]];
        %stretching tensor???
        U = transpose(F)*F;
        %principal stretch
        [~,diagU] = eig(U); 
        %logarithmic strain in the principal coordinate system
        LE_calc = 0.5 * [[log(diagU(1,1)),0];[0,log(diagU(2,2))]];
        %Logarithmic Cumulative Plastic strain
        LEp(i) = ( 2./3. * ( LE_calc(1,1)^2 + LE_calc(2,2)^2 + (-LE_calc(1,1)-LE_calc(2,2))^2 ))^0.5 ;
    end

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
    axis([min(x) max(x) min(y) max(y) min(LEp) max(LEp) min(LEp) max(LEp)])
    set(gcf,'Units','Normalized','Outerposition',[0 0 1 1]); 
    hold off

    disp('Draw a box around the zone you want to study (ie highest strain zone)');
    crop = ginput(2)
    
    
    Xmin = min(crop(:,1));  Xmax = max(crop(:,1));  Ymin = min(crop(:,2)); Ymax = max(crop(:,2)); 

    box_x = [Xmin , Xmax , Xmax , Xmin , Xmin]';
    box_y = [Ymin , Ymin , Ymax , Ymax , Ymin]';
    box_z = [max(LEp) , max(LEp) , max(LEp) , max(LEp) , max(LEp)]';
    close
end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%2-Calculate strain values in the scanned area.%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Stage Time Force data...Create the 10 points we'll highlight in Dr K's classic figure
STF = load(sprintf('%stime_force.dat',PATH));                    %Read in and load stage time force data.
STF((last+2):length(STF),:)=[];                                  %In case I exported  more un-cleaned Aramis files than the last stage PRIOR to failure
                                                                 %Remember - last is the last STAGE NUMBER.  So length(STF) = last+1.  So we want to cut last+2 and beyond.
STF(1,:)=[];                                                     %Get rid of the first line (stage zero) so that indexing is consistent

%Initialize a few things before looping and calculating every stage
M=[];
m=[];



%Cycle through the stages
  
if exp(G) == 31;
    itarray = [1:628 630:last];
else
    itarray = [1:last];
end;

for k = itarray
    k
    col_count=0;                                                  %Count the number of point columns that we'll calculate 
    scot_count=0;                                                     %Count the number of points in the Scott method
    clear prefix A LEp NEx NEy NExy gamma;
    if exp(G) == 31
        prefix = ls(sprintf('%s*stg%d.dat',PATH,k));
    elseif exp(G) == 30
        prefix = ls(sprintf('%s*SS6%d.dat',PATH,k));
    else
        prefix = ls(sprintf('%s*_%d.dat',PATH,k));
    end;
    name = sprintf('%s%s',PATH,prefix);             %name of the file
    A = load(name);
    A(:,3:5)=A(:,3:5);      % Output file is in mm.   Converts to inches for all further calculations
    A(:,6:9)=A(:,6:9);
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
            for q = 1 : length(AramisYIndices);
                F = [[A(RowsInAInTheBox(uniqeXloc(q)),10),A(RowsInAInTheBox(uniqeXloc(q)),11)] ;[A(RowsInAInTheBox(uniqeXloc(q)),12),A(RowsInAInTheBox(uniqeXloc(q)),13)]]; %transformation gradient F=RU
                U = transpose(F)*F; %stretching tensor???
                [~,diagU] = eig(U); %principal stretch
                LE_calc = 0.5 * [[log(diagU(1,1)),0];[0,log(diagU(2,2))]];  %logarithmic strain in the principal coordinate system
                colLEp(q) = ( 2./3. * ( LE_calc(1,1)^2 + LE_calc(2,2)^2 + (-LE_calc(1,1)-LE_calc(2,2))^2 ))^0.5;   %Logarithmic Cumulative Plastic strain
                %R = F * U^(-0.5);  %Rotation tensor
                %NE_calc = R.'*(U^0.5 - eye([2 2]))*R;
                NE_calc = (U^0.5 - eye([2 2]));     
                %Technical strain in the x',y' material coordinate system
                colNEx(q) = NE_calc(1,1);       %% col prefix implies column, as in NEx for each point in this column currently being calculated
                colNEy(q) = NE_calc(2,2);
                colNExy(q) = NE_calc(1,2);
                colgamma(q) = atan(colNExy(q)/(1+colNEx(q))) + atan(colNExy(q)/(1+colNEy(q)));
                col_NEy_NEW(q)=F(2,2)-1;
                col_gamm_NEW(q)=atan(F(1,2)/F(2,2));
            end 
        [LEp(col_count) , locLEp] = max(colLEp);  %Max LEp in the current colum (the one that passed the if test) of computation points and its location in AramisYIndices
        NEy(col_count) = colNEy(locLEp);          %Notice maxLEp and things below it are being kept as ROW VECTORS, indexing with col_count
        gamma(col_count) = colgamma(locLEp);      %So all of the strain data and coordinates of the max LEp point in each column will be saved
        NEy_NEW(col_count)=col_NEy_NEW(locLEp);
        Gamma_NEW(col_count)=col_gamm_NEW(locLEp);
        end                                         %End of "if column of X-index is unbroken"
    end                                             %End of  "for uniqueX=unique(AramisXIndices)'"
    %%%% Keep strain data for only those points satisfying criteria
    %%%% Criteria is to omit points whose NEy/gamma ratio doesn't fall within half of a StdDev of the mean of this ratio
    
    ratio = NEy ./ gamma;                           %Find the ratio of the NEy/gamma each scanned column
    ratioAvg = nanmean(ratio);                      %Avg ratio
    ratioSDEV = nanstd(ratio);                      %StdDev of ratio
    passed = find(ratio >= ratioAvg - 0.5 * ratioSDEV & ratio <= ratioAvg + 0.5 * ratioSDEV);    %Locations in LEp, NEx, etc vectors of the points that passed.  A row vector.
    
    LEp_pass=LEp(passed); NEy_pass=NEy(passed); gamma_pass=gamma(passed); %Keep just those that passed.
    NEy_NEW_pass=NEy_NEW(passed);Gamma_NEW_pass=Gamma_NEW(passed);
    [~, index_max] = max(LEp_pass);                 %Find where the max LEp_pass was and export.  These will end up being last x 7 matrices
    
    M =  [M; abs(gamma_pass(index_max)) NEy_pass(index_max) abs(Gamma_NEW_pass(index_max)) NEy_NEW_pass(index_max) ];
    m = [m; abs(nanmean(gamma_pass)) nanmean(NEy_pass) abs(nanmean(Gamma_NEW_pass)) nanmean(NEy_NEW_pass)];
end

save(sprintf('%seps-gam_mean.dat',savepath),'m','-ascii')
save(sprintf('%seps-gam_max.dat',savepath),'M','-ascii')

figure
plot(m(:,1),m(:,2),'k:')
hold on
plot(m(:,3),m(:,4),'k')
end

%{
 
d30=load('eps-gam-mean-30.dat');
d24=load('eps-gam-mean-24.dat');
d12=load('eps-gam-mean-12.dat');
d08=load('eps-gam-mean-08.dat');

figure 
%plot(d08(:,1),d08(:,2),'k:')
plot(d08(:,3),d08(:,4),'k','linewidth',2)
hold on
%plot(d12(:,1),d12(:,2),'m:')
plot(d12(:,3),d12(:,4),'m','linewidth',2)
%plot(d24(:,1),d24(:,2),'b:')
plot(d24(:,3),d24(:,4),'b','linewidth',2)
%plot(d30(:,1),d30(:,2),'r:')
plot(d30(:,3),d30(:,4),'r','linewidth',2)
axis([0 1.3 0 .45])
l=legend('2.5','1.5','0.75','0.5')
set(l,'fontsize',14)
xlabel('\gamma','fontsize',18)
ylabel('\epsilon','fontsize',18,'rot',0)
%}