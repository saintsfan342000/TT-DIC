function [AramIJ_max,stg,LEp_max]=TT2_DIC_Analysis_MaxPoint();

clear;
close all;

%This function will identify the Aramis IJ indices of the point with the
%max LEp in the final stage and export this as AramIJ_max.  It will then
%go through each stage and identify the stages in which this point was not
%computed, and name these stages in "stg".
%Finally, this information will be saved to a dat file

%%%%%% PRELIMINARY DATA to ENTER %%%%%% 
    
%%%%    Define the path 
    exp=31;
%path for the files that the PERL program produces
    patha = sprintf('E:\\Martin_Experiments\\AAA_TensionTorsion\\TT2-%d\\AramisExport_MissingRemoved',exp); 
    prefix='TT2-31_FS19SS6-';
    last=350;
%path where MaxPoint.dat and Max.dat are
    pathb=sprintf('E:\\Martin_Experiments\\AAA_TensionTorsion\\TT2-%d\\TT2-%d_MatlabResults',exp,exp);
%Save path 
%prefix of the name
    %prefix = sprintf('TT2-%d_DC%d_FS%d_SS%d',TT2,tube,Facet_size,Step_size);                %MUST CHANGE don't include underscore after linear
    limts = load([pathb '\box_limts.dat']);
    Xmin=limts(1); Xmax=limts(2);Ymin=limts(3);Ymax=limts(4);
    

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%2-Calculate strain values in the scanned area.%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%Cycle through the stages
for k = last
    col_count=0;                                                  %Count the number of point columns that we'll calculate 
    %clear A LEp NEx NEy NExy gamma;
    name = sprintf('%s\\%s%d.dat',patha,prefix,k);                  %Name of the stage file to open
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
            %clear colLEp colNEx colNEy colNExy colgamma;
            col_count = col_count+1;
            %Initialize / pre-allocate length of the vectors we'll create in the for loop below
                colLEp=zeros(1,length(AramisYIndices)); colNEx=zeros(1,length(AramisYIndices)); colNEy=zeros(1,length(AramisYIndices)); 
                colNExy=zeros(1,length(AramisYIndices)); colgamma=zeros(1,length(AramisYIndices)); 
            % Calculate strain for every point in this column
            for m = 1 : length(AramisYIndices);
                F = [[A(RowsInAInTheBox(uniqeXloc(m)),9),A(RowsInAInTheBox(uniqeXloc(m)),10)] ;[A(RowsInAInTheBox(uniqeXloc(m)),11),A(RowsInAInTheBox(uniqeXloc(m)),12)]]; %transformation gradient F=RU
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
                col_I(m) = A(RowsInAInTheBox(uniqeXloc(m)),1);    
                col_J(m) = A(RowsInAInTheBox(uniqeXloc(m)),2);
            end 
        [LEp(col_count) , locLEp] = max(colLEp);  %Max LEp in the current colum (the one that passed the if test) of computation points and its location in AramisYIndices
        NEy(col_count) = colNEy(locLEp);          %Notice maxLEp and things below it are being kept as ROW VECTORS, indexing with col_count
        gamma(col_count) = colgamma(locLEp);      %So all of the strain data and coordinates of the max LEp point in each column will be saved
        aramI(col_count) = col_I(locLEp);        %After all iterations of uniqueX each of these will be a row vector with the same length (though this won't necessarily be the same length as, say, AramisXIndices since we're only scanning columns whose j-indices are unbroken
        aramJ(col_count) = col_J(locLEp);
        end                                         %End of "if column of X-index is unbroken"
    end                                             %End of  "for uniqueX=unique(AramisXIndices)'"
    %%%% Keep strain data for only those points satisfying criteria
    %%%% Criteria is to omit points whose NEy/gamma ratio doesn't fall within half of a StdDev of the mean of this ratio
    ratio = NEy ./ gamma;                           %Find the ratio of the NEy/gamma of the point that had max cum.plast.stn each scanned column
    ratioAvg = nanmean(ratio);                      %Avg ratio
    ratioSDEV = nanstd(ratio);                      %StdDev of ratio
    passed = find(ratio >= ratioAvg - 0.5 * ratioSDEV & ratio <= ratioAvg + 0.5 * ratioSDEV);    %Locations in LEp, NEx, etc vectors of the points that passed.  A row vector.
    LEp_pass=LEp(passed);
    aramI_pass=aramI(passed); 
    aramJ_pass=aramJ(passed); %Keep just those that passed.
    [LEp_max, index_max] = max(LEp_pass);                 %Find where the max LEp_pass 
    AramIJ_max=[aramI_pass(index_max) aramJ_pass(index_max)];
end    %%%  ENDS THE BIG "for k=1:last" loop.  Onto the next stage
 
%Now determine which stages this computation point shows up in
stg=[];
for k=[1:last];
    name = sprintf('%s\\%s%d.dat',patha,prefix,k); %Name of the stage file to open
    A = load(name);
    good=find(A(:,1)==AramIJ_max(1) & A(:,2)==AramIJ_max(2));   %Find will return an empty matrix there aren't any matches
    if ~isempty(good)   %We want a list of the stages that do have the point, hence the ~
        stg=[stg; k good];
    else
        sprintf('Not here %d',k)
    end;
    clear A good
end;

fclose all;

%output = fopen(sprintf('%s\\box_limts.dat',pathb),'w');    
%fprintf(output,'%% Xmin Xmax Ymin Ymax \n');
%fprintf(output,'%.15f %.15f %.15f %.15f',[Xmin Xmax Ymin Ymax]);
%fclose(output);


%%{

fid=fopen(sprintf('%s\\MaxPoint.dat',pathb),'w');
fprintf(fid,'First Row: Aramis I,J index of max point\nSecond Row and beyond: (1)Stage in which point exists (2)Row in that AramisExport_MissingRemoved file where that point is\n');
fprintf(fid,'%.0f %.0f\n',AramIJ_max');
fprintf(fid,'%.0f %.0f\n',stg');
fclose all
%}

 % If the last-stg max point doesn't appear in a sufficient number of
 % stages, consider finding the 2nd-highest point by using:
%  v1=sort(LE_pass,'descend');
%  v2=find(LE_pass==v1(2));   %Now we have the index in LE_pass of the second highest!
% secmaxXY=find(A(:,3)==xcoord_pass(index_max) & A(:,4)==ycoord_pass(index_max));
% AramIJ_2ndmax=A(secmaxXY,1:2);