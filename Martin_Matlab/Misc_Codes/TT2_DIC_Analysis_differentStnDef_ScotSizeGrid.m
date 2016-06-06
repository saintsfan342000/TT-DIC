function TT2_DIC_Analysis_differentStnDef_ScotSizeGrid;
clear;
close all;

c={[0 0 0],[255 130 171]/255,[156 102 31]/255,[0 0 1],[1 0 0],[0 201 87]/355,[139 87 66]/255,[255 165 79]/255,[178 58 238]/255,[0 .95 0],[0 220 255]/255,[0 0 0],[255 130 171]/255,[156 102 31]/255,[0 0 1],[1 0 0]};

%Specialty script for making epsilon-v-gamma plots
% This works a lot like Nico's original script, but it employs Dr.K's definition of epsilon and gamma
% And exports eps and gamma for (1) The mean of all points and (2)The max LEp point in each stage (not the same point throughout, necessarily)

key = xlsread('E:\Martin_Experiments\AAA_TensionTorsion\TT-Summary.xlsx');

exp = [31];
%exp = 20;

Size_Scott = 1/16;
Facet_size = 19;
Step_size = 6;
MMtoIn = 1 / 25.4;

curdir = pwd;
addpath(sprintf('%s\\MATLAB\\extras',curdir(1:2)));
figure
for G = 1:length(exp)
    clear PATH savepath last itarray sig tau boxlims Xmin Xmax Ymin Ymax m locm loct loc strains min_disp
    
    last = key(key(:,1)==exp(G),end);
    alpha= key(key(:,1)==exp(G),2);
    thickness = key(key(:,1)==exp(G),5);
    
    if exp(G) ~= 17
        PATH = sprintf('E:\\Martin_Experiments\\AAA_TensionTorsion\\TT2-%d\\AramisExport_MissingRemoved\\',exp(G));
    else
        PATH = sprintf('E:\\Martin_Experiments\\AAA_TensionTorsion\\TT2-%d\\AramisExport_MissingRemoved_expanded\\',exp(G));
    end;
    
    savepath = sprintf('E:\\Martin_Experiments\\AAA_TensionTorsion\\TT2-%d\\TT2-%d_MatlabResults\\',exp(G),exp(G));
    
    % % Don't open STF, open max.dat for LV stresses
    fid = fopen(sprintf('%s\\max.dat',savepath),'r');
    td = cell2mat(textscan(fid,'%f %f %f %f %f %f %f %f','delimiter',' ','Headerlines',1));
    fclose(fid);
    sig = td(:,3);
    tau = td(:,4);
    clear td
        
    % Open box_limits.  Box limits are in inches!
    fid = fopen(sprintf('%s\\box_limts.dat',savepath),'r');
    boxlims = cell2mat(textscan(fid,'%f %f %f %f','delimiter',' ','Headerlines',1));
    fclose(fid);
    Xmin = boxlims(1); Xmax = boxlims(2) ; Ymin = boxlims(3); Ymax = boxlims(4);
    
    %Need to open up last stage and calc min_disp
                    if exp(G) == 31;
                        prefix = ls(sprintf('%s*stg%d.dat',PATH,last));
                    elseif exp(G) == 30;
                        prefix = ls(sprintf('%s*SS6%d.dat',PATH,last));
                    else;
                        prefix = ls(sprintf('%s*_%d.dat',PATH,last));
                    end;
                    name = sprintf('%s%s',PATH,prefix);             
                    A = load(name);          
                    x=A(:,3) .* MMtoIn;     
                    y=A(:,4) .* MMtoIn;
                    z=A(:,5) .* MMtoIn;

                    for i = 1:length(x); 
                        F = [[A(i,10),A(i,11)] ;[A(i,12),A(i,13)]];
                        U = transpose(F)*F;
                        [~,diagU] = eig(U); 
                        LE_calc = 0.5 * [[log(diagU(1,1)),0];[0,log(diagU(2,2))]];
                        LEp(i) = ( 2./3. * ( LE_calc(1,1)^2 + LE_calc(2,2)^2 + (-LE_calc(1,1)-LE_calc(2,2))^2 ))^0.5 ;
                    end

                    min_disp = 1000;
                    for k = [1:round(length(x)/30)];
                        for q= [1:round(length(x)/30)];
                            dist= ( (x(k)-x(q))^2 + (y(k)-y(q) )^2 + ( z(k)-z(q) )^2 )^0.5 ; 
                            if dist==0
                                dist = 11000;
                            end
                            min_disp=min(min_disp,dist);
                        end
                    end
    
    Facet_size_th = Facet_size / Step_size * (min_disp/thickness);
    size_av = Size_Scott - thickness * Facet_size_th;
    
    %Initialize strains
    strains = [];
    
    %Cycle through the stages
    if exp(G) == 31;
        itarray = [1:628 630:last];
    else
        itarray = [1:last];
    end;
    
    for k = itarray
        if any(k==[0:25:1500])
            [exp(G) k last]
        end;
        col_count=0;                                                  %Count the number of point columns that we'll calculate
        scot_count=0;                                                     %Count the number of points in the Scott method
        clear prefix A LEp NEx NEy NExy gamma diagU R xcoord y coord RowsInAInTheBox AramisXIndices;
        if exp(G) == 31
            prefix = ls(sprintf('%s*stg%d.dat',PATH,k));
        elseif exp(G) == 30
            prefix = ls(sprintf('%s*SS6%d.dat',PATH,k));
        else
            prefix = ls(sprintf('%s*_%d.dat',PATH,k));
        end;
        name = sprintf('%s%s',PATH,prefix);             %name of the file
        A = load(name);
        A(:,3:5)=A(:,3:5) / 25.4;      % Output file is in mm.   Converts to inches for all further calculations
        A(:,6:9)=A(:,6:9)  / 25.4;
        %Need to find the unbroken vertical colums of points in the scan zone, unbroken meaning there aren't any missing points in each column of computation points
        %RowsInAInTheBox will give us the matlab row index in A in which there's a computation point that lies inside the box we drew.  It returns a COLUMN VECTOR of varying length.
        RowsInAInTheBox = find( A(:,3) >= Xmin & A(:,3) <= Xmax & A(:,4) >= Ymin & A(:,4) <= Ymax);
        %NICO: "Find the index i that are used to calculate the strain ie without index j missing"
        AramisXIndices=A(RowsInAInTheBox,1); %AramisXIndices=A(RowsInAInTheBox,1) will give us the aramis x-axis indices of those points that lie in the box.  COLUMN VECTOR since RowsInAInTheBox is a column
        %Loop through each unique aramis x-axis index that lies in the box
        clear uniqueX
        for uniqueX = unique(AramisXIndices)';                  %<-- transpose because this must be a row vector for the for loop to incrememnt correctly!
            clear uniqeXloc AramisYIndices 
            uniqeXloc = find( AramisXIndices == uniqueX );      %uniqueXloc is location in AramisXIndices, which is also the row in A, in which the particular value of uniqueX lies.  Usually a column vector.
            AramisYIndices = A(RowsInAInTheBox(uniqeXloc),2);   %AramisYIndices is a vector of all the different Aramis y-indexes (column 2 of A) that are paired with the particular uniqueX. Same dim as uniqueXloc.
            %This checks whether the y-indexes in AramisYIndices are contiguous (that is, if no index number is skipped), so if this "column" of particular Aramis x indices is unbroken then we'll calculate its points and "scan
            if (length(AramisYIndices) == max(AramisYIndices) - min(AramisYIndices) +1)
                clear colLEp colNEx colNEy colNExy colgamma diagU col_x col_y colgamma colNExy colNEy colNEx NE_calc;
                col_count = col_count+1;
                %Initialize / pre-allocate length of the vectors we'll create in the for loop below
                colLEp=zeros(1,length(AramisYIndices)); colNEx=zeros(1,length(AramisYIndices)); colNEy=zeros(1,length(AramisYIndices));
                colNExy=zeros(1,length(AramisYIndices)); colgamma=zeros(1,length(AramisYIndices)); diagU = zeros(2,length(AramisYIndices));
                col_x = zeros(1,length(AramisYIndices));col_y = zeros(1,length(AramisYIndices));
                % Calculate strain for every point in this column
                for q = 1 : length(AramisYIndices);
                    clear F U LE_calc R Uxyz N 
                    F = [[A(RowsInAInTheBox(uniqeXloc(q)),10),A(RowsInAInTheBox(uniqeXloc(q)),11)] ;[A(RowsInAInTheBox(uniqeXloc(q)),12),A(RowsInAInTheBox(uniqeXloc(q)),13)]]; %transformation gradient F=RU
                    U = sqrtm(transpose(F)*F);
                    diagU(:,q) = eig(U); %principal stretch
                    LE_calc = [[log(diagU(1,q)),0];[0,log(diagU(2,q))]];  %logarithmic strain in the principal coordinate system
                    colLEp(q) = ( 2./3. * ( LE_calc(1,1)^2 + LE_calc(2,2)^2 + (-LE_calc(1,1)-LE_calc(2,2))^2 ))^0.5;   %Logarithmic Cumulative Plastic strain
                    NE_calc = (U - eye([2 2]));
                    %Technical strain in the x',y' material coordinate system
                    colNEx(q) = NE_calc(1,1);       %% col prefix implies column, as in NEx for each point in this column currently being calculated
                    colNEy(q) = NE_calc(2,2);
                    colNExy(q) = NE_calc(1,2);
                    colgamma(q) = atan(colNExy(q)/(1+colNEx(q))) + atan(colNExy(q)/(1+colNEy(q)));
                    col_x(q) = A(RowsInAInTheBox(uniqeXloc(q)),3);
                    col_y(q) = A(RowsInAInTheBox(uniqeXloc(q)),4);
                end
                [LEp(col_count) , locLEp] = max(colLEp);
                xcoord(col_count) = col_x(locLEp);        
                ycoord(col_count) = col_y(locLEp);
            end                                         
        end                                             
        
        for i = 1:length(xcoord); 
            clear RowsInAInScotGrid F U RotU diagU LE_calc scotLEp NE_calc scot_x scot_y scotNEx scotNEy scotNExy scotgamma colNEyS colgammaS;
            scot_count = scot_count + 1;     
            RowsInAInScotGrid = find( A(:,3)<=(xcoord(i)+size_av/2) & A(:,3)>=(xcoord(i)-size_av/2) & A(:,4)<=(ycoord(i)+size_av/2) &  A(:,4)>=(ycoord(i)-size_av/2));
            for m = 1:length(RowsInAInScotGrid);
                clear F U LE_calc
                F = [[A(RowsInAInScotGrid(m),10),A(RowsInAInScotGrid(m),11)] ;[A(RowsInAInScotGrid(m),12),A(RowsInAInScotGrid(m),13)]];     %transformation gradient F=RU
                colNEyS(m) = F(2,2) - 1;
                colgammaS(m) = atan(F(1,2)/F(2,2));
                U = transpose(F)*F;                         %stretching tensor???
                [~,diagU] = eig(U);                      %principal stretch
                LE_calc = 0.5 * [[log(diagU(1,1)),0];[0,log(diagU(2,2))]];  %logarithmic strain in the principal coordinate system
                scotLEp(m) = ( 2./3. * ( LE_calc(1,1)^2 + LE_calc(2,2)^2 + (-LE_calc(1,1)-LE_calc(2,2))^2 ))^0.5 ;   %Logarithmic Cumulative Plastic strain
                %R = F * U^(-0.5);                          %Rotation tensor
                NE_calc = (U^0.5 - eye([2 2]));             %Technical strain calculation in the stretching coordinate system
                scotNEx(m) = NE_calc(1,1);
                scotNEy(m) = NE_calc(2,2);
                scotNExy(m) = NE_calc(1,2);
                scotgamma(m) = atan(scotNExy(m)/(1+scotNEx(m))) + atan(scotNExy(m)/(1+scotNEy(m)));
             end
            NEySnew(scot_count) = nanmean(colNEyS); %For this one DIC point, we average all the m Scott-grid points that were near it and store this average strain data
            gammaSnew(scot_count) = nanmean(colgammaS);   %scot_count(k) will be equal to length(xcoord), but not the same as col_count. Especialy true in later stages when facets drop out, so fewer columns will be contiguous
            NEySold(scot_count) = nanmean(scotNEy);
            gammaSold(scot_count) = nanmean(scotgamma);
        end
        
        strains = [strains; abs(nanmean(gammaSold)) nanmean(NEySold) abs(nanmean(gammaSnew)) nanmean(NEySnew)];
        clear gammaSold NEySold gammaSnew NEySnew
        
    end
    
    save(sprintf('%s\\eps-gam_GridSize.dat',savepath),'strains','-ascii');
 
    [~,locm] = max(sig);
    [~,loct] = max(tau);
    loc = max(locm,loct);
    
    if exp(G) == 31;
        stgs = load(sprintf('%s\\StagesToInclude.txt',savepath));
        stgs(1) = [];
        strains = strains(stgs,:);
        sig = sig(stages);
        tau = tau(stages);
        [~,locm] = max(sig);
        [~,loct] = max(tau);
        loc = max(locm,loct);
    end;
       
    plot(strains(:,3),strains(:,4),'color',c{G});
    plot(strains(loc,3),strains(loc,4),'rs','MarkerFaceColor','r');
    hold on
    
    %xlswrite('E:\Martin_Experiments\AAA_TensionTorsion\Paper and Conference\epsilon_v_gamma\Eps_v_gamma_LinearPaths_ScotSize.xlsx',{'ColA = Gamma, ColB = Eps.  FirstRow = Limit Load'},sprintf('TT2-%d, alpha = %.2f',exp(G),alpha),'A1:A1')
    %xlswrite('E:\Martin_Experiments\AAA_TensionTorsion\Paper and Conference\epsilon_v_gamma\Eps_v_gamma_LinearPaths_ScotSize.xlsx',strains(loc,[3 4]),sprintf('TT2-%d, alpha = %.2f',exp(G),alpha),'A2')
    %xlswrite('E:\Martin_Experiments\AAA_TensionTorsion\Paper and Conference\epsilon_v_gamma\Eps_v_gamma_LinearPaths_ScotSize.xlsx',strains(:,[3 4]),sprintf('TT2-%d, alpha = %.2f',exp(G),alpha),'A3')
end

