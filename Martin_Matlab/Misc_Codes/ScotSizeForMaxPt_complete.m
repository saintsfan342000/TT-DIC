function ScotSizeForMaxPt_complete()
clear;
close all;

    
%Change mm->in
    MMtoIn = 1/25.4;
%Size of the side of the Area of calculation of Strain : Scott (grid method)
    Size_Scott = 1/16; %[in]
    
    Facet_size = 19;
    Step_size = 6;

alpha=fliplr([0.25 .5 0.5 .75 1 1 1.5 2 2.5 3 3.5 4 Inf]);
exp=fliplr([31 20 22 24 27 17 12 9 8 30 16 15 7]);

ttpath = sprintf('F:\\Martin_Experiments\\AAA_TensionTorsion\\'); 

key = xlsread('F:\Martin_Experiments\AAA_TensionTorsion\TT-Summary.xlsx');
%(1) Exp (2) Alpha (4) Rad (5) thick (end) Last stage

Mf =[]; %Max stn at failure
mf = []; %Mean stn at failure
sf = []; %Scott-size grid at failure
axsh = [];
T = [];
NewScott = [];
storealpha = [];

curdir = pwd;
addpath(sprintf('%s\\MATLAB\\extras',curdir(1:2)));

k = 1;

for G = 1:length(exp)
    
    clear PATH savepath last itarray boxlims Xmin Xmax Ymin Ymax m locm loct loc matpath fidmaxpt maxIJ maxdat maxlast 
        
    if exp(G) == 17
        PATH = sprintf('F:\\Martin_Experiments\\AAA_TensionTorsion\\TT2-%d\\AramisExport_MissingRemoved_expanded',exp(G));    
    else
        PATH = sprintf('F:\\Martin_Experiments\\AAA_TensionTorsion\\TT2-%d\\AramisExport_MissingRemoved',exp(G));
    end

    matpath = sprintf('%s\\TT2-%d\\TT2-%d_MatlabResults',ttpath,exp(G),exp(G));
    
    last = key(key(:,1)==exp(G),end);
    alpha= key(key(:,1)==exp(G),2);

    thickness = key(key(:,1)==exp(G),4);
    
     if exp(G) == 31
         prefix = ls(sprintf('%s\\*stg%d.dat',PATH,last));
     elseif exp(G) == 17;
        prefix = ls(sprintf('%s\\*%d.dat',PATH,last));
     elseif exp(G) == 30
         prefix = ls(sprintf('%s\\*SS6%d.dat',PATH,last));
     else
         prefix = ls(sprintf('%s\\*_%d.dat',PATH,last));
     end;
     
     name = sprintf('%s\\%s',PATH,prefix);             %name of the file
     
    A = load(name);
    A(:,3:9) = A(:,3:9) * MMtoIn;
    x=A(:,3);
    y=A(:,4);
    z=A(:,5);

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
    
    Step_size_th = (min_disp/thickness);
    Facet_size_th = Facet_size / Step_size * (min_disp/thickness);
    size_av = Size_Scott - thickness * Facet_size_th;

    fidmaxpt = fopen(sprintf('%s\\MaxPoint.dat',matpath),'r');
    maxdat = cell2mat(textscan(fidmaxpt,'%d %d','delimiter',' ','Headerlines',2));
    fclose(fidmaxpt);
    maxIJ = maxdat(1,:);
    maxlast = maxdat(end,:);
        
    MaxPt = A(maxlast(2),:);
    clear RowsInAInScotGrid F U RotU diagU LE_calc scotLEp NE_calc scot_x scot_y scotNEx scotNEy scotNExy scotgamma;
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
            NEx_Scott = nanmean(scotNEx); %For this one DIC point, we average all the m Scott-grid points that were near it and store this average strain data
            NEy_Scott = nanmean(scotNEy);   %scot_count(k) will be equal to length(xcoord), but not the same as col_count. Especialy true in later stages when facets drop out, so fewer columns will be contiguous
            gamma_Scott = nanmean(scotgamma);
            LEp_Scott = nanmean(scotLEp);
            OurSize_Scott = Facet_size_th * thickness + mean( [max(scot_x) - min(scot_x) , max(scot_y) - min(scot_y)] );  %Calculate our Scott-grid size based on the points we're averaging over
            NewScott = [NewScott; nanmean(OurSize_Scott),length(RowsInAInScotGrid),nanmean(NEx_Scott),nanmean(NEy_Scott),abs(nanmean(gamma_Scott)),nanmean(LEp_Scott)];
            
            clear A
            
     %Load max.dat
        fid = fopen(sprintf('%s\\max.dat',matpath),'r');
        A = cell2mat(textscan(fid,'%f %f %f %f %f %f %f %f','Delimiter',' ','Headerlines',1,'MultipleDelimsAsOne',1));
        fclose(fid);
        %Grab the stresses at last
        axsh = A(last,[3 4]);
        %Grab epeq at failure
        Mf = [Mf ; A(last,end)];
        
        clear A
        
    %Load mean.dat
        fid = fopen(sprintf('%s\\mean.dat',matpath),'r');
        A = cell2mat(textscan(fid,'%f %f %f %f %f %f %f %f','Delimiter',' ','Headerlines',1,'MultipleDelimsAsOne',1));
        fclose(fid);
        %Grab epeq at failure
        mf = [mf ; A(last,end)];
        
        clear A
        
    %Load scot.dat
        fid = fopen(sprintf('%s\\like_Scott.dat',matpath),'r');
        A = cell2mat(textscan(fid,'%f %f %f %f %f %f %f %f','Delimiter',' ','Headerlines',1,'MultipleDelimsAsOne',1));
        fclose(fid);
        %Grab epeq at failure;
        sf = [sf ; A(last,end)];
        
        clear A
        
        %Convert axsh to triax
        T = [T;(axsh(1)/2) / sqrt( (3/4)*(axsh(1)^2 + 4*axsh(2)^2) )];
        
        storealpha = [storealpha;alpha];
end;
        
%        (1)Exp (2) Alpha (3)T failure (4) Max DIC (5)Mean DIC (6)DIC-Scott
%        (7)DIC-NewScott 
    
fid = fopen(sprintf('%sCompleteFailureSummary.dat',ttpath),'w');
fprintf(fid,'(1)Exp (2) Alpha (3)T failure (4) Max DIC (5)Mean DIC (6)DIC-Scott (7)DIC-NewScott\n');
fprintf(fid,'%d %.2f %f %f %f %f %f\n',[exp' storealpha T Mf mf sf NewScott(:,end)]');
fclose(fid)

fid = fopen(sprintf('%sNewScottSummary.dat',ttpath),'w');
fprintf(fid,'(1)Exp (2) Alpha (3)Size Scott Grid (4)NbrPtsAveragedOver (5)NEx (6)NEy (7)Gamma (8)LEp\n');
fprintf(fid,'%d %.2f %f %d %f %f %f %f\n',[exp' storealpha NewScott]');
fclose(fid)

        
plot(T,sf,'o','MarkerFaceColor',[1 0 0],'MarkerEdgeColor',[1 0 0]);
hold on
plot(T,NewScott(:,end),'o','MarkerEdgeColor',[0 0 1]);
axis([0 .6 0 1.8]);
xlabel('Triaxiality: \sigma_m / \sigma_e_q','Fontsize',14);
ylabel('e^p_e_q   ','Fontsize',16,'rot',0);
title('TT2 Failure Summary');
l = legend('Averaging All Points','Maximum Point Only');
set(l,'Location','Northeast');