function TT2_DIC_Analysis_ThirdPrincipalStretch;
clear;
close all;

%Specialty script for making epsilon-v-gamma plots
% This works a lot like Nico's original script, but it employs Dr.K's definition of epsilon and gamma
% And exports eps and gamma for (1) The mean of all points and (2)The max LEp point in each stage (not the same point throughout, necessarily)

key = xlsread('E:\Martin_Experiments\AAA_TensionTorsion\TT-Summary.xlsx');
key = key(:,[1 2 end]);

exp = [7 8 9 12 15 16 17 20 22 24 27 30];
%exp = 20;

curdir = pwd;
addpath(sprintf('%s\\MATLAB\\extras',curdir(1:2)));

for G = 1:length(exp)
    clear PATH savepath last itarray sig tau boxlims Xmin Xmax Ymin Ymax m locm loct loc
    
    last = key(key(:,1)==exp(G),end);
    alpha= key(key(:,1)==exp(G),2);

    PATH = sprintf('E:\\Martin_Experiments\\AAA_TensionTorsion\\TT2-%d\\AramisExport_MissingRemoved\\',exp(G));
    savepath = sprintf('E:\\Martin_Experiments\\AAA_TensionTorsion\\TT2-%d\\TT2-%d_MatlabResults\\',exp(G),exp(G)); 
    
% % Don't open STF, open max.dat for LV stresses
    fid = fopen(sprintf('%s\\max.dat',savepath),'r');
    td = cell2mat(textscan(fid,'%f %f %f %f %f %f %f %f','delimiter',' ','Headerlines',1));
    fclose(fid);
    sig = td(:,3);
    tau = td(:,4);
    clear td

% 
% % Open box_limits.  Box limits are in inches!
%     fid = fopen(sprintf('%s\\box_limts.dat',savepath),'r');
%     boxlims = cell2mat(textscan(fid,'%f %f %f %f','delimiter',' ','Headerlines',1));
%     fclose(fid);
%     Xmin = boxlims(1); Xmax = boxlims(2) ; Ymin = boxlims(3); Ymax = boxlims(4);
% 
% %Cycle through the stages
%     if exp(G) == 31;
%         itarray = [1:628 630:last];
%     else
%         itarray = [1:last];
%     end;
% 
%     avgL3=[];
% 
%     for k = itarray
%         col_count=0;                                                  %Count the number of point columns that we'll calculate 
%         scot_count=0;                                                     %Count the number of points in the Scott method
%         clear prefix A LEp NEx NEy NExy gamma diagU R;
%         if exp(G) == 31
%             prefix = ls(sprintf('%s*stg%d.dat',PATH,k));
%         elseif exp(G) == 30
%             prefix = ls(sprintf('%s*SS6%d.dat',PATH,k));
%         else
%             prefix = ls(sprintf('%s*_%d.dat',PATH,k));
%         end;
%         name = sprintf('%s%s',PATH,prefix);             %name of the file
%         A = load(name);
%         A(:,3:5)=A(:,3:5) / 25.4;      % Output file is in mm.   Converts to inches for all further calculations
%         A(:,6:9)=A(:,6:9)  / 25.4;
%         %Need to find the unbroken vertical colums of points in the scan zone, unbroken meaning there aren't any missing points in each column of computation points
%         %RowsInAInTheBox will give us the matlab row index in A in which there's a computation point that lies inside the box we drew.  It returns a COLUMN VECTOR of varying length.
%         RowsInAInTheBox = find( A(:,3) >= Xmin & A(:,3) <= Xmax & A(:,4) >= Ymin & A(:,4) <= Ymax);
%         %NICO: "Find the index i that are used to calculate the strain ie without index j missing"
%         AramisXIndices=A(RowsInAInTheBox,1); %AramisXIndices=A(RowsInAInTheBox,1) will give us the aramis x-axis indices of those points that lie in the box.  COLUMN VECTOR since RowsInAInTheBox is a column
%         %Loop through each unique aramis x-axis index that lies in the box
%         for uniqueX = unique(AramisXIndices)';                  %<-- transpose because this must be a row vector for the for loop to incrememnt correctly!
%             uniqeXloc = find( AramisXIndices == uniqueX );      %uniqueXloc is location in AramisXIndices, which is also the row in A, in which the particular value of uniqueX lies.  Usually a column vector.
%             AramisYIndices = A(RowsInAInTheBox(uniqeXloc),2);   %AramisYIndices is a vector of all the different Aramis y-indexes (column 2 of A) that are paired with the particular uniqueX. Same dim as uniqueXloc.
%             %This checks whether the y-indexes in AramisYIndices are contiguous (that is, if no index number is skipped), so if this "column" of particular Aramis x indices is unbroken then we'll calculate its points and "scan 
%             if (length(AramisYIndices) == max(AramisYIndices) - min(AramisYIndices) +1)
%                 clear colLEp colNEx colNEy colNExy colgamma diagU diagUxyz U Uxyz;
%                 col_count = col_count+1;
%                 %Initialize / pre-allocate length of the vectors we'll create in the for loop below
%                     colLEp=zeros(1,length(AramisYIndices)); colNEx=zeros(1,length(AramisYIndices)); colNEy=zeros(1,length(AramisYIndices)); 
%                     colNExy=zeros(1,length(AramisYIndices)); colgamma=zeros(1,length(AramisYIndices)); diagU = zeros(2,length(AramisYIndices));
%                     colLy = zeros(1,length(AramisYIndices));
%                 % Calculate strain for every point in this column
%                 for q = 1 : length(AramisYIndices);
%                     clear F U LE_calc R Uxyz N
%                     F = [[A(RowsInAInTheBox(uniqeXloc(q)),10),A(RowsInAInTheBox(uniqeXloc(q)),11)] ;[A(RowsInAInTheBox(uniqeXloc(q)),12),A(RowsInAInTheBox(uniqeXloc(q)),13)]]; %transformation gradient F=RU
%                     U = sqrtm(transpose(F)*F); 
%                     diagU(:,q) = eig(U); %principal stretch
%                     LE_calc = [[log(diagU(1,q)),0];[0,log(diagU(2,q))]];  %logarithmic strain in the principal coordinate system
%                     colLEp(q) = ( 2./3. * ( LE_calc(1,1)^2 + LE_calc(2,2)^2 + (-LE_calc(1,1)-LE_calc(2,2))^2 ))^0.5;   %Logarithmic Cumulative Plastic strain
%                     R = F * inv(U);  %Rotation tensor
%                     Uxyz = R*U*R';
%                     colLy(q) = Uxyz(2,2);                
%                     NE_calc = (U - eye([2 2]));     
%                     %Technical strain in the x',y' material coordinate system
%                     colNEx(q) = NE_calc(1,1);       %% col prefix implies column, as in NEx for each point in this column currently being calculated
%                     colNEy(q) = NE_calc(2,2);
%                     colNExy(q) = NE_calc(1,2);
%                     colgamma(q) = atan(colNExy(q)/(1+colNEx(q))) + atan(colNExy(q)/(1+colNEy(q)));
%                 end 
%             [LEp(col_count) , locLEp] = max(colLEp);  %Max LEp in the current colum (the one that passed the if test) of computation points and its location in AramisYIndices
%             NEy(col_count) = colNEy(locLEp);          %Notice maxLEp and things below it are being kept as ROW VECTORS, indexing with col_count
%             gamma(col_count) = colgamma(locLEp);      %So all of the strain data and coordinates of the max LEp point in each column will be saved
%             L1L2(:,col_count) = diagU(:,locLEp);
%             Ly(col_count) = colLy(locLEp);
%                     end                                         %End of "if column of X-index is unbroken"
%         end                                             %End of  "for uniqueX=unique(AramisXIndices)'"
%         %%%% Keep strain data for only those points satisfying criteria
%         %%%% Criteria is to omit points whose NEy/gamma ratio doesn't fall within half of a StdDev of the mean of this ratio
% 
%         ratio = NEy ./ gamma;                           %Find the ratio of the NEy/gamma each scanned column
%         ratioAvg = nanmean(ratio);                      %Avg ratio
%         ratioSDEV = nanstd(ratio);                      %StdDev of ratio
%         passed = find(ratio >= ratioAvg - 0.5 * ratioSDEV & ratio <= ratioAvg + 0.5 * ratioSDEV);    %Locations in LEp, NEx, etc vectors of the points that passed.  A row vector.
% 
%         LEp_pass=LEp(passed);
%         L1L2_pass = L1L2(:,passed);
%         L3 = 1 ./ (L1L2_pass(1,:).*L1L2_pass(2,:)) ;
%         Ly_pass = Ly(passed);
% 
% 
%         %m = [m; sig tau  nanmean(L3) sig/nanmean(L3) tau/nanmean(L3) nanmean(L3xyz) sig/nanmean(L3xyz) tau/nanmean(L3xyz)];
%         avgL3 = [avgL3;nanmean(L3) nanmean(Ly_pass)];
% 
%     end
% 
%     m = [sig(itarray) tau(itarray)  avgL3(:,1) sig(itarray)./avgL3(:,1) tau(itarray)./avgL3(:,1) avgL3(:,2) sig(itarray).*avgL3(:,2) tau(itarray).*avgL3(:,2)];

    fid = fopen(sprintf('%s\\StressPath.dat',savepath),'r');
    %fprintf(fid,'(1) Nom Sig (2) Nom Tau (3)Lambda3-Traditional (4)Tru-Sig (5)Tru-Tau (6)Lamda2-Nico (7)Tru-Sig-Nico (8)Tru-Tau-Nico\n');
    %fprintf(fid,'%f %f %f %f %f %f %f %f\n',m');
    m = cell2mat(textscan(fid,'%f %f %f %f %f %f %f %f','delimiter',' ' ,'headerlines',1));
    fclose(fid);clear fid;
    
    if exp(G) == 31;
        stgs = load(sprintf('%s\\StagesToInclude.txt',savepath));
        stgs(1) = [];
        m = m(stgs,:);
    end;

%copyfile(sprintf('%s\\StressPath.dat',savepath),sprintf('C:\\Users\\admin-local\\Desktop\\StressPathsForNico\\TT2-%d_StressPath.dat',exp(G)));

%     [~,locm] = max(sig(itarray));
%     [~,loct] = max(tau(itarray));
%     loc = max(locm,loct);
     [~,locm] = max(sig);
     [~,loct] = max(tau);
     loc = max(locm,loct);


    figure
    a1 = plot(m(:,5),m(:,4),'r:');%,'linewidth',2);
    hold on
    a2 = plot(m(:,8),m(:,7),'k:');%,'Linewidth',2);
    a3 = plot(m(:,2),m(:,1),':','Color',[0 .7 0]);%,'linewidth',2);
    plot(m(loc,5),m(loc,4),'ro','MarkerFaceColor','r');
    a4 = plot(m(loc,8),m(loc,7),'ko','MarkerFaceColor','k');
    plot(m(loc,2),m(loc,1),'o','Color',[0 .7 0],'MarkerFaceColor',[0 .7 0]);
    plot(m(last,5),m(last,4),'r^','MarkerFaceColor','r');
    a5 = plot(m(last,8),m(last,7),'k^','MarkerFaceColor','k');
    plot(m(last,2),m(last,1),'^','Color',[0 .7 0],'MarkerFaceColor',[0 .7 0]);
    l = legend([a1 a2 a3 a4 a5],{'True-\lambda_I_I_I','True-\lambda_a_x','Nominal','LL','Fail'});
    set(l,'location','southeast')
    title(sprintf('TT2-%d Stress Path.  \\alpha = %.2f',exp(G),alpha),'fontsize',14);
    xlabel('\tau','fontsize',14);
    ylabel('\sigma','fontsize',14,'rot',0);
    %axis equal
    axis([0 1.1*max(max(m(:,[2 5 8]))) 0 1.1*max(max(m(:,[1 4 7])))])
    %autoArrangeFigures;
     print(gcf,'-dpng',sprintf('%s\\StressPath',savepath));
     close
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