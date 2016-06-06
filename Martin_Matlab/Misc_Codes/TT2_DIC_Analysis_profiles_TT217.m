function TT2_DIC_Analysis_profiles();
clear;
close all;
curdir=pwd;

% This script generates equivalent-plastic strain profiles 
% Necessary inputs are listed below
% Of highest importance is entering the Aramis i index of the point in the last stage with max epeq
% The program wil read in the expt's time_force.dat file, and automaticaly generate profiles
%   at 10 equally-spaced stages from load maximum onwards.  
% If you prefer, though, you can manually enter the stages
% The data is saved into an excel file, with each profile stored in a worksheet
    

%%%%%% PRELIMINARY DATA to ENTER %%%%%% 
    TT2=17;             %Expt number
    last = 1070;         %Last stage in which specimen is not fractures
%     profMaxlocA = 84;   %Determined to be the i index of the max stn point in the last stage 84 is best so far
%     prefix = 'AramisExport_MissingRemoved_expanded\TT2-17_DC15_FS19_SS6_';                %MUST CHANGE don't include underscore after linear
     profMaxlocA = 196;   %Determined to be the i index of the max stn point in the last stage 84 is best so far
     prefix = 'AramisExport_MissingRemoved_expanded\TT2-17_expanded_FS19SS6_';                %MUST CHANGE don't include underscore after linear
    
savestuff = 1;
    
PATH = sprintf('%s\\Martin_Experiments\\AAA_TensionTorsion\\TT2-%d\\',curdir(1:2),TT2); %MUST CHANGE

addpath(sprintf('%s\\Matlab\\extras',curdir(1:2)));     %Adds export_fig



%Stage Time Force data...Create the 10 points we'll highlight in Dr K's classic figure
STF = load(sprintf('%s\\AramisExport_MissingRemoved\\time_force.dat',PATH));                    %Read in and load stage time force data.
STF((last+2):length(STF),:)=[];                                  %In case I exported  more un-cleaned Aramis files than the last stage PRIOR to failure
                                                                 %Remember - last is the last STAGE NUMBER.  So length(STF) = last+1.  So we want to cut last+2 and beyond.
STF(1,:)=[];                                                     %Get rid of the first line (stage zero) so that indexing is consistent
%STF(:,3)=STF(:,3)*cal/stsfactor;                                 % Convert voltage to stress
[~,locf] = max(STF(:,3));                                     %Need the index of the max force
%[~,locs] = max(STF(:,4));
%locf = max([locs locf]);
prof_num = 10;                                                   %Number of profile points we'll highlight in our figures
prof_inc = round((last-locf) / (prof_num-1));                    %Determine the equally-spaced stage increment between each of these points
profStages = [locf:prof_inc:locf+(prof_num-2)*prof_inc,last]';   %profStages is the stage numbers of our 10 pts.  Column vector for printf purposes 

profStages = [ 200   320  444   468   492   516   540   564   588   609]';
prof_num = length(profStages);

%Initialize a few things before looping and calculating every stage

prof_count=1;
profLEp{prof_num}=[];

%Cycle through the stages
p=1;
for k = [profStages]';
    col_count=0;                                                  %Count the number of point columns that we'll calculate 
    scot_count=0;                                                     %Count the number of points in the Scott method
    clear A LEp NEx NEy NExy gamma;
    name = sprintf('%s%s%d.dat',PATH,prefix,k);                  %Name of the stage file to open
    A = load(name);
    profColum = find( A(:,1) == A(profMaxlocA,1));  %Find all the other points that share this point's Aramis X index, and are less than 1.5 t away
   for j = 1:length(profColum) 
        %Create a cell array, each array being a profCount stage. Each array has j rows, j=# points in the Aramis X-index column
        %COLUMN 1 - Calculate Y coordinate layer :  Undef y-coord + V displacement normalized by wall thickness
        profLEp{prof_count}(j,1) = (A(profColum(j),4)+A(profColum(j),7))/25.4/.0382;     
        F = [[A(profColum(j),10),A(profColum(j),11)] ;[A(profColum(j),12),A(profColum(j),13)]]; %transformation gradient F=RU
        U = transpose(F)*F;             %stretching tensor???
        diagU = eig(U);                 %principal stretch
        LE_calc = 0.5 * log(diagU);     %logarithmic strain in the principal coordinate system
        %COLUMN 2 - Log Cum Plastic strain Column 2
        profLEp{prof_count}(j,2) = ( 2./3. * (  LE_calc(1)^2 + LE_calc(2)^2 + (-LE_calc(1)-LE_calc(2))^2) )^0.5 ;      
    end
    prof_count=prof_count+1;        %Incrememnt prof_count
end;
[~,loce]=max(profLEp{end}(:,2));
Z=profLEp{end}(loce,1);
for i=1:length(profLEp);
     profLEp{i}(profLEp{i}(:,1) > (Z+2),:)=[];
     profLEp{i}(profLEp{i}(:,1) < (Z-2),:)=[];
end;


%%{
%Y coord vs e^p profiles
c={[238 201 0]/255,[0 201 87]/355,[0 0 1],[139 58 58]/255,[0 1 0],[238 106 167]/255,[0 1 1],[255 127 36]/255,[0 0 0],[154 50 205]/255};

figure
% subplot(2,1,1)
for i = 1 : length(profLEp)
    hold on
    plot(profLEp{i}(:,1),profLEp{i}(:,2),'Color',c{i},'Linewidth',2);
end
%axis([-4 2 0 .8])
title(sprintf('TT2-%d',TT2))
xlabel('y/t','Fontsize',14)
ylabel('e^p','Fontsize',14,'Rotation',0)
set (gca,'Fontsize',14)
set(gcf, 'color', [1 1 1] );
l = legend(strsplit(num2str((profStages'))));
set(l,'fontsize',8,'location','eastoutside')
hold off

% subplot(2,1,2)
% plot(STF(:,3));
% hold on
% plot(profStages,STF(profStages,3),'ro')
% xlabel('Stage')
% ylabel('Load')
% set(gcf,'PaperPositionMode','auto')
% set(gcf, 'Position', [0 0 1.2*500 1.2*900])

if savestuff == 1

     for i=1:length(profLEp);
         xlswrite(sprintf('%sTT2-%d_MatlabResults\\StrainProfiles_ColIndex%dNew.xlsx',PATH,TT2,profMaxlocA),profLEp{i},sprintf('Stage%d',profStages(i)));
     end;

     export_fig(sprintf('%sTT2-%d_MatlabResults\\StrainProfilesNew.png',PATH,TT2),'-r200')
 
end

%print(gcf,'-dpng',sprintf('%s/ep_profile',folder))
%export_fig(sprintf('%s\\ep_profile.png',folder));
%{
for i=1:prof_num                                        %PROFILE DATA
  clear out;
  name = sprintf('%s/export_Y_LEp_time_%6.3f_stress_%6.3f.dat',folder2,STF(profStages(i),2),STF(profStages(i),3));
  output = fopen(name,'w');
  fprintf(output,'%s Y/th [] True Cumulative plastic strain [] \n','%');
  out = profLEp{i}(:,:);
  fprintf(output,'%f %f',out');
  save(name,'out','-ASCII');
end
 
% Key of the different columns in the files we're reading in
%A(:,1) index i [] corresponding to the x axis. Note Aramis i,j index seems to always start in bottom left corner of facet field, though this is not he origin of the coordinate system
%A(:,2) index j [] corresponding to the y axis
%A(:,3:5) undeformed coordinates X Y Z [mm] 
%A(:,6:9) displacement U V W ur=w [mm] 
%A(:,10:13) Surface gradient transformation F(0,0) F(0,1) F(1,0) F(1,1) []

fclose all;
    %}