function TT2_DIC_Analysis_profiles_all;

clear;
close all;

%Specialty script for making a profile plot of the last stage max point for all expts

key = xlsread('E:\Martin_Experiments\AAA_TensionTorsion\TT-Summary.xlsx');
key = key(:,[1 2 end]);

exp = [8 9 15 16 17 32 20 22 24 27 30];
%exp = 20;

profLEp={[]};

curdir = pwd;
addpath(sprintf('%s\\MATLAB\\extras',curdir(1:2)));

for G = 1:length(exp)
    exp(G)
    
    clear PATH savepath last itarray sig tau boxlims Xmin Xmax Ymin Ymax m locm loct loc
    
    last = key(key(:,1)==exp(G),end);
    alpha= key(key(:,1)==exp(G),2);

    if exp(G) == 17;
        PATH = sprintf('E:\\Martin_Experiments\\AAA_TensionTorsion\\TT2-%d\\AramisExport_MissingRemoved_expanded\\',exp(G));
    else
        PATH = sprintf('E:\\Martin_Experiments\\AAA_TensionTorsion\\TT2-%d\\AramisExport_MissingRemoved\\',exp(G));
    end
    
    matpath = sprintf('E:\\Martin_Experiments\\AAA_TensionTorsion\\TT2-%d\\TT2-%d_MatlabResults\\',exp(G),exp(G)); 
    
    % % Open maxpt.dat
    fid = fopen(sprintf('%s\\MaxPoint.dat',matpath),'r');
    td = cell2mat(textscan(fid,'%f %f','delimiter',' ','Headerlines',2));
    fclose(fid);
    
    k = last;
    p=1;
    
    clear A LEp NEx NEy NExy gamma;
    if exp(G) == 31
        prefix = ls(sprintf('%s*stg%d.dat',PATH,k));
    elseif exp(G) == 30
        prefix = ls(sprintf('%s*SS6%d.dat',PATH,k));
    else
        prefix = ls(sprintf('%s*_%d.dat',PATH,k));
    end;
    
    name = sprintf('%s%s',PATH,prefix,k);                  %Name of the stage file to open
    name(end)=[];
    A = load(name);
    profColum = find( A(:,1) == A(td(1,1)));  %Find all the other points that share this point's Aramis X index, and are less than 1.5 t away
   for j = 1:length(profColum) 
        %Create a cell array, each array being a profCount stage. Each array has j rows, j=# points in the Aramis X-index column
        %COLUMN 1 - Calculate Y coordinate layer :  Undef y-coord + V displacement normalized by wall thickness
        profLEp{G}(j,1) = A(profColum(j),4)+A(profColum(j),7);     
        F = [[A(profColum(j),10),A(profColum(j),11)] ;[A(profColum(j),12),A(profColum(j),13)]]; %transformation gradient F=RU
        U = transpose(F)*F;             %stretching tensor???
        diagU = eig(U);                 %principal stretch
        LE_calc = 0.5 * log(diagU);     %logarithmic strain in the principal coordinate system
        %COLUMN 2 - Log Cum Plastic strain Column 2
        profLEp{G}(j,2) = ( 2./3. * (  LE_calc(1)^2 + LE_calc(2)^2 + (-LE_calc(1)-LE_calc(2))^2) )^0.5 ;      
    end
    p=p+1;        %Incrememnt prof_count
end

c={[238 201 0]/255,[0 201 87]/355,[0 0 1],[139 58 58]/255,[0 1 0],[238 106 167]/255,[0 1 1],[255 127 36]/255,[0 0 0],[154 50 205]/255};

for i=1:length(profLEp);
    alpha(i)= key(key(:,1)==exp(i),2);
    [~,loce]=max(profLEp{i}(:,2));
    Z=profLEp{i}(loce,1);
    profLEp{i}(profLEp{i}(:,1) > (Z+2),:)=[];
    profLEp{i}(profLEp{i}(:,1) < (Z-2),:)=[];
    profLEp{i}(:,1) = profLEp{i}(:,1) - Z;
    plot(profLEp{i}(:,1),profLEp{i}(:,2),'color',c{i});
    hold on
end;

set(gcf,'color','w')

legend(strsplit(num2str(alpha)),'location','northwest')
title('e^p_e Across the Neck; Max Point')