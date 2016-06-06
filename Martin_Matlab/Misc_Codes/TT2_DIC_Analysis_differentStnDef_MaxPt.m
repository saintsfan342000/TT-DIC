function TT2_DIC_Analysis_differentStnDef_MaxPt;
clear;
close all;
curdir = 'E:';
%This script generates eps-v-gamma for a single point in all stages
%It is built to read in the MaxPt.dat file, but any aramis IJ could be specified

key = xlsread(sprintf('%s\\Martin_Experiments\\AAA_TensionTorsion\\TT-Summary.xlsx',curdir));
key = key(:,[1 end]);

exp = [17];
for G = 1:length(exp)
    close all
    last = key(key(:,1)==exp(G),2);

    %path for the files that the PERL program produces
        if exp(G) == 17
            PATH = sprintf('%s\\Martin_Experiments\\AAA_TensionTorsion\\TT2-%d\\AramisExport_MissingRemoved_expanded\\',curdir,exp(G));
        else
            PATH = sprintf('%s\\Martin_Experiments\\AAA_TensionTorsion\\TT2-%d\\AramisExport_MissingRemoved\\',curdir,exp(G));
        end;

        savepath = sprintf('%s\\Martin_Experiments\\AAA_TensionTorsion\\TT2-%d\\TT2-%d_MatlabResults\\',curdir,exp(G),exp(G)); 

        fid = fopen(sprintf('%s\\Martin_Experiments\\AAA_TensionTorsion\\TT2-%d\\TT2-%d_MatlabResults\\MaxPoint.dat',curdir,exp(G),exp(G)));
        maxpt = cell2mat(textscan(fid,'%f %f','Delimiter',' ','Headerlines',2));
        IJ = maxpt(1,:);
        maxpt = maxpt(2:end,:);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%   1-Determination of the area to scan.    %%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%   Make the plot in which you'll draw a box %%%%%%%%%%%%%%%%%%%% 

    M=[];
    m=[];

    %Cycle through the stages

    if exp(G) == 31;
        itarray = [1:628 630:last];
    else
        itarray = [1:last];
    end;

    for k = 1:length(maxpt(:,1))
        k;
        clear prefix A LEp NEx NEy NExy gamma;
        if exp(G) == 31
            prefix = ls(sprintf('%s*stg%d.dat',PATH,maxpt(k,1)));
        elseif exp(G) == 30
            prefix = ls(sprintf('%s*SS6%d.dat',PATH,maxpt(k,1)));
        else
            prefix = ls(sprintf('%s*_%d.dat',PATH,maxpt(k,1)));
        end;
        name = sprintf('%s%s',PATH,prefix);             %name of the file
        A = load(name);
        A=A(maxpt(k,2),:);
        if A(1) ~= IJ(1) && A(2) ~= IJ(2)
            sprintf('SPORTS!')
        end
        F = [A(10) A(11);A(12) A(13)];
        U = transpose(F)*F; %stretching tensor???
        NE_calc = (U^0.5 - eye([2 2]));     
        NEx = NE_calc(1,1);       %% col prefix implies column, as in NEx for each point in this column currently being calculated
        NEy = NE_calc(2,2);
        NExy = NE_calc(1,2);
        gamma = atan(NExy/(1+NEx)) + atan(NExy/(1+NEy));
        NEy_NEW=F(2,2)-1;
        Gamma_NEW=atan(F(1,2)/F(2,2));
        M =  [M; gamma NEy abs(Gamma_NEW) NEy_NEW ];
    end

    %save(sprintf('%sEps-Gam-MaxPt.dat',savepath),'M','-ascii')

    figure
    plot(M(:,1),M(:,2),'k:')
    hold on
    plot(M(:,3),M(:,4),'k')
end
