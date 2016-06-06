function EndSxnAnalysis();
clear all
close all
clc
tic

savefile=1;

frompath='E:\zzMISC\zzARAMIS Export Files\TT2-35_DC16_BC';     %Folder where the files you're reading from are
fromprefix='TT2-35_point';     %Prefix of the inidividual file names
filename='E:\Martin_Experiments\AAA_TensionTorsion\TT2-35\TT2-35_BCInfo.xlsx'; %Path and name of the excel file we'll save to

num = 19;   %Number of points(minus 1) since first is ID zero
            %Numbers 0-9 must be above text sxn
            %Numbers 10-19 must be below test sxn

for i=1:1:num+1;    
    opentext=sprintf('%s\\%s%d.txt',frompath,fromprefix,i-1);
    fidfrom=fopen(opentext);  %Open the file you're reading from
    data=textscan(fidfrom,'%f %f %f %f %f %f %f',...
        'Delimiter',',','EndOfLine','\n','CommentStyle','#','MultipleDelimsAsOne',1);
    d{i}=cell2mat(data);
    clear data
    index_xy(i,:)=d{i}(1,1:2);   %copy aramis indices
    undef(i,:)=d{i}(1,3:5);     %Copy undef coordinates
    d{i}(1,:)=[];               %Delete top row (only contains aramix index_x,y, undef XYZ)
    if i==1;
        STF(:,1:2)=d{i}(:,1:2);       %Copy STF.  Only do it once since it's same for all points
%%%% MUST BE CAREFUL HERE.  IF ALPHA >= 3.25 then we're outputting torque not force!  And coversion is 2000
        STF(:,3)=d{i}(:,3)*2000;    %Copy force and convert voltage to force
    end;
    d{i}(:,1:3)=[];             %Delete STF info. d{i} now should have 7 columns
    fclose(fidfrom);        %Close the file we read from
end;

%   d{i} columns
%   (1)DispX (2)DispY (3)DispZ (4)RotationZ

%Plotting
if savefile~=1;
    %Generate plot markers
    marker_type = ['.';'.';'.';'.'];
    marker_color = ['r';'g';'b';'m';'k'];
    marker_style(num,:) = '  ';
    m0=1; n0=1; i=1;  %Initialize
    while i <= num+1 %prof_num is the number of profile stages we're highlighting (10)
        marker_style(i,:) = [marker_color(m0) marker_type(n0)]; 
        i=i+1;
        if m0<length(marker_color)
            m0=m0+1;
        else
            m0=1;
        end;
        if n0<length(marker_type)
            n0=n0+1;
        else
            n0=1;
        end;
    end

    %Show initial z coord is ~same for each point 
    figure
    for i=1:10;
            plot(i,undef(i,3),marker_style(i,:),'Linewidth',1);
            hold on
    end;
    k=1;
    for i=11:20;
        plot(k,undef(i,3),marker_style(i,:),'Linewidth',1);
        hold on
        k=k+1;
    end;
    hold off

    %ok, sometimes one of the points may disappear for some stages
    %This finds which points have stages missing
    elim=[];
    for i=1:num;
        if length(d{i})~=length(STF(:,1));
            elim=[elim, i];
        end;
    end;

    %Show that none of the rotations of top and bottom deviate from the rest
    for k=[1 2 3 4];
        figure
        inc=1:1:num;
        inc(elim)=[];   %Notice now we incrememnt only thru the points that aren't missing stages
        for i=inc;
            plot(STF(:,2),d{i}(:,k),marker_style(i,:),'Linewidth',1);
            hold on
            xlabel('Time (sec)')
            if k==1;
                ylabel('X Disp (Horizontal) (deg)')
            elseif k==2; %%These deviate b/c the Y axis points "into the image" - not radial coordinate
                ylabel('Y Disp - LESS MEANINGFUL')
            elseif k==3;
                ylabel('Z Disp - Tube Axial')
            elseif k==4;
                ylabel('Rotation from Initial (deg)')
            end
        end;
    end;
end;

% Get initial Z Displacement
ZInitial(1,1)=mean(undef(1:10,3));  %Average inital Z coord of upper row
ZInitial(2,1)=mean(undef(11:20,3));  %Average initial Z coord of lower row


%d{i} columns:   (1)DispX (2)DispY (3)DispZ (4)RotationZ
clear ForMeans
%%% Z DISP
    %Upper Points
        k=1;
        for i=1:10;
            if length(d{i})==length(STF(:,1));
                ForMeans(k,:)=d{i}(:,3)';
                k=k+1;
            end;
        end;
        ZdispMean(1,:)=mean(ForMeans);  %Mean Z Displacement of the TOP row of points
        clear ForMeans
    %Lower Points
        k=1;
        for i=11:20;
            if length(d{i})==length(STF(:,1));
                ForMeans(k,:)=d{i}(:,3)';
                k=k+1;
            end;
        end;
        ZdispMean(2,:)=mean(ForMeans);  %Mean Z Displacement of the BOTTOM row of points
        ZDist=abs((ZInitial(1,1)+ZdispMean(1,:))-(ZInitial(2,1)+ZdispMean(2,:)));

%d{i} columns:   (1)DispX (2)DispY (3)DispZ (4)RotationZ
clear ForMeans
%%% Rotation
    %Upper Points
        k=1;
        for i=1:10;
             if length(d{i})==length(STF(:,1));
                ForMeans(k,:)=d{i}(:,4)';
                k=k+1;
            end;
        end;
        RotMean(1,:)=mean(ForMeans);  %Mean Rotation of the TOP row of points
        clear ForMeans
    %Lower Points
        k=1;
        for i=11:20;
            if length(d{i})==length(STF(:,1));
                ForMeans(k,:)=d{i}(:,4)';
                k=k+1;
            end;
        end;
        RotMean(2,:)=mean(ForMeans);  %Mean Rotation of the BOTTOM row of points
        RelativeRot=abs(RotMean(2,:)-RotMean(1,:));
    plot(STF(:,2),RotMean(2,:))
    hold on
    plot(STF(:,2),RotMean(1,:))

% 
if savefile==1;
xlswrite(filename,{'Stage'},'Sheet1','A1:A1');
xlswrite(filename,{'Time (sec)'},'Sheet1','B1:B1');
xlswrite(filename,{'Aramis Force (lb)'},'Sheet1','C1:C1');
xlswrite(filename,{'Axial Disp on Top(mm)'},'Sheet1','D1:D1');
xlswrite(filename,{'Axial Disp on Bottom(mm)'},'Sheet1','E1:E1');
xlswrite(filename,{'Axial Separation (mm)'},'Sheet1','F1:F1');
xlswrite(filename,{'Rotation on Top (deg)'},'Sheet1','G1:G1');
xlswrite(filename,{'Rotation on Bottom (deg)'},'Sheet1','H1:H1');
xlswrite(filename,{'Relative Rotation (deg)'},'Sheet1','I1:I1');
xlswrite(filename,[STF ZdispMean' ZDist' RotMean' RelativeRot'],'Sheet1','A2');
end;

%{
%d{i} columns:   (1)DispX (2)DispY (3)DispZ (4)RotationZ
clear ForMeans
%%% X DISP 
    %Upper Points
        k=1;
        for i=1:10; %Now the points that don't give data for every stage are excluded
            if length(d{i})==length(STF(:,1));
                ForMeans(k,:)=d{i}(:,1)';
                k=k+1;
            end;
        end;
        XdispMean(1,:)=mean(ForMeans);  %Mean X Displacement of the TOP row of points
        clear ForMeans
    %Lower Points
        k=1;
        for i=11:20;
            if length(d{i})==length(STF(:,1));
                ForMeans(k,:)=d{i}(:,1)';
                k=k+1;
            end;
        end;
        XdispMean(2,:)=mean(ForMeans);  %Mean X Displacement of the BOTTOM row of points
%}