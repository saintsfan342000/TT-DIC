function AramisScan();

% Modified for columns delimited by commas
% Accepts # as comment lines so it ignores them

%Read in using textscan
%Write using fprintf
%Seems to be fastest method

%Not general
%Requires an Aramis export file in which each row has 13 colums, columns
%delimited by commas
%First row, first col = stage
%First row, 2nd col = time
%First row, third col = force
%First row, other columns contain meaningless commas
%2nd row: Computation point data begins 

clear all;
fclose all;

frompath='E:\zzMISC\zzARAMIS Export Files\TT2-31_FS19SS6\';     %Folder where the files you're reading from are
fromprefix='TT2-31_FS19SS6-Stage-0-';     %Prefix of the inidividual file names

topath='E:\Martin_Experiments\AAA_TensionTorsion\TT2-31';    %Where you want to create the folder of the processed files
toprefix='TT2-31_FS19SS6-';

LVfid=fopen(sprintf('%s\\TT2-35_Oct915_a0p375.dat',topath));      % FOr interpolating force from the Labview file


newfoldername=sprintf('%s\\AramisExport_MissingRemoved',topath);
mkdir(newfoldername);
web(newfoldername,'-browser')

%Number of final stage
last=350;

k=1;
for i=0:1:last;
%for i=last
    tic
    opentext=sprintf('%s\\%s%d.txt',frompath,fromprefix,i);
    fidfrom=fopen(opentext);  %Open the file you're reading from
    %Scan in the data; each row of the file reading from should contain 13 comma-separated numbers
    %textscan will read in the data and make 13 column vectors, one for
    %each column of data in the file.  Empty spots will contain NaNs.  
    %data will be a 1x13 cell, each cell being a column vector
    data=textscan(fidfrom,'%f %f %f %f %f %f %f %f %f %f %f %f',...
        'Delimiter',',','EndOfLine','\n','CommentStyle','#','MultipleDelimsAsOne',1);

        %Sadly, Aramis doesn't write a EOL character at the end of the
    %final line point in the export file, so textscan doesnt read in the very
    %last data point, so the column vector that textscan creates is missing
    %one element.  Therefore, the 13th column vectors in the cell array is
    %one element shorted than the other 12, so the cell is non-rectangular
    %and can't be converted directly to matrix.
    %Indeed, direct from the Mathworks help page on textscan: 
    %"When processing numeric data, textscan also ignores trailing white space."
    %So I truncate each of the
    %other arrays in the cell.  This has the end  effect of removing the
    %last line of data (which is usually missing data anyways)
    L=length(data{12});
    for z=1:11;
        if length(data{z})~=L;
            data{z}(L)=[];
        end;
    end;
        data=cell2mat(data);    %Convert the cell to matrix
    fclose(fidfrom);        %Close the file we read from
    data(any(isnan(data),2),:)=[]; 
    %Removes rows that contained NaN due to empty spot in Aramis export file 
    %If mx is a p-by-q matrix 
    %any(mx,1) returns a 1-by-q vector with 1 in the locations where that column contained a non-zero.
    %any(mx,2) returns a p-by-1 vector with 1 in locations corresponding to rows in MX that contained a non-zero.  
    %any(mx,3) returns a matrix of same size as mx, with ones in the locations where mx  non-zero
    stage_time_force(k,:)=data(1,1:3);  %Pull out stage/time/force data from the first row
    data(1,:)=[];  %Now delete the first row
    savestring=sprintf('%s\\%s%d.dat',newfoldername,toprefix,i);    %File name with stage number
    fidto=fopen(savestring,'w');   %a+ simply means we append data to that already in the file, which in this case is none.
    fprintf(fidto,'%d,%d,%.6f,%.6f,%.6f,%.6f,%.6f,%.6f,%.6f,%.6f,%.6f,%.6f\r\n',data');
    fclose(fidto);
    clear data savestring
    k=k+1;
    toc
    end;
% 
%LV=cell2mat(textscan(LVfid,'%f %f %f %f %f','Delimiter','\t','Headerlines',1));
%fclose(LVfid);
%V(:,1)=LV(:,1)/(2*pi()*0.0381*0.8937^2)/1000; %Convert labview force to ksi stress, and torque to ksi stress
%LV(:,3)=LV(:,3)/(2*pi()*0.0381*0.8937)/1000; %Convert labview force to ksi stress, and torque to ksi stress
%LVinterp=interp1(LV(:,5),LV(:,[3,1]),stage_time_force(:,2)); %Col1 = force, Col2 = torque


%stage_time_force=[stage_time_force(:,1:2) LVinterp];

%savestring2=sprintf('%s\\time_force.dat',newfoldername);    %Name of time_force
%fidforce=fopen(savestring2,'a+');   %Here again a+ ins't imortant
%printf(fidforce,'%d %.4f %.6f %.6f\r\n',stage_time_force');
%fclose(fidforce);
