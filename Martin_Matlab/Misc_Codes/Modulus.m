    function mods = Modulus();

ex = [8972.2       3723.8];
ab = [8917 3893]    
    
(ex - ab)./ex * 100 = [0.61554    -4.544];

%%% This function calculates the modulus response of each TT experiment
%%% It is set up to open up BCInfo xls sheets given only exp ID #
%%% TT-Summary.xlsx is very important here and tells script which stage is
%%% linear upperbound.

clear;
close all;

curdir = pwd;

curdir = curdir(1:2);

key = xlsread([curdir '\Martin_Experiments\AAA_TensionTorsion\TT-Summary.xlsx']);
key = key(:,[1 2 4 5 6 7 end-1 end]);
    xpts = key(:,1);
    alpha = key(:,2);
    Rm = key(:,3);
    tav = key(:,4);
    sgyld = key(:,5);
    tauyld = key(:,6);
    stgyld = key(:,7);
    last = key(:,end);

exp = [8 9 15 16 17 32 20 24 30 32];

mods = [];
ylds=[];

for G = 1:length(exp);
    
    R = Rm(xpts == exp(G));
    t = tav(xpts == exp(G));
    
    lub = stgyld(xpts ==exp(G));
    
    exppath = sprintf('%s\\Martin_Experiments\\AAA_TensionTorsion\\TT2-%d\\',curdir,exp(G)); 
    matpath = sprintf('%s\\Martin_Experiments\\AAA_TensionTorsion\\TT2-%d\\TT2-%d_MatlabResults\\',curdir,exp(G),exp(G)); 
    
    del=xlsread(sprintf('%sTT2-%d_BCInfo.xlsx',exppath,exp(G)),'F:F');
    %Note that this is not the displacement but the "axial separation".
    % Convert to disp/Lg by subtracting off initial value then diving thru by initial value
    initD = del(1);
    del=(del-del(1))/initD;
    del(1)=[];
    
    phi=xlsread(sprintf('%sTT2-%d_BCInfo.xlsx',exppath,exp(G)),'I:I') * pi()/180 ;
    phi(1)=[];
    
    fid = fopen( sprintf( '%smax.dat', matpath) , 'r');
    d = cell2mat(textscan(fid,'%f %f %f %f %f %f %f %f','headerlines',1,'delimiter',' '));
    fclose(fid);
       
    E = polyfit(del(1:lub),d(1:lub,3),1);
    E = E(1);

    
    mu = polyfit(phi(1:lub)*(R+t/2)/(initD/25.4),d(1:lub,4),1);
    mu = mu(1);
    
    ylds = [ylds; exp(G) alpha(xpts==exp(G)) d(lub,3) sgyld(xpts==exp(G)) d(lub,4) tauyld(xpts==exp(G)) ];
    
    mods = [mods ;exp(G) alpha(xpts==exp(G)) E mu];
    
end;    
    
mods = sortrows(mods,2);

mean(mods(:,3:4))

plot(mods(:,2),mods(:,3),'o')
figure()
plot(mods(:,2),mods(:,4),'o')

% Epos = s12*(2*a*e12 - e11 + e22 - abs(2*a*e12 - e11 + e22))/(2*e12*(a*e12 + e22));
% Eneg = s12*(2*a*e12 - e11 + e22 + abs(2*a*e12 - e11 + e22))/(2*e12*(a*e12 + e22));
% s22pos = s12*(-e11 + e22 + abs(2*a*e12 - e11 + e22))/(2*e12);
% s22neg = s12*(-e11 + e22 - abs(2*a*e12 - e11 + e22))/(2*e12);
% nu = Epos * e12 / s12 - 1 ; 


    







%exp = [7 8 9 10 11 12 13 14 15 16 17 18 20 21 22 23 24 25 26 27 30 31 32];
% for G = 1:length(exp);
%     
%     matpath = sprintf('%s\\Martin_Experiments\\AAA_TensionTorsion\\TT2-%d\\TT2-%d_MatlabResults\\',curdir,exp(G),exp(G)); 
%     
%     fid = fopen( sprintf( '%smax.dat', matpath) , 'r');
%     d = cell2mat(textscan(fid,'%f %f %f %f %f %f %f %f','headerlines',1,'delimiter',' '));
%     fclose(fid)
%     if (alpha(xpts == exp(G)) >= 3.0) || isnan(alpha(xpts == exp(G))) ;
%         plot(d(:,3))
%         title('\SIGMA');
%     else
%         plot(d(:,4))
%         title('\TAU');
%     end
%         
%     %set(gcf,'position',[300 100 600 1000])
%     
%     clicks = ginput(1);
%     stgno = round(clicks(1));
%          
%     ylds = [ylds;stgno d(stgno,3:4) ];
%         
% end;