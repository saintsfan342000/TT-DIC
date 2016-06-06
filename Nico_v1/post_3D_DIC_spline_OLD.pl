#!/usr/bin/env perl

#NAME OF THE FILES
#######################################################"
#INPUT_FILE
$prefixe_in = 'E:\Martin\DIC results\TT2_9_DC14_FS19_SS6_spline\TT2_9_DC14_FS19_SS6_spline-Stage-0-';
$suffixe_in = '.txt';
#OUTPUT_FILE
mkdir ('C:\Students\Martin\Martin_Experiments\TT2-9_June27_a2p0\TT2-9_19-6_spline\post') or die ("Erreur creation repertoire\n");
$prefixe_ou = 'C:\Students\Martin\Martin_Experiments\TT2-9_June27_a2p0\TT2-9_19-6_spline\post\TT2-9_19-6_spline_';
$suffixe_ou = '.dat';
open(OUT_force, '>C:\Students\Martin\Martin_Experiments\TT2-9_June27_a2p0\TT2-9_19-6_spline\post\time_force.dat');
open(UNIT, '>C:\Students\Martin\Martin_Experiments\TT2-9_June27_a2p0\TT2-9_19-6_spline\post\units.dat');



#Number of pictures
$nb_im =467; 



#calibration factor
#######################################################
#Force
$cal=5000; #1V=5000lbf
    
#Programme
########################################################

print UNIT "time_force.dat \n stage [] time [s] Force [lbf] \n\n\n\n\n 
$prefixe_ou$nb_im$suffixe_ou \n
index i [] index [j] undeformed coordinates X Y Z [mm] displacement U V W ur=w [mm] logarithmic plastic strain [] 

X,Y,Z,U,V,W refer to the global coordinate system
ur refer to the undeformed  local coordinate system
e^p";



    
for ($i=0 ; $i <= $nb_im ; $i++ )
{
#    $im = substr("0" x $nb . "$i", -$nb);  #Rajoute des zéros en début de nombre s'il en manque
    $arg =  "$prefixe_in$i$suffixe_in \n";
    print "$arg\n";
    open(DONNEE,"$arg") || die ("Erreur d'ouverture de DONNEE") ;
    open(SORTIE, ">$prefixe_ou$i$suffixe_ou");
    
    while (<DONNEE>){
	if (/^#\s{1,}Time\s{1,}:\s{1,}undeformt:\s{1,}([\+\-\.Ee0-9]{1,})\s{1,}\[seconds\s{1,}since\s{1,}01\.01\.1970\]$/)
	{
#	    print "bouep \n";
	    $t0 = $1;
	}
	elsif (/^#\s{1,}deformt:\s{1,}([\+\-\.Ee0-9]{1,})\s{1,}\[seconds\s{1,}since\s{1,}01\.01\.1970\]$/)
	{
	    $t = $1-$t0;
	    print OUT_force "$i $t $F\n"
	}
	elsif (/^#\s{1,}AD\-0:\s{1,}([\+\-\.Ee0-9]{1,})\s{1,}([\+\-\.Ee0-9]{1,})\s{1,} V$/)
	{
	    $F = $2*$cal;
	}
	elsif (/^\s{1,}(\d{1,})\s{1,}(\d{1,})\s{1,}([\+\-\.Ee0-9]{1,})\s{1,}([\+\-\.Ee0-9]{1,})\s{1,}([\+\-\.Ee0-9]{1,})\s{1,}([\+\-\.Ee0-9]{1,})\s{1,}([\+\-\.Ee0-9]{1,})\s{1,}([\+\-\.Ee0-9]{1,})\s{1,}([\+\-\.Ee0-9]{1,})$/)
	{
	    $indexi = $1;
	    $indexj = $2;
	    $X = $3;
	    $Y = $4;
	    $Z = $5;
	    $u = $6;
	    $v = $7;
 	    $w = $8;
	    $EP = $9;
	    print SORTIE "$indexi $indexj $X $ Y $Z $u $v $w $EP\n";
	}
    }
}
