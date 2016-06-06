#!/usr/bin/env perl

# perl "C:\Students\Martin\Martin_Experiments\DIC processing files\post_3D_DIC_linear_new.pl"

#NAME OF THE FILES
#######################################################"
#INPUT_FILE
$prefixe_in = 'F:\zzMISC\TORSION_PE_unchanged_from_Nico\EXPERIMENT\DC15_TT2_17_alpha_1.0\Rupture\TT2_17_DC15_FS_19_SS6-Stage-0-';
$suffixe_in = '.txt';
#OUTPUT_FILE
mkdir ('F:\Martin_Experiments\TT2-17\post') or die ("That folder path does not exist\n");
$prefixe_ou = 'F:\Martin_Experiments\TT2-17\post\TT2-17_OLD_';
$suffixe_ou = '.dat';
open(OUT_force, '>F:\Martin_Experiments\TT2-17\post\time_force.dat');
open(UNIT, '>F:\Martin_Experiments\TT2-17\post\units.dat');



#Number of pictures
$nb_im =609; 



#calibration factor
#######################################################
#Force
$cal=5000; #1V=5000lbf
    
#Programme
########################################################

print UNIT "time_force.dat \n stage [] time [s] Force [lbf] \n\n\n\n\n 
$prefixe_ou$nb_im$suffixe_ou \n
index i [] index [j] undeformed coordinates X Y Z [mm] displacement U V W ur=w [mm] Surface gradient transformation F(0,0) F(0,1) F(1,0) F(1,1) [] 

X,Y,Z,U,V,W refer to the global coordinate system
ur refer to the undeformed  local coordinate system
F refer to the local coordinate system with z' normal to the surface, x' orthognal to z' and in the plane parallel to XY plane,  y' orthonormal coordinate system";



    
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
	elsif (/^\s{1,}(\d{1,})\s{1,}(\d{1,})\s{1,}([\+\-\.Ee0-9]{1,})\s{1,}([\+\-\.Ee0-9]{1,})\s{1,}([\+\-\.Ee0-9]{1,})\s{1,}([\+\-\.Ee0-9]{1,})\s{1,}([\+\-\.Ee0-9]{1,})\s{1,}([\+\-\.Ee0-9]{1,})\s{1,}([\+\-\.Ee0-9]{1,})\s{1,}([\+\-\.Ee0-9]{1,})\s{1,}([\+\-\.Ee0-9]{1,})\s{1,}([\+\-\.Ee0-9]{1,})\s{1,}([\+\-\.Ee0-9]{1,})$/)
	{
	    $indexi = $1;
	    $indexj = $2;
	    $X = $3;
	    $Y = $4;
	    $Z = $5;
	    $u = $6;
	    $v = $7;
 	    $w = $8;
	    $ur = $9;
	    $F00 = $10;
	    $F01 = $11;
	    $F10 = $12;
	    $F11 = $13;
	    print SORTIE "$indexi $indexj $X $ Y $Z $u $v $w $ur $F00 $F01 $F10 $F11\n";
	}
    }
}
