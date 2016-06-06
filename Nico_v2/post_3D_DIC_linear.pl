# This was Nico's most recent Perl file to take the raw Aramis output files, remove missing lines, and create step_time_force.
# My function "AramisScan.m" does this now.
# This requires the use of Nico's export file format

#!/usr/bin/env perl

# perl "E:\MATLAB\post_3D_DIC_linear.pl"

#NAME OF THE FILES
#######################################################"
#INPUT_FILE
$prefixe_in = 'E:\zzARAMIS Export Files\TT2-25_DC16_FS32_SS8\TT2-25_DC16_FS32_SS8-Stage-0-';
$suffixe_in = '.txt';
#OUTPUT_FILE
mkdir ('E:\Martin_Experiments\TT2-25_Dec3_a3p25\Anisotropy_perl_post') or die ("The path where you're trying to put your new folder doesn't exist\n Or it already exists");
$prefixe_ou = 'E:\Martin_Experiments\TT2-25_Dec3_a3p25\Anisotropy_perl_post\TT2-25_DC16_FS32_SS8_';
$suffixe_ou = '.dat';
open(OUT_force, '>E:\Martin_Experiments\TT2-25_Dec3_a3p25\Anisotropy_perl_post\time_force.dat');
open(UNIT, '>E:\Martin_Experiments\TT2-25_Dec3_a3p25\Anisotropy_perl_post\units.dat');

#$prefixe_ou = 'E:\Martin\TORSION_PE_unchanged_from_Nico\RESULTS_FOR_FAILURE\DC15_TT2-20_alpha_0.5\post\TT2_20_DC15_FS19_SS6_';
#$suffixe_ou = '.dat';
#open(OUT_force, ">post/time_force.dat");
#open(UNIT, ">post/units.dat");


#Number of pictures
$nb_im =170; 



#calibration factor
#######################################################
#Force
$cal=1; #Do not put any calibration factor, it is done after during the matlab file
    
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
