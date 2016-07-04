#!/usr/bin/perl -w

# Copyright (c) 2006,2007,2008,2009,2012 
# by Christian Panse <Christian.Panse@fgcz.uzh.ch>, Bertran Gerrits <gerritsb@gmail.com>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.

=head1 NAME 

mascotDat2RData.pl - mascot dat file to RData exporter

=head1 SYNOPSIS

mascotDat2RData.pl [-m=<mascot mod-file>] [-d=<mascot dat-file>] 

=head1 COPYRIGHT

Copyright (c) 2012 Christian Panse.

GNU General Public License Version 3

=head1 AUTHORS

Christian Panse <cp@fgcz.ethz.ch>

Bertran Gerrits

=head1 DESCRIPTION

The program exports Matrix Science (http://www.matrixscience.com/) 
mascot dat files to an R (http://www.r-project.org/) object.

The program reqires R install and is testet on a debian linux system.

The program is part of the peakplot code. 

=head1 OPTIONS

=head2 -m, --modfilename

Sets the location of the Matrix Science Mascot Server dat 'mod_file'
Usually it is located under '/usr/local/mascot/mascot-2_2/config/mod_file' Mascot Server.

=head2 -d, --datfilename

Defines the location of the Mascot Server 'dat' file.

=cut

use strict;
use warnings;
use File::Copy;
use File::Basename;
use File::Path;
use POSIX qw( strftime );
use POSIX;
# use Search::Binary;

my $version;
$version="1.00; Wed Nov 28 14:10:08 CET 2012";

my $datfilename="";
my $datfilename_="";
my $modfilename="";
my $peakplotTitle="";

my $extractScanNumber=0;
my $extractNS2 = 1;

my $outputdir="/tmp/peakplot";
my $patternFile="";


my %FIXEDMOD;
my %DAT;

my @VARIAMOD;
my @DATFILE;
my @MODFILE;
my @MZPATTERN;

my @LABEL=("1","2","3","4","a","a*","a_0","b","b*","b_0","c","a'","y","y*","y_0","z+1","z+2","b-98","b*-98","b_0-98","y-98","y*-98","y_0-98","c+58", "z", "z-57", "z+1-57", "z+2-57");

my @LABELPRINT=("","","","","a","a*","a_0","b","b*","b_0","c","a'","y","y*","y_0","z+1","z+2","","","","","","","c+58", "z", "z-57", "z+1-57", "z+2-57");

my @LABELPRINT_CID=("","","","","a","a*","a_0","b","b*","b_0","","","y","y*","y_0","","","","","","","","","", "", "", "", "");
my @LABELPRINT_ETD=("","","","","","","","","","","c","","y","","","z+1","z+2","","","","","","","", "z", "", "", "");

my @LABELPRINT_ETD_ISOASP=("","","","","","","","","","","c","","","","","z+1","z+2","","","","","","","c+58", "z", "z-57", "z+1-57", "z+2-57");
my @LABELPRINT_CID_PHOSPHO=("","","","","a","a*","a_0","b","b*","b_0","","a'","y","y*","y_0","","","b-98","b*-98","b_0-98","y-98","y*-98","y_0-98","", "", "", "", "");


# 26 Rules
my @BYRULES=(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0);
my @M =(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0);

my $C_; 
my $N_; 
my $O_; 
my $H_; 

# dat file
my $N_term_;
my $C_term_;

# mod_file
my $ProteinNterm=0;
my $ProteinCterm=0;

# Neutral loss for Phosphorylation
my $NeutralLossP=1;
my $P_=30.9737620;
my $Pavg_=30.97;
my $H3PO4_;

#sub fabs {
#my $f = shift;
#if ($f < 0) {
#    return -$f;
#}
#return $f;
#}

sub round {
    my($number) = shift;
    $number*=100;
    my $rNumber=int($number + .5);
    $rNumber/=100;

    return $rNumber;
}


sub peaklistsort {
    my @a_ = split /:/, $a;
    my @b_ = split /:/, $b;
    $a_[0] <=> $b_[0]
}

sub byTableCmp0 {
    my @a_ = split /\ /, $a;
    my @b_ = split /\ /, $b;
    $a_[8] <=> $b_[8]
}

sub byTableCmp1 {
    my @a_ = split /\ /, $a;
    my @b_ = split /\ /, $b;
    $a_[3] cmp $b_[3] || $a_[2] <=> $b_[2] || $a_[9] <=> $b_[9]
}

sub readMZPATTERN(){
    my $fn=shift || die "no filename given $!\n";
    print "-- reading '".$fn."' ...";
    open (FILE, $fn) || die "could not open file '$fn' $!\n:";
    while(<FILE>){
        chomp;
        if (/^(\d+\.\d+)\s*$/){
            push @MZPATTERN, $_;
        }
    }
    print "DONE\n";
    close(FILE);
    if ($#MZPATTERN < 1){
        die "no mods defined $!\n";
    }
    print "-- found ".($#MZPATTERN+1)." patterns\n";
    return(1); 
}


sub init() {
    my @M;
    my $i;

print  "-- reading ".$datfilename.".\n";
open (FILE, $datfilename);
@DATFILE=<FILE>;
foreach(@DATFILE) {
    chomp;
    s/\\$//g;
    if (/^([A-Z])=(\d+\.\d+)$/){
        $FIXEDMOD{$1}=$2;
    }
    elsif (/^(.+)=(.+)$/){
        if(!exists($DAT{$1})) {
            $DAT{$1}=$2;
        }
    }
}
close(FILE);

    print  "-- setting rule set.\n";
    print "-- ";
    foreach (split/,/, $DAT{"RULES"}) {
        $BYRULES[$_]=1;
    }
    foreach ($i=1; $i<$#BYRULES; $i++) {
        printf  (" %2d", $i);
    }
    print "\n-- ";
    foreach ($i=1; $i<$#BYRULES; $i++) {
        print  "  ".$BYRULES[$i];
    }
    print "\n";

print  "-- reading ".$modfilename.".\n";
    if ( -s $modfilename ){
        open (FILE, $modfilename) || die "could not open modfile $!\n";
        @MODFILE=<FILE>;
        close(FILE);
    }

print  "-- setting variable modification\n";
if (exists($DAT{"IT_MODS"})) {
    @M=split /,/, $DAT{"IT_MODS"};
    for ($i=0; $i <= $#M; $i++) {
        &get_varia_modification($M[$i]);
    }
}

if (exists($DAT{"MODS"})) {
    print  "-- setting fixed modification\n";
    @M=split /,/, $DAT{"MODS"};
    for ($i=0; $i<=$#M; $i++) {
        &get_fixed_modification($M[$i]);
    }
}

print  "-- setting Carbon, Nitrogen, Hydrogen, Oxygen\n";
$C_=$DAT{"Carbon"} || die "could not find Carbon definion.\n";
$N_=$DAT{"Nitrogen"};
$H_=$DAT{"Hydrogen"};
$O_=$DAT{"Oxygen"};

$C_term_=$DAT{"C_term"};
$H3PO4_= (3 * $H_) + $P_ + (4 * $O_);
}

sub get_varia_modification() {
my $m=shift || die "no argument given $!.\n";
my %h;
my $i;
my $modname="";

for ($i=0; $i <= $#MODFILE; $i++) {
    if ($MODFILE[$i] =~ /^Title:(.+)/) {
        if ($1 eq $m) {
            print  "\t->".$1."\n";
            $modname=(split /[\s\_]/, $1)[0];
            for (my $j=0; $MODFILE[$i+$j] !~ /^\*/; $j++) {
                if ($MODFILE[$i+$j] =~ /^Residues[a-zA-Z]*:([A-Z])\ (-*\d+\.\d+)\ (-*\d+\.\d+)/) {
                    print  "\t\t->".$1." ".$2."\n";
                    $h{$1}=$2;
                    $h{$1."_name"}=$modname;
                }
                elsif ($MODFILE[$i+$j] =~ /^ProteinNterm:(-*\d+\.\d+)\ (-*\d+\.\d+)/) {
                    print  "\t\t->".$1." ".$2."\n";
                    $ProteinNterm=$1;
                }
                elsif ($MODFILE[$i+$j] =~ /^ProteinCterm:(-*\d+\.\d+)\ (-*\d+\.\d+)/) {
                    print  "\t\t->".$1." ".$2."\n";
                    $ProteinCterm=$1;
                }
            }
            if (keys(%h)) {
                push @VARIAMOD, { %h };
                return;
            } 
            else {
                push @VARIAMOD, { %h };
                return;
            }
        }
    }
}
}

sub get_fixed_modification() {
my $m=shift || die "no argument given $!.\n";
my $i;

for ($i=0; $i <= $#MODFILE; $i++) {
    if ($MODFILE[$i] =~ /^Title:(.+)/) {
        if ($1 eq $m) {
            print  "\t->".$1."\n";
            for (my $j=0; $MODFILE[$i+$j] !~ /^\*/; $j++) {
                if ($MODFILE[$i+$j] =~ /Residues:([A-Z])\ (\d+\.\d+)\ (\d+\.\d+)/) {
                    print  "\t\t->".$1." ".$2."\n";
                    $FIXEDMOD{$1}=$2;
                }
            }
            return;
        }
    }
}
}

sub getTitle() {
    my $q=shift || die "no argument given $!.\n";
    my ($key, $info);
    my $title="NA";
    my ($scannumber, $mzXMLFilename) = (-1, "NA");
    if ( $q=~/q(\d+)_p\d+/ || $q=~/(\d+)/ ) {
        my $n=$DAT{"query".$1}; # compute line number
        for (my $i=0; $DATFILE[$n+$i] !~ /^--/; $i++) {
            if ($DATFILE[$n+$i]  =~ /(^[a-z]+)=(.+)$/){
                $key = $1;
                $info = $2;
                chomp $key;
                chomp $info;
                if ( $key =~ /title/ ){
                    $info =~ s/%2d/_/g;
                    $info =~ s/%5c/\//g; 
                    $info =~ s/%2e/\./g;
                    $info =~ s/\\/\//g;
                    $info =~ s/%5d//g;
                    $info =~ s/%3a/:/g; 
                    $info =~ s/%20/\ /g;
                    $info =~ s/%28/(/g; 
                    $info =~ s/%3d/=/g; 
                    $info =~ s/%29/)/g; 
                    $info =~ s/%5b//g;
                    return ($info);
                }
            }
        }
    }
    return ($title);
}
sub getRtinseconds() {
    my $q=shift || die "no argument given $!.\n";
    my ($key, $info);
    my ($scannumber, $mzXMLFilename) = (-1, "NA");
    if ( $q=~/q(\d+)_p\d+/ || $q=~/(\d+)/ ) {
        my $n=$DAT{"query".$1}; # compute line number
        for (my $i=0; $DATFILE[$n+$i] !~ /^--/; $i++) {
            if ($DATFILE[$n+$i]  =~ /(^rtinseconds)=(.+)$/){
                return($2);
            }
        }
    }
    return($scannumber);
}
sub getScans() {
    my $q=shift || die "no argument given $!.\n";
    my ($key, $info, $scans);
    my ($scannumber, $mzXMLFilename) = (-1, "NA");
    if ( $q=~/q(\d+)_p\d+/ || $q=~/(\d+)/ ) {
        my $n=$DAT{"query".$1}; # compute line number
        for (my $i=0; $DATFILE[$n+$i] !~ /^--/; $i++) {
            if ($DATFILE[$n+$i]  =~ /(^scans)=(.+)$/){
                $scans=$2;
                $scans =~ s/-/:/g;
                return($scans);
            }
        }
    }
    return($scannumber);
}
sub get_scanNumber() {
    my $q=shift || die "no argument given $!.\n";
    my $f=shift || die "no argument given $!.\n";
    my ($key, $info);
    my ($scannumber, $mzXMLFilename) = (-1, "NA");
    if ( $q=~/q(\d+)_p\d+/ || $q=~/(\d+)/ ) {
        my $n=$DAT{"query".$1}; # compute line number
        for (my $i=0; $DATFILE[$n+$i] !~ /^--/; $i++) {
            if ($DATFILE[$n+$i]  =~ /(^[a-z]+)=(.+)$/){
                $key = $1;
                $info = $2;
                chomp $key;
                chomp $info;
                if ( $key =~ /title/ ){
                    $info =~ s/%2d/_/g;
                    $info =~ s/%5c/\//g; 
                    $info =~ s/%2e/\./g;
                    $info =~ s/\\/\//g;
                    $info =~ s/%5d//g;
                    $info =~ s/%3a/:/g; 
                    $info =~ s/%20/\ /g;
                    $info =~ s/%28/(/g; 
                    $info =~ s/%3d/=/g; 
                    $info =~ s/%29/)/g; 
                    $info =~ s/%5b//g;
                    $mzXMLFilename = &basename($info);
                    if($info =~ /.+Scan\s([0-9]+)\s.+/){
                        $scannumber=$1;
                    }
                } elsif ( $key =~ /scans/) {
                    $scannumber=$info;
                }
#elsif ( $key =~ /index/ && $scannumber < 0) {
#                   $scannumber="INDEX=".$info;
#               }
            }
        }
        print basename($f)."\t".$q."\t".$scannumber."\t".$mzXMLFilename."\n";
    }
    return($scannumber."\t".$mzXMLFilename);
}

sub get_peaklist_mz() {
    my $q=shift || die "no argument given $!.\n";
    my @MZ;
    foreach (&get_peaklist($q)){
        my ($mz, $intensity) = split /:/;
        push @MZ, $mz;
    }
    return @MZ;
}
sub get_peaklist() {
    my $q=shift || die "no argument given $!.\n";

    if ( $q=~/q(\d+)_p\d+/ || $q=~/(\d+)/ ) {
        my $n=$DAT{"query".$1} || die "can not find spectrum index for query $1 \n"; # compute line number

        for (my $i=0; $DATFILE[$n + $i] !~ /^--/; $i++) {
            if ($DATFILE[$n+$i] =~ /^Ions1=(\d+.+)$/) {
                return sort peaklistsort split /,/, $1;
            }
        }
    }
}

sub compute_by_phosphotylation {
my $b = shift;
my $y = shift;
my $a = ($b - $C_ - $O_);
my @R;

# rules 25
if ($NeutralLossP == 1) {
    push @R, ($b-$H3PO4_); 
} else {
    push @R, 0; 
}
# rules 26
if ($NeutralLossP == 1) {
    push @R, ($b - $N_ - (3 * $H_) - $H3PO4_); 
} else {
    push @R, 0; 
}

# rules 27
if ($NeutralLossP == 1) {
    push @R, ($b - $O_ - (2 * $H_) - $H3PO4_); 
} else {
    push @R, 0; 
}

# rules 28
if ($NeutralLossP == 1) {
    push @R, ($y - $H3PO4_); 
} else {
    push @R, 0; 
}

# rules 29
if ($NeutralLossP == 1) {
    push @R, ($y - $N_ - (3 * $H_) - $H3PO4_); 
} else {
    push @R, 0; 
}
# rules 30
if ($NeutralLossP == 1) {
    push @R, ($y - $O_ - (2 * $H_) - $H3PO4_); 
} else {
    push @R, 0; 
}
return @R;
}


sub compute_by {
my $b = shift;
my $y = shift;

my @R;

my $a = $b - $C_ - $O_;
my $c = $b + $N_ + (3 * $H_);

my $z = $y - ($N_ + (3 * $H_));

# rules 1-3
push @R, 0;
push @R, 0;
push @R, 0;

# rules 4
push @R, 0;

# rules 5 
push @R, $a;

# rules 6
push @R, ($a - $N_ - (3 * $H_));

# rules 7
push @R, ($a - $O_ - (2 * $H_));

# rules 8
push @R, $b; 

# rules 9 
push @R, ($b - $N_ - (3 * $H_));

# rules 10 
push @R, ($b - $O_ - (2 * $H_));

# rules 11 (c)
push @R, $c;

# rules 12.1 (a°)
push @R, ($b + ($N_ + 3 * $H_) - 44);

# rules 13
push @R, $y;

# rules 14
push @R, ($y - $N_ - (3 * $H_));

# rules 15
push @R, ($y - $O_ - (2 * $H_));

# rule 21 (z+1)
push @R, ($z  + $H_);

# rule 25 (z+2)
push @R, ($z + (2*$H_));

# rules 25
if ($NeutralLossP == 1) {
    push @R, ($b-$H3PO4_); 
} else {
    push @R, 0; 
}
# rules 26
if ($NeutralLossP == 1) {
    push @R, ($b - $N_ - (3 * $H_) - $H3PO4_); 
} else {
    push @R, 0; 
}

# rules 27
if ($NeutralLossP == 1) {
    push @R, ($b - $O_ - (2 * $H_) - $H3PO4_); 
} else {
    push @R, 0; 
}

# rules 28
if ($NeutralLossP == 1) {
    push @R, ($y - $H3PO4_); 
} else {
    push @R, 0; 
}

# rules 29
if ($NeutralLossP == 1) {
    push @R, ($y - $N_ - (3 * $H_) - $H3PO4_); 
} else {
    push @R, 0; 
}
# rules 30
if ($NeutralLossP == 1) {
    push @R, ($y - $O_ - (2 * $H_) - $H3PO4_); 
} else {
    push @R, 0; 
}

# rule 31 (c+58)  C2H2O2
push @R, ($c + ($C_*2 + $H_*2 + $O_*2));    

# rule 32 (z)
push @R, ($z);    

# rule 33 (z-57) C2HO2
push @R, ($z - ($C_* 2 + $H_ + $O_ * 2));    

# rule 34 (z+1-57) C2HO2
push @R, ($z + $H_ - ($C_*2 + $H_ + $O_ * 2));    

# rule 35 (z+2-57) C2HO2
push @R, ($z + ($H_ * 2) - ($C_*2 + $H_ + $O_ * 2));    


return @R;
}

sub writeBY() {
my @BY=sort byTableCmp0 @{$_[0]};
my $query=$_[1];

my %h;

my ($k, $k2, $label, $insilicomass, $mass, $intensity);
my @BY1;
my @BY2;

foreach (@BY){
    chomp;
    $label = (split /\ /)[2];
    $insilicomass = (split /\ /)[3];

    $mass = (split /\ /)[6];
    $intensity = (split /\ /)[7];

    $k=$label."_".$insilicomass;
    $k2=$mass."_".$intensity;

    if (exists $h{$k} || exists $h{$k2}) {
        push @BY1, "F ".$_;
    } else {
        $h{$k}=$_;
        $h{$k2}=$_;
        push @BY1, "T ".$_;
    }
}

@BY2 = sort byTableCmp1 @BY1;

open(FILE, ">".$outputdir."/by.".$query.".txt");
print FILE "'VALID'\t'AA'\t'LPOS'\t'LABEL'\t'INSILICOMASS'\t'CHARGE'\t'DELTA'\t'MASS'\t'INTENSITY'\t'ERROR'\t'QUERY'\n";
foreach (@BY2){
    chomp;
    s/\ /\t/g;
    print FILE $_;
    print FILE "\n";
}
close(FILE);
}

sub MonoisotopicAAmass4R(){
    my ($s, $modification, $aa, $rule);
    $s = shift || die "no peptide string given.\n";
    $modification = shift || die "no modification string fiven\n";

    my @A=split //, $s;
    my @M1=split //, $modification;

    my $aaWeight=0.0;
    my $res="c(";

    for (my $i=0; $i <= $#A; $i++) {

        $aa=$A[$i];
        $aaWeight = 0.0;

        if ( $i < $#M1 && $M1[$i+1] =~ /[1-9]/ ){

            $rule = $M1[$i+1] - 1;

            if( exists ($VARIAMOD[$rule]{$aa}) ){
                    $aaWeight += $VARIAMOD[$rule]{$aa} - $FIXEDMOD{$A[$i]};
            }else{

                    print $s . "\t" . $modification . "\n";
                    print "error: VARIAMOD[" . $rule . "]{" . $aa . "} could not be found. check mascot dat and mod_file!\n";
            }

        } else {
            $aaWeight = $FIXEDMOD{$A[$i]};
        }
        $res .= $aaWeight;
        if ($i < $#A){
                $res.=",";
        }
    }
    $res.=")";
    return ($res);
}



=head2 searchPattern()

To find the mZ pattern we iterrate over each MS2 scan and save the precorusor mass in a hash list.
In a second run we try to match the input pattern in each spectrum. On each MS2 pattern match we consider all
all list items in the 



Peakplot retrospectively labels 
the spectra from a peptide sequence assignments by the Mascot 
search algorithm with the appropriate fragment ion labels. The 
application uses Perl for the manipulation of data and R for 
label heuristics and plotting.


The searchPattern function iterates throw all MS2 type peaklists and find a match, 
e.g. a newline sep. file of mZ values Glyco Pattern 

./peakplot_glyco.pl -d=F165616.dat -m=mod_file --pattern=p1130--pattern.txt

NOTE: we assume we have the same number of MS2 as in the MGF even 
if there is no usefull pepseq assignment;
because we can do the 'assignment by a query shift' 

2012-01-20 by PN,CP
TODO: 

* add SEEN HASH

* move the PATTERFILE stuff into the ms2 class

INPUT:

* Mass Spectrometric Measurement

* MS2 m/z pattern

OUTPUT:

* mgf fie containing the MS2 having the given pattern included plus the next MS2 scan 

* all the peakplist and composed peptide information for the R peakplot plotting code

=cut back to the main
sub mascotDAT2R() {

    my $RFILE = shift || die "no file descriptor given.\n";
    my $datfilename = shift || "no data structur name given.\n";

    my $mPC;
    my @MSM;
    my $n=$DAT{'queries'};

    my ($peptideSequence, $mascotScore, $modification, $proteinInformation, $start, $end);


    ## COMPUTING
    for my $i(1 .. $n){
        $mPC=&getPrecursor($i);
        my $ms2 = eval { new FGCZMS2(); } or die "calling constructor failed  $?\n"; 
        $ms2->setTitle(&getTitle($i));
        $ms2->setCharge(&getCharge($i));
        $ms2->setScan(&getScans($i));
        $ms2->setRtinseconds(&getRtinseconds($i));
        $ms2->setPeaklist([ &get_peaklist($i) ]);
        $ms2->setPepMass($mPC);
        $ms2->setQuery($i);

        # push @MSM, $ms2;
        
        print RFILE $datfilename."[[".$i."]] <- list(\n";
        my $query="q".$i."_p1";          

        $peptideSequence="NA";
        $mascotScore="NA";
        $proteinInformation="NA";
        $proteinInformation="NA";
        $modification="NA";
       
        if ( exists ($DAT{$query}) ) {
            my @Q=split /,/, $DAT{$query};

            if ($#Q > 1){
                $peptideSequence=$Q[4];
                $mascotScore=$Q[7];
                $modification=$Q[6];
                $proteinInformation=(split /\"/, $Q[$#Q])[1] || die "split proteinInformation failed $!\n";

                ##### TODO
                # my $proInfo=split /;/, $DAT{$query};

                $proteinInformation =~ s/"/\\"/g;
            }
        }

        print RFILE "id=".$i.",\n";
        print RFILE "peptideSequence='".$peptideSequence."',\n";
        print RFILE "mascotScore=".$mascotScore.",\n";
        print RFILE "searchEngine='mascot',\n";
        print RFILE "modification='".$modification."',\n";
        print RFILE "MonoisotopicAAmass=".&MonoisotopicAAmass4R($peptideSequence, $modification).",\n";
        print RFILE "proteinInformation='".$proteinInformation."'";
        if ($extractNS2 == 1){
            print RFILE ",\n";
            $ms2->printR(\*RFILE);
        }
        print RFILE ")\n\n";

    } # END foreach query
}

sub getCharge(){
    my $queryNumer=shift || die "no argument given $!\n";
    my  $qexp=(split /=/, $DAT{"qexp".$queryNumer})[0];

    my ($m,$c) = split /,/, $qexp;

    $c =~ 's/\+//g';
    return ($c);
}


sub getPrecursor(){
    my $queryNumer=shift || die "no argument given $!\n";
    my  $qexp=(split /=/, $DAT{"qexp".$queryNumer})[0];
    
    my ($m,$c) = split /,/, $qexp;
   
    return ($m);
}

sub usage {
    die "usage:\n$0 -d=<mascot dat file> -m=<mascot mod_file>\n\n";
}

sub main {
    my $now = strftime("%Y%m%d-%H%M%S", localtime);
    print "-- starting ".$0." at ".$now."\n";

    my $a;
    while ($a = shift @ARGV) {
        chomp $a;

        if ($a =~ /-d=(.+)/ || $a =~ /--datfilename=(.+)/) {
            $datfilename=$1;
            my @TMP=split /\//, $1;
            $datfilename_=$TMP[$#TMP];
        }
        elsif ($a =~ /-m=(.+)/ || $a =~ /--modfilename=(.+)/) {
            $modfilename=$1;
        }
        else {
            &usage;
        }
    }

    if (-f $modfilename && -f $datfilename) {
        if (! -d $outputdir ){
            mkdir $outputdir || die "could not create ".$outputdir.".\n";
        }

        my $RdataFileName=&basename($datfilename,".dat");
        open (RFILE, ">/tmp/dat.R") || die "could not open file for writing ...$!\n";
        #open (RFILE, "| tee /tmp/dump.R | R --no-save --slave") || die "could not open file for writing ...$!\n";

        &init;

        my $modcount=0;
        my ($modDesc, $modMass);

       for my $varmod(@VARIAMOD){

               $modcount++;

                for my $rule(keys %$varmod){
                    print $modcount . "\t" . $varmod->{$rule} . "\t" . $rule ."\n";

                    if ($varmod->{$rule} =~ /\d+\.\d+/){
                        if ($modcount > 1){
                        $modMass .= ", " ;
                        }
                        $modMass .= $varmod->{$rule};
                    }else{
                        if ($modcount > 1){
                            $modDesc .= ", " ;
                        }
                        $modDesc .= "'".$varmod->{$rule}."'" ;
                    }
                }
       }

        print RFILE "#R\n". $RdataFileName .".modification.description <- c('', ". $modDesc . ");\n";
        print RFILE $RdataFileName .".modification.mass <- c(0.0,". $modMass . ");\n";

        print RFILE $RdataFileName ." <- list()\n";
        &mascotDAT2R(\*RFILE, $RdataFileName);

        print RFILE "save(list=c('" . $RdataFileName . "','" . $RdataFileName . ".modification.mass', '" . $RdataFileName . ".modification.description') , file='" . $RdataFileName . ".RData', compress=TRUE)\n";

        close(RFILE);
    }
    else {
        &usage;
    }
    my $finish = strftime("%Y%m%d-%H%M%S", localtime);
    print "-- ending ".$0." at ".$finish."\n";
}

&main;

exit 0;

package FGCZMS2;

use strict;
use POSIX;

sub new {
        my ($class) = @_;
        my $self = {
            _title => undef,
            _pepmass => undef,
            _charge => undef,
            _scan => undef,
            _rtinseconds => undef,
            _patternsumintensity => 0.0,
            _pattern => [ ],
            _peaklist => undef,
            _query => undef,
            _mz => [ ],
            _intensity => [ ]
        };
        bless $self, $class;
        return $self;
}

sub printR {
        my ($self, $fh) = @_;
        my $count = 0;


        print $fh "title=\"".$self->{_title}."\",\n";
        print $fh "pepmass=".$self->{_pepmass}.",\n";
        print $fh "charge=".$self->{_charge}.",\n";
        print $fh "scans=".$self->{_scan}.",\n" if (defined ($self->{_scan}));
        print $fh "rtinseconds=".$self->{_rtinseconds}.",\n" if (defined ($self->{_rtinseconds}));

        my @INTENSITY;
        my @MZ;
        my $n;

        foreach my $peak(@{$self->{_peaklist}}){
            chomp $peak;
            # print $peak . "\n";
            $count++;
            my ($mZ, $intensity) = split /:/, $peak;

            if (!($peak =~ m/^\w+$/)){
                push @INTENSITY, $intensity;
                push @MZ, $mZ;
            }
        }

        # exit(0);
        $n = $#INTENSITY;

        print $fh "mZ=c(";
        for (my $i=0; $i < $n; $i++){
            print $fh $MZ[$i].","; 
        }
        print $fh $MZ[$n]."),\n"; 

        print $fh "intensity=c(";
        for (my $i=0; $i < $n; $i++){
            print $fh $INTENSITY[$i].","; 
        }
        print $fh $INTENSITY[$n].")\n"; 

}

sub setTitle {
        my ( $self, $x ) = @_;
        $self->{_title} = $x if (defined($x));
        return ($self->{_title});
}

sub setPepMass{
        my ( $self, $x ) = @_;
        $self->{_pepmass} = $x if (defined($x));
        return ($self->{_pepmass});
}

sub getPepMass{
        my ( $self ) = @_;
        return ($self->{_pepmass});
}

sub setCharge {
        my ( $self, $x ) = @_;
        $self->{_charge} = $x if (defined($x));
        $self->{_charge} =~ s/\+//g;
        return ($self->{_charge});
}

sub getScan {
        my ( $self ) = @_;
        return ($self->{_scan});
}

sub setScan {
        my ( $self, $x ) = @_;
        $self->{_scan} = $x if (defined($x));
        return ($self->{_scan});
}

sub getRtinseconds {
        my ( $self ) = @_;
        return ($self->{_rtinseconds});
}

sub setRtinseconds {
        my ( $self, $x ) = @_;
        $self->{_rtinseconds} = $x if (defined($x));
        return ($self->{_rtinseconds});
}

sub getPatternsumintensity{
        my ( $self ) = @_;
        return ($self->{_patternsumintensity});
}

sub getPeak {
        my ( $self, $idx ) = @_;
        my $peak=@{$self->{_peaklist}}[$idx];
        my ($mZ, $intensity) = split /:/, $peak;
        return (($mZ, $intensity));
}

sub setPatternsumintensity{
        my ( $self, $x ) = @_;
        $self->{_patternsumintensity} = $x if (defined($x));
        return ($self->{_patternsumintensity});
}

sub setPeaklist {
        my ( $self, $x ) = @_;
        my ($mZ, $intensity);
        if (defined($x)){
            $self->{_peaklist} = $x;
            foreach my $peak(@{$self->{_peaklist}}){
                ($mZ, $intensity) = split /:/, $peak;
                push @{$self->{_mz}}, $mZ;
                push @{$self->{_intensity}}, $intensity;
            }
        }
        return ($self->{_peaklist});
}

sub pushPattern {
        my ( $self, $x ) = @_;
        push @{$self->{_pattern}},  $x if (defined($x));
        return ($self->{_pattern});
}

sub setQuery {
        my ( $self, $x ) = @_;
        $self->{_query} = $x if (defined($x));
        return ($self->{_query});
}

sub getQuery {
        my ( $self ) = @_;
        return ($self->{_query});
}

1;

=head1 NAME
package FGCZMGF2MS2 
./peakplot_glyco.pl -d=F165616.dat -m=mod_file --pattern=p1130--pattern.txt
=cut
