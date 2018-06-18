#!/usr/bin/perl -w

# Copyright (c) 2013
# by Christian Panse <cp@fgcz.ethz.ch>
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

protViz_mgf2RData.pl - mascot generic  file to RData exporter

=head1 SYNOPSIS

cat <mascot generic file> | mgf2RData.pl -n=<name if R data object>  

=head1 COPYRIGHT

Copyright (c) 2013 Christian Panse.

GNU General Public License Version 3

=head1 AUTHORS

Christian Panse <cp@fgcz.ethz.ch>

=head1 DESCRIPTION

The program exports Matrix Science (http://www.matrixscience.com/) 
mascot dat files to an R (http://www.r-project.org/) object.

The program reqires R install and is testet on a debian linux system.

The program is part of the protViz R package on CRAN

=head1 OPTIONS

=head2 -m

xxx

=cut

use strict;
use warnings;
use Pod::Usage;

my $message_text  = "This text precedes the usage message.";
my $exit_status   = 2;          ## The exit status to use
my $verbose_level = 0;          ## The verbose level to use
my $filehandle    = \*STDERR;   ## The filehandle to write to

#my $MS2PEAKPATTERN="/^(\d+\.\d+)\s(\d[-\+eE0-9\.]*)$/";

sub main(){
    
    my $Rdata = "mgf";
    $Rdata = shift || die "no data name provided";
    my $i = 1;
    my ($_title, $_pepmass, $_charge, $_scan, $_rtinseconds, @mZ, @intensity);


    open (RFILE, " | tee /tmp/dump.R | R --no-save") || die "could not open file for writing ...$!\n";

    print RFILE $Rdata." <- list()\n";

    while (<>){
        s/\r\n/\n/;
        s/\\/\//g;
        chomp;

        if (/^BEGIN IONS/){
            $_scan='NA';
            $_charge='NA';
            $_pepmass='NA';
            $_title=$i;
            $_rtinseconds='NA';
        }elsif (/^END IONS/){
            #$_charge = "NA";

            print RFILE "\n". $Rdata . "[[" . $i . "]] <- list(\n";
            print RFILE "title=\"" . $_title ."\",\n";
            print RFILE "rtinseconds=" . $_rtinseconds.",\n";
            print RFILE "charge=" . $_charge.",\n";
            print RFILE "scan=" . $_scan.",\n";
            print RFILE "pepmass=" . $_pepmass.",\n";

            push @mZ, 0;
            push @intensity, 0;
            print RFILE "mZ=c(";
            for (my $ii=0; $ii < $#mZ; $ii++){
                print RFILE $mZ[$ii];                    

                print RFILE ", " if ($ii < ($#mZ - 1)); # print a , but at the end
            }
            print RFILE "),\n";

            print RFILE "intensity=c(";
            for (my $ii=0; $ii < $#intensity; $ii++){
                print RFILE $intensity[$ii];                    
                print RFILE ", " if ($ii < $#intensity-1);
            }
            print RFILE ")\n";

            print RFILE ")\n";
            undef $_title;
            undef $_scan;
            undef $_charge;
            undef $_rtinseconds;
            undef $_pepmass;
            undef @mZ;
            undef @intensity;
            $i = $i +1;

        }elsif (/^TITLE=(.+)/){
            $_title = $1;
        }elsif (/^PEPMASS=(\d+\.\d+)\s([-\+eE0-9\.]+)$/){
            $_pepmass = $1;
        }elsif (/^CHARGE=(\d+)./){
            $_charge = $1;
        }elsif (/^SCANS=(\d+)/){
            $_scan = $1;
        }elsif (/^RTINSECONDS=(.+)/){
            $_rtinseconds = $1;
        }elsif (/^(\d+\.\d+)\s(\d[-\+eE0-9\.]*)(\s.+){0,1}$/){
            push @mZ, $1;
            push @intensity, $2;
        }
   }
   print RFILE "save($Rdata , file='" . $Rdata . ".RData', compress=TRUE)\n";
   close(RFILE);
}

## MAIN

my $a;
while ($a = shift @ARGV) {
    chomp $a;
    if ($a =~ /-n=(.+)/ || $a =~ /--name=(.+)/) {
        &main($1);
        exit(0);
    }
}
pod2usage($message_text);

exit(0);

# UNIT TEST
##BEGIN IONS
#TITLE=File26603 Spectrum7633 scans: 6829 [s:\p1000\Proteomics\QEXACTIVE_3\ctrachse_20131203_Hyal_digest\20131203_02_Hyal_trp_digst.raw]
#PEPMASS=916.39435 63427192.00000
#CHARGE=2+
#RTINSECONDS=1670
#SCANS=6829
#120.08123 308270
#126.05535 2.02338E+06
#136.07600 340035
#138.05524 7.37022E+06
#139.05844 382919
#144.06575 537940
#168.06578 2.91089E+06
#175.11929 258069
#186.07643 1.71301E+06
#187.08005 169914
#204.08702 5.36984E+06
#205.09042 338055
#217.09813 34703.1
#243.10995 38117
#283.14462 870741
#284.14822 57171.2
#311.13953 300598
#318.12369 158268
#335.15051 152774
#366.13980 451333
#367.14273 40376.9
#382.17639 52479.7
#390.18008 49542.2
#403.92993 42889.8
#419.17188 204617
#436.19861 548369
#508.22101 171272
#536.21649 310240
#549.28271 182880
#621.30670 48497.6
#649.30060 50656.2
#663.32648 54566.5
#666.32678 150397
#703.32202 80006.2
#720.34485 159669
#774.35883 173619
#791.38373 354543
#792.38403 51205.3
#874.42590 51064.1
#923.42566 229054
#994.46338 423313
#995.47021 247284
#1084.48743 399649
#1085.48718 281613
#1101.51599 1.43512E+06
#1102.51855 684589
#1103.51953 72945.7
#1184.55273 1.0719E+06
#1185.55884 583711
#1197.53857 71304.6
#1213.57800 67498
#1260.57275 216926
#1261.57178 83385.7
#1286.58398 530865
#1287.58398 357716
#1288.57861 225659
#1304.59595 1.02399E+07
#1305.59851 6.31469E+06
#1306.60022 1.2285E+06
#1346.61853 66599
#1507.67017 248417
#1508.67334 247386
#1669.72925 139324
#1670.72620 69618.9
#END IONS

