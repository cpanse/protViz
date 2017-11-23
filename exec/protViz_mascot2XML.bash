#!/bin/bash -v

#$ -N mascot_RData_export
#$ -S /bin/bash
#$ -q PRX@fgcz-r-018

# Christian Panse <cp@fgcz.ethz.ch>
# Sun Apr  7 11:30:05 CEST 2013

# $HeadURL: http://fgcz-svn.uzh.ch/repos/fgcz/stable/bfabric/sgeworker/bin/fgcz_sge_mascot2RData $
# $Id: fgcz_sge_mascot2RData 8341 2017-05-25 21:29:55Z cpanse $
# $Date: 2017-05-25 23:29:55 +0200 (Thu, 25 May 2017) $

MASCOTCGI=/usr/local/mascot/mascot_2_5/cgi/

test -d MASCOTCGI || { echo "ERROR: '$MASCOTCGI' directory is  not available."; exit 1; }

function mascot_export() {
    DAT=$1
    XML=$2
    BASENAME=`basename $DAT .dat`
    
    [ -s $XML ] && { echo "file '$XML' is already there."; return; }
    touch $XML || { echo "can not create $XML file"; exit 1; }

    EXPORTOPTIONS="_minpeplen=5 _server_mudpit_switch=0.000000001 _showsubsets=1 _sigthreshold=0.05 do_export=1 export_format=XML group_family=1 pep_calc_mr=1 pep_delta=1 pep_end=1 pep_exp_mr=1 pep_exp_mz=1 pep_exp_z=1 pep_expect=1 pep_isbold=1 pep_isunique=1 pep_miss=1 pep_query=1 pep_rank=1 pep_scan_title=1 pep_score=1 pep_seq=1 pep_start=1 pep_var_mod=1 peptide_master=1 prot_acc=1 prot_cover=1 prot_desc=1 prot_empai=1 prot_hit_num=1 prot_len=1 prot_mass=1 prot_matches=1 prot_pi=1 prot_score=1 prot_seq=1 protein_master=1 query_master=1 query_params=1 query_peaks=1 query_qualifiers=1 query_raw=1 query_title=1 search_master=1 show_format=1 show_header=1 show_masses=1 show_mods=1 show_params=1 show_pep_dupes=1 use_homology=1 user=command line"

    test -s $DAT \
    || test -s $DAT.gz \
    && { zcat $DAT.gz > $SCRATCH/`basename $DAT .gz` \
        && DAT=$SCRATCH/`basename $DAT .gz`; \
        echo "unzip dat file $DAT."; }

    cd $MASCOTCGI \
    && ./export_dat_2.pl $EXPORTOPTIONS file=$DAT \
    > $XML

    [ $? -eq 0 ] \
    || { echo "error: export_dat_2.pl failed: $?"; exit 1; }
}

function mascot_xml2RData(){
    cd $1 \
&& R --no-save <<EOF
#R
library(XML)
for (fn in dir(pattern="F*.xml")){ 
    S <- XML::xmlToList(XML::xmlParse(fn))

    class(S) <-  c("list", "mascot")

    for (i in 1:length(S\$queries)){
        class(S\$queries[[i]]) <- c("list", "mascot_query")
    }                

    assign(gsub(".xml", "", fn), S)
}
remove(i)
remove(S)
remove(fn)
save.image(file='dump.RData', compress=TRUE)
EOF

    [ $? -eq 0 ] || { echo "xml parsing failed."; exit 1; } 
}

### ___MAIN___
test -f $1 \
  && mascot_export $1 $2 

exit 0
