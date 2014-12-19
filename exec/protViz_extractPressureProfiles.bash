#!/bin/bash
# <cp@fgcz.ethz.ch>, Bernd Roschitzki <bernd.roschitzki@fgcz.uzh.ch>
# 2012-08-28 CP,BR
# 2012-09-06 CP,BR
# 2012-09-18 CP,CT

# USAGE
# goto the directory and just run the script by typping
# bash plotme_pressure_profile.bash

# some preprocessing
# we itereate over all sorted files ending with *.txt in the current dir.
#
# this is what we excpect as column names as input
# time(min)       Signal  Reference       Qa(nL/min)      Qb(nL/min)      Pc(psi) AUX_A/D


test -f data.txt && rm -iv data.txt

if [ ! -f data.txt ];
then
    find . -type f -name "*.txt" \
    |sort \
    | while read i;
    do
	sed 's//\n/g' $i \
        | awk -v fn=$i 'NF>6&&$1~/^[0-9]/{print fn"\t"$1"\t"$(NF-2)"\t"$(NF-1)"\t"$(NF-3)}' ;
    done \
    > data.txt
fi

# finnaly we have a file calles data.txt containg
# filename\tRT\t5thcolumn\t6th column



# now comes the fun
# R --no-save <<EOF
R --no-save <<EOF

# load data
s<-read.table("data.txt")


# open pdf device
pdf("pressure_profile_plot.pdf",19,12);

# two plot on one page
op<-par(mfrow=c(2,1))

# iterate over all files group files having similar rt range
dd<-round(tapply(s[,2], s[,1], max))

for (i in unique(sort(dd))){

        s.filter<-as.data.frame(s[s[,1] %in% names(dd[dd==i]),])
    s.filter.names<-as.character(unique(s.filter[,1]))

    cm<-rainbow(length(s.filter.names))

        plot(V5~V2,type="n",cex=0.1,
                data=s.filter,
		#ylim=c(0,500),
                xlab='time', ylab='Qb(nL/min)' , main="flow profile")

    for (j in c(1:length(s.filter.names))){
        lines(V3~V2,data=s.filter[s.filter[,1]==s.filter.names[j],], col=cm[j])
        lines(V5~V2,data=s.filter[s.filter[,1]==s.filter.names[j],], col=cm[j])
    }

    idx<-c(1:length(s.filter.names))
        legend("bottomleft", s.filter.names[idx], col=cm[idx],pch=22,cex=0.5)

###

        plot(V4~V2,type="n",cex=0.1,
                data=s.filter,
                xlab='time', ylab='Pc(psi)' , main="pressure profile")

    for (j in c(1:length(s.filter.names))){
        lines(V4~V2,data=s.filter[s.filter[,1]==s.filter.names[j],], col=cm[j])
    }

    idx<-c(1:length(s.filter.names))
        legend("bottomleft", s.filter.names[idx], col=cm[idx],pch=22,cex=0.5)
}
dev.off()
EOF


exit 0
