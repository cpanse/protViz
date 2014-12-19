#!/bin/bash
# $HeadURL: http://fgcz-svn.unizh.ch/repos/fgcz/testing/proteomics/R/protViz/exec/protViz_fun.bash $
# $Id: protViz_fun.bash 6089 2014-02-11 08:02:53Z cpanse $


find ../ -type f \
| egrep "Rd|Rnw" \
| grep -v svn \
| while read l;
do
    cat $l;
done \
| tr -cs "[:alpha:]" "\n" \
| awk '{h[$1]++}END{for (k in h){print k" "h[k]}}'  \
> /tmp/word.freq.txt


R --no-save <<EOF
library(tm)
library(wordcloud)
s<-read.table("/tmp/word.freq.txt",sep=" ")
s.filter<-s[!s[,1] %in% stopwords(c("en","SMART")),]
pdf("/tmp/wordcloud-ProtViz.pdf", 10,10)
wordcloud(s.filter[,1],s.filter[,2])
dev.off()
EOF
