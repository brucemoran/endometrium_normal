#################################################
#Process CHERs to overlap with features annotated
#################################################
#! /usr/bin/bash
tssb="/store2/bmoran/cardio/chris_meDip/Chris_Watson_5MeC_Array_data_posted/human_hg18_2010_11_10/hg18_gene-cpg-tss.bed"
miRlnc="/home/bmoran/bin/hg18/mirBase_lncipedia_3_1_hc-hg18.merge.bed"
rs="/store2/bmoran/cardio/chris_meDip/Chris_Watson_5MeC_Array_data_posted/ringo-starr_meDIP"

perl -ane 'if($F[0] eq "name"){next;}@sp=split(/\;/,$F[6]);$o=substr($F[-2],0,5);$t=substr($F[-1],0,5);print "$sp[1]\t$F[7]\t$F[8]\t$o;$t\n";' $rs"/_hypoxia_chers.txt" | sort -k1,1 -k2,2n > $rs"/hypo_chers.bed"

perl -ane 'if($F[0] eq "name"){next;}@sp=split(/\;/,$F[6]);$o=substr($F[-2],0,5);$t=substr($F[-1],0,5);print "$sp[1]\t$F[7]\t$F[8]\t$o;$t\n";' $rs"/_normoxia_chers.txt" | sort -k1,1 -k2,2n > $rs"/norm_chers.bed"

##annotate versus lncipedia, miRBase~/bin/hg18/~/bin/hg18/mirBase_lncipedia_3_1_hc-hg18.merge.bed
##looking for 'closest', so second line test if 2kb distance (should be - so upstream...)
bedtools closest -D "ref" -iu -a $rs"/hypo_chers.bed" -b $miRlnc | perl -ane 'if($F[-1]<2000){print $_;}'  > $rs"/reannotate/hypo-mirBase_lncipedia_3_1_hc-2kb.bed"

bedtools closest -D "ref" -iu -a $rs"/norm_chers.bed" -b $miRlnc | perl -ane 'if($F[-1]<2000){print $_;}'  > $rs"/reannotate/norm-mirBase_lncipedia_3_1_hc-2kb.bed"

bedtools closest -D "ref" -iu -a $rs"/hypo_chers.bed" -b $tssb | perl -ane 'if($F[-1]<2000){print $_;}'  > $rs"/reannotate/hypo-genes-cpg-tss_2kb.bed"

bedtools closest -D "ref" -iu -a $rs"/norm_chers.bed" -b $tssb | perl -ane 'if($F[-1]<2000){print $_;}'  > $rs"/reannotate/norm-genes-cpg-tss_2kb.bed"

##no limits!
bedtools closest -D "ref" -iu -a $rs"/hypo_chers.bed" -b $miRlnc > $rs"/reannotate/hypo-mirBase_lncipedia_3_1_hc.bed"

bedtools closest -D "ref" -iu -a $rs"/norm_chers.bed" -b $miRlnc > $rs"/reannotate/norm-mirBase_lncipedia_3_1_hc.bed"

bedtools closest -D "ref" -iu -a $rs"/hypo_chers.bed" -b $tssb  > $rs"/reannotate/hypo-genes-cpg-tss.bed"

bedtools closest -D "ref" -iu -a $rs"/norm_chers.bed" -b $tssb > $rs"/reannotate/norm-genes-cpg-tss.bed"

##headers
echo -e "CHR\tSTART \tEND\tMAX_LEVEL;SCORE\tANNO_CHR\tANNO_START\tANNO_END\tANNO_GENE-CpG\tREFSEQ_ID\tINT_IDS\tANNO_DISTANCE" > cpgenes.head
echo -e "CHR\tSTART \tEND\tMAX_LEVEL;SCORE\tANNO_CHR\tANNO_START\tANNO_END\tANNO_lncMir\tANNO_DISTANCE" > lncmir.head

for x in $rs/reannotate/*genes-cpg*bed;do
  outnm=$(echo $x | sed 's/bed/anno.bed/')
  cat cpgenes.head > $outnm;
  perl -ane 'chomp;
            %cpgene=();
            %nms=();
            %ids=();
            @s=split(/\;/,$F[7]);
            foreach $k (@s){
              if($k=~m/^ID=/){$k=~s/ID=//;$ids{$k}=$k;}
              if($k=~m/^Accession=/){$k=~s/Accession=//;$nms{$k}=$k;}
              if(($k=~m/^Name=/) || ($k=~m/^CpG:/)){$k=~s/Name=//;$cpgene{$k}=$k;}
            }
            @cpgenes=();@nmss=();@idss=();
            while(my($k,$v)=each(%ids)){
              push(@idss,$k);
            }
            while(my($k,$v)=each(%nms)){
              push(@nmss,$k);
            }
            while(my($k,$v)=each(%cpgene)){
              push(@cpgenes,$k);
            }
            print join("\t",@F[0..6]) . "\t" . join(";",@cpgenes[0..$#cpgenes]) . "\t" . join(";",@nmss[0..$#nmss]) . "\t" . join(";",@idss[0..$#idss]) . "\t$F[8]\n";' $x | \
            perl -ane 'chomp;if(($F[7]=~m/^CpG:/) && ($F[7]!~m/\;/)){
              print join("\t",@F[0..7]) . "\t-\t-\t$F[8]\n";}else{print $_ . "\n";
            }' >> $outnm
done

for x in $rs/reannotate/*lnc*bed;do
  outnm=$(echo $x | sed 's/bed/anno.bed/')
  cat lncmir.head > $outnm;
  cat $x >> $outnm;
done

##custom-perl-script alert!
perl ~/perl/textToExcel.pl "norm-genes-cpg-tss_2kb.anno.bed" "norm-mirBase_lncipedia_3_1_hc-2kb.anno.bed" "norm-genes-cpg-tss.anno.bed" "norm-mirBase_lncipedia_3_1_hc.anno.bed" "normoxia_hyper.meDIP.annotated"
perl ~/perl/textToExcel.pl "hypo-genes-cpg-tss_2kb.anno.bed" "hypo-mirBase_lncipedia_3_1_hc-2kb.anno.bed" "hypo-genes-cpg-tss.anno.bed" "hypo-mirBase_lncipedia_3_1_hc.anno.bed" "hypoxia_hyper.meDIP.annotated"
