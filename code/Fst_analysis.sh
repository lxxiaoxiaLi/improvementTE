#Fst
vcftools --vcf all.te.species.outgroup.SUPP.no247.vcf --weir-fst-pop im.txt --weir-fst-pop lan.txt --out im_lan/im_lan 
vcftools --vcf all.te.species.outgroup.SUPP.no247.vcf --weir-fst-pop Osi_im.txt --weir-fst-pop Osi_lan.txt --out Osi_im_lan/Osi_im_lan 
vcftools --vcf all.te.species.outgroup.SUPP.no247.vcf --weir-fst-pop Osj_im.txt --weir-fst-pop Osj_lan.txt --out Osj_im_lan/im_lan

#top 5% TE variations
awk -F'\t' '$3 > 0' Osi_im_lan.weir.fst|awk -F'\t' '$3 != "-nan"'|sort -k3,3rn|head -2193 >Osi_im_lan.top5.fst
awk -F'\t' '$3 > 0' im_lan.weir.fst|awk -F'\t' '$3 != "-nan"'|sort -k3,3rn|head -712 >Osj_im_lan.top5.fst
awk -F'\t' '$3 > 0' im_lan.weir.fst|awk -F'\t' '$3 != "-nan"'|sort -k3,3rn|head -2293 >im_lan.top5.fst
cat */*.top5.fst|cut -f2|sort -u >all.top5.fst 

#gene assiocated with top 5% TE variations 
for m in {Osi_im_lan,Osj_im_lan,im_lan};do bedtools intersect -a NIP.gff3 -b ${m}/${m}.top5.fst.bed -wb|grep LOC_Os|cut -f4|sed "s/ID=//g"|sed "s/_downstream2k//g"|sed "s/_gene//g"|sed "s/_upstream2k//g"|sed "s/_intro/\t/g"|sed "s/\./\t/g"|cut -f1|sort -u >${m}/${m}.top5.fst.gene;done;done


#top5 eGenes between im vs lan
awk  'NR==FNR{a[$0]}NR>FNR{ if(($1 in a)) print $0}' Osiimlan.allgene Osi_im_lan.top5.fst.gene >Osi_im_lan.top5.fst.siggene
awk  'NR==FNR{a[$0]}NR>FNR{ if(($1 in a)) print $0}' allimlan.allgene im_lan.top5.fst.gene >im_lan.top5.fst.siggene
awk  'NR==FNR{a[$0]}NR>FNR{ if(($1 in a)) print $0}'  Osjimlan.allgene Osj_im_lan.top5.fst.gene >Osj_im_lan.top5.fst.siggene
