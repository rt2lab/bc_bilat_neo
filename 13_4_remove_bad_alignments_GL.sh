for i in {1..6}; do
full_name=`echo $line | cut -d$' ' -f 2`;
filename="${full_name##*/}" ; # get filename
dirname="${full_name%/*}"; # get directory/path name
echo $full_name 
full_name=./results/bc_bilat_ewok/bc_bilat/P$i\_GL/Variant_Calling/HC/P$i\_GL_uniq.nodup.onTarget.q20.realign.recal.hc.vcf
filename="${full_name##*/}" ; # get filename
dirname="${full_name%/*}"; # get directory/path name
echo $filename
echo $dirname 
echo `wc -l $full_name`;
new_filename=`echo $filename | sed s/\.vcf/_filter\.vcf/g`;
cp $full_name $dirname/$new_filename
while IFS= read -r pattern; do sed -i "/$pattern/d" $dirname/$new_filename; done < bilat_pattern_to_eliminate_bis.csv; done