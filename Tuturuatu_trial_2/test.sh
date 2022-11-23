for sample in ${datadir}*aug_R1.fq.gz
do
    base=$(basename $sample _R1.fq.gz)
    #base = Xxx_apr_aug
    name=$(echo $base | sed 's/_apr//g')
    echo $name
    rename "s/${base}/${name}/g" ${datadir}/${base}*
done




Trial relevant
        name=$(echo ${base} | sed 's/_S[0-9]/_S/g')
        name1=$(echo ${name} | sed 's/_S[0-9]//g')
        name2=$(echo ${name1} | sed 's/_S//g')