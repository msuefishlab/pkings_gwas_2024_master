 var="CHR SNP BP A1 A2 TEST AFF UNAFF P"

awk 'NR <=1 || $5 =="REC" {print $0}' $1 >  $1.rec
awk '{split($2,a,"_"); print a[1],$2,a[2],$3,$4,$5,$6,$7,$8}' $1.rec > $1.rec.gwas
sed -i "1s/.*/$var/" $1.rec.gwas

awk 'NR <=1 || $5 =="GENO" {print $0}' $1 >  $1.geno
awk '{split($2,a,"_"); print a[1],$2,a[2],$3,$4,$5,$6,$7,$8}' $1.geno > $1.geno.gwas
sed -i "1s/.*/$var/" $1.geno.gwas

awk 'NR <=1 || $5 =="DOM" {print $0}' $1 >  $1.dom
awk '{split($2,a,"_"); print a[1],$2,a[2],$3,$4,$5,$6,$7,$8}' $1.dom > $1.dom.gwas
sed -i "1s/.*/$var/" $1.dom.gwas

awk 'NR <=1 || $5 =="ALLELIC" {print $0}' $1 >  $1.allelic
awk '{split($2,a,"_"); print a[1],$2,a[2],$3,$4,$5,$6,$7,$8}' $1.allelic > $1.allelic.gwas
sed -i "1s/.*/$var/" $1.allelic.gwas
