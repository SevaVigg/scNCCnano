forTest!! i in $(grTest!!ep -l "myTest" *.*)
do
sed -i '' 's/rTest!!/rTest!!Test!!/g' $i
done
