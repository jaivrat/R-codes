mkdir -p data
for file in `cat datafiles.txt`
do
cd data
wget $file
cd ../
done
