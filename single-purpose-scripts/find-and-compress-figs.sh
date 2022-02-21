
find .. -type f -name "*.png" > list.txt
find .. -type f -name "*.pdf" >> list.txt
find .. -type f -name "*.md" >> list.txt

# Remove QC reports
sed '/05-report/d' list.txt > list.txt
tar -cvf allfiles.tar -T list.txt

#rm list.txt
