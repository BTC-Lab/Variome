### Steps for extracting annotated variants summary table:
#### Step 1:

```python
#!/bin/bash

status="status/status_step1.txt"

# filter files based on repeat maseker regions
for i in {1..22}; do
echo "chr$i start" >> $status
    file_name="maf_VEP_chr${i}.txt" # get file name
    #  varinats in  repeat maseker regions
    bcftools view -O v -o Table4_repeatmasker_byChr/Table4_repeatmasker_chr${i}_In.vcf --targets-file rmsk.bed  fully_annotated/maf_VEP_chr${i}.txt
    
    #  varinats not in  repeat maseker regions
    bcftools view -O v -o Table4_repeatmasker_byChr/Table4_repeatmasker_chr${i}_out.vcf --targets-file ^rmsk.bed  fully_annotated/maf_VEP_chr${i}.txt

    echo "$file_name is done"  >> $status
done 


# the same for chr X,Y,M
bcftools view -O v -o Table4_repeatmasker_byChr/Table4_repeatmasker_chrX_In.vcf --targets-file rmsk.bed fully_annotated/maf_VEP_chrX.txt # in repeatmasker
bcftools view -O v -o Table4_repeatmasker_byChr/Table4_repeatmasker_chrX_out.vcf --targets-file ^rmsk.bed  fully_annotated/maf_VEP_chrX.txt # not in repeatmasker

bcftools view -O v -o Table4_repeatmasker_byChr/Table4_repeatmasker_chrY_In.vcf --targets-file rmsk.bed  /egp-data/research-g42-khalifa-pfs/maf_g42/fully_annotated/maf_Illumina_VEP_chrY.txt # in repeatmaske
bcftools view -O v -o Table4_repeatmasker_byChr/Table4_repeatmasker_chrY_out.vcf --targets-file ^rmsk.bed  /egp-data/research-g42-khalifa-pfs/maf_g42/fully_annotated/maf_Illumina_VEP_chrY.txt # in repeatmasker

bcftools view -O v -o Table4_repeatmasker_byChr/Table4_repeatmasker_chrM_In.vcf --targets-file rmsk.bed fully_annotated/maf_VEP_chrM.txt # in repeatmasker
bcftools view -O v -o Table4_repeatmasker_byChr/Table4_repeatmasker_chrM_out.vcf --targets-file ^rmsk.bed fully_annotated/maf_VEP_chrM.txt # not in repeatmasker
echo "files are split by in and out of repeat masker regions!"

###################################################################################################################################################################################################################
###################################################################################################################################################################################################################
###################################################################################################################################################################################################################

# check for novel varaints
for i in {1..22}; do
    start_time=$(date +"%H:%M:%S")
    echo "chr$i start time: $start_time" >> "$status"
    file_name="maf_VEP_chr${i}.txt"
    
    #  NOVEL varinats not in  repeat maseker regions
    filter_vep -i Table4_repeatmasker_byChr/Table4_repeatmasker_chr${i}_out.vcf -o Table4_repeatmasker_byChr/Table4_repeatmasker_novel/Table4_repeatmasker_chr${i}_not_noRS.vcf --filter "not Existing_variation"

    #  NOVEL varinats in  repeat maseker regions
    filter_vep -i Table4_repeatmasker_byChr/Table4_repeatmasker_chr${i}_In.vcf -o Table4_repeatmasker_byChr/Table4_repeatmasker_novel/Table4_repeatmasker_chr${i}_In_noRS.vcf --filter "not Existing_variation"
    echo "$file_name is done" 
    end_time=$(date +"%H:%M:%S")
    echo "chr$i end time: $end_time" >> "$status"
done 

# the same for chr X,Y,M
grep -E Table4_repeatmasker_byChr/Table4_repeatmasker_chrM_out.vcf | filter_vep -o Table4_repeatmasker_byChr/Table4_repeatmasker_novel/Table4_repeatmasker_chrM_not_noRS.vcf -filter "not Existing_variation"
grep -E Table4_repeatmasker_byChr/Table4_repeatmasker_chrM_In.vcf | filter_vep -o Table4_repeatmasker_byChr/Table4_repeatmasker_novel/Table4_repeatmasker_chrM_In_noRS.vcf -filter "not Existing_variation"
echo "M done!"
grep -E Table4_repeatmasker_byChr/Table4_repeatmasker_chrX_out.vcf | filter_vep -o Table4_repeatmasker_byChr/Table4_repeatmasker_novel/Table4_repeatmasker_chrX_not_noRS.vcf -filter "not Existing_variation"
grep -E Table4_repeatmasker_byChr/Table4_repeatmasker_chrX_In.vcf | filter_vep -o Table4_repeatmasker_byChr/Table4_repeatmasker_novel/Table4_repeatmasker_chrX_In_noRS.vcf -filter "not Existing_variation"
echo "X done!"
grep -E Table4_repeatmasker_byChr/Table4_repeatmasker_chrY_out.vcf | filter_vep -o Table4_repeatmasker_byChr/Table4_repeatmasker_novel/Table4_repeatmasker_chrY_not_noRS.vcf -filter "not Existing_variation"
grep -E Table4_repeatmasker_byChr/Table4_repeatmasker_chrY_In.vcf | filter_vep -o Table4_repeatmasker_byChr/Table4_repeatmasker_novel/Table4_repeatmasker_chrY_In_noRS.vcf -filter "not Existing_variation"
echo "Y done!"

echo "files are split by NOVEL in and out of repeat masker regions!"


# check for known varaints
for i in {1..22}; do
    start_time=$(TZ="Asia/Dubai" date +"%H:%M:%S")
    echo "chr$i Known start time: $start_time" >> "$status"

    file_name="maf_VEP_chr${i}.txt"
    #  KNOWN varinats not in  repeat maseker regions
    filter_vep -i ../Table4_repeatmasker_byChr/Table4_repeatmasker_chr${i}_out.vcf -o ../Table4_repeatmasker_byChr/Table4_repeatmasker_Known/Table4_repeatmasker_chr${i}_not_RS.vcf --filter "Existing_variation"
    endy_time=$(TZ="Asia/Dubai" date +"%H:%M:%S")
    echo "chr$i done not in RR: $endy_time" >> "$status"

    #  KNWON varinats in  repeat maseker regions
    filter_vep -i ../Table4_repeatmasker_byChr/Table4_repeatmasker_chr${i}_In.vcf -o ../Table4_repeatmasker_byChr/Table4_repeatmasker_Known/Table4_repeatmasker_chr${i}_In_RS.vcf --filter "Existing_variation"
    echo "$file_name is done" 
    end_time=$(TZ="Asia/Dubai" date +"%H:%M:%S")
    echo "chr$i Known end time: $end_time" >> "$status"
done 

# the same for chr X,Y,M
grep -E ../Table4_repeatmasker_byChr/Table4_repeatmasker_chrM_out.vcf | filter_vep -o ../Table4_repeatmasker_byChr/Table4_repeatmasker_Known/Table4_repeatmasker_chrM_not_RS.vcf -filter "Existing_variation"
grep -E ../Table4_repeatmasker_byChr/Table4_repeatmasker_chrM_In.vcf | filter_vep -o ../Table4_repeatmasker_byChr/Table4_repeatmasker_Known/Table4_repeatmasker_chrM_In_RS.vcf -filter "Existing_variation"
echo "M done!"
grep -E ../Table4_repeatmasker_byChr/Table4_repeatmasker_chrX_out.vcf | filter_vep -o ../Table4_repeatmasker_byChr/Table4_repeatmasker_Known/Table4_repeatmasker_chrX_not_RS.vcf -filter "Existing_variation"../
grep -E ../Table4_repeatmasker_byChr/Table4_repeatmasker_chrX_In.vcf | filter_vep -o Table4_repeatmasker_byChr/Table4_repeatmasker_Known/Table4_repeatmasker_chrX_In_RS.vcf -filter "Existing_variation"
echo "X done!"
grep -E ../Table4_repeatmasker_byChr/Table4_repeatmasker_chrY_out.vcf | filter_vep -o ../Table4_repeatmasker_byChr/Table4_repeatmasker_Known/Table4_repeatmasker_chrY_not_RS.vcf -filter "Existing_variation"
grep -E ../Table4_repeatmasker_byChr/Table4_repeatmasker_chrY_In.vcf | filter_vep -o ../Table4_repeatmasker_byChr/Table4_repeatmasker_Known/Table4_repeatmasker_chrY_In_RS.vcf -filter "Existing_variation"
echo "Y done!"

echo "files are split by KNOWN in and out of repeat masker regions!"

```

#### Step 2:

```python
#!/bin/bash

# Create a CSV file to store the results
output_file="output_files/step1_allVariants_noFiler.csv"
status="status/status_step2.txt"

# Array of keywords
IFS=$'\n' keywords=($(cat "impacts_table_1less.txt"))
IFS=$'\n' keywords_all=($(cat "impacts_table.txt"))

header="transcript_ablation"
for keyword in "${keywords[@]}"
do
    header="$header,$keyword"
done
echo "$header" >> "$output_file"

val=("X" "Y" "M")
# Loop through each file
for i in {1..22}
do
    file="fully_annotated/maf_VEP_chr${i}.txt"
    file_chr="${file:76}"
    start_time=$(date +"%H:%M:%S")
    echo "$file_chr start time: $start_time" >> "$status"
    
    # Initialize an array to store counts
    counts=("$file_chr")  
    counts1=()
    # Loop through the keywords and count occurrences in the file
    
    for keyword in "${keywords_all[@]}"
    do
        count=("$(grep -o -i -w -c "$keyword" "$file")")
        counts+=("$count")
        counts1=($(echo "${counts[*]}" | tr '\n' ','))
    done
    counts1=${counts1::-1}
    
    # Append the results to the CSV file
    echo "${counts1[*]}" >> "$output_file"
    end_time=$(date +"%H:%M:%S")
    echo "$file_chr end time: $end_time" >> "$status"
done
```


#### Step 3:

```python
#!/bin/bash

# Create a CSV file to store the results
output_file="output_files/step3_novelVariants_notRepeat.csv"
status="status/status_step3.txt"

# Array of keywords
IFS=$'\n' keywords=($(cat "impacts_table_1less.txt"))
IFS=$'\n' keywords_all=($(cat "impacts_table.txt"))


header="transcript_ablation"
for keyword in "${keywords[@]}"
do
    header="$header,$keyword"
done
echo "$header" >> "$output_file"


# Loop through each file
for i in {1..22}
do
    file="../Table4_repeatmasker_byChr/Table4_repeatmasker_novel/Table4_repeatmasker_chr${i}_not_noRS.txt"
    file_chr="${file:75}"
    echo "start ${file_chr}"
    start_time=$(date +"%H:%M:%S")
    echo "$file_chr start time: $start_time" >> "$status"
    # Initialize an array to store counts
    counts=("$file_chr")  
    counts1=()
    # Loop through the keywords and count occurrences in the file
    for keyword in "${keywords_all[@]}"
    do
        count=("$(grep -o -i -w -c "$keyword" "$file")")
        counts+=("$count")
        counts1=($(echo "${counts[*]}" | tr '\n' ','))
    done
    counts1=${counts1::-1}
    # Append the results to the CSV file
    #echo "${counts1[*]}" >> "$output_file"
    echo "${counts1[*]}" >> "$output_file"
    end_time=$(date +"%H:%M:%S")
    echo "$file_chr end time: $end_time" >> "$status"
done
```


#### Step 4:

```python
#!/bin/bash

# Create a CSV file to store the results
output_file="output_files/step4_novelVariants_inRepeat.csv"
status="status/status_step4.txt"
# Array of keywords
IFS=$'\n' keywords=($(cat "impacts_table_1less.txt"))
IFS=$'\n' keywords_all=($(cat "impacts_table.txt"))


header="transcript_ablation"
for keyword in "${keywords[@]}"
do
    header="$header,$keyword"
done
echo "$header" >> "$output_file"

# Loop through each file
for file in Table4_repeatmasker_byChr/Table4_repeatmasker_novel/Table4_repeatmasker_chr*_in_noRS.txt
do
    file_chr="${file:81}"
    echo "start ${file_chr}"
    start_time=$(date +"%H:%M:%S")+4
    echo "$file_chr start time: $start_time" >> "$status"
    
    # Initialize an array to store counts
    counts=("$file_chr")  
    counts1=()
    # Loop through the keywords and count occurrences in the file
    for keyword in "${keywords_all[@]}"
    do
        count=("$(grep -o -i -w -c "$keyword" "$file")")
        counts+=("$count")
        counts1=($(echo "${counts[*]}" | tr '\n' ','))
    done
    counts1=${counts1::-1}
    # Append the results to the CSV file
    #echo "${counts1[*]}" >> "$output_file"
    echo "${counts1[*]}" >> "$output_file"
    end_time=$(date +"%H:%M:%S")+4
    echo "$file_chr end time: $end_time" >> "$status"
done
```


#### Step 5:

```python
#!/bin/bash

# Create a CSV file to store the results
output_file="../output_files/step6_KnownInRepeat.csv"
status="status/status_step5.txt"

# Array of keywords
IFS=$'\n' keywords=($(cat "impacts_table_1less.txt"))
IFS=$'\n' keywords_all=($(cat "impacts_table.txt"))


header="transcript_ablation"
for keyword in "${keywords[@]}"
do
    header="$header,$keyword"
done
echo "$header" >> "$output_file"

# Loop through each file
for file in Table4_repeatmasker_byChr/Table4_repeatmasker_chr*_out.vcf
do
    file_chr="${file:46}"
    echo "start ${file_chr}"
    start_time=$(TZ="Asia/Dubai" date +"%H:%M:%S")
    
    echo "$file_chr start time: $start_time" >> "$status"
    
    # Initialize an array to store counts
    counts=("$file_chr")  
    counts1=()
    # Loop through the keywords and count occurrences in the file
    for keyword in "${keywords_all[@]}"
    do
        count=("$(grep -o -i -w -c "$keyword" "$file")")
        counts+=("$count")
        counts1=($(echo "${counts[*]}" | tr '\n' ','))
    done
    counts1=${counts1::-1}
    # Append the results to the CSV file
    #echo "${counts1[*]}" >> "$output_file"
    echo "${counts1[*]}" >> "$output_file"
    end_time=$(TZ="Asia/Dubai" date +"%H:%M:%S")
    echo "$file_chr end time: $end_time" >> "$status"
done
```


#### Step 6:

```python
#!/bin/bash

# Create a CSV file to store the results
output_file="../output_files/step6_notRepeat_Known.csv"
status="status/status_step6.txt"

# Array of keywords
IFS=$'\n' keywords=($(cat "impacts_table_1less.txt"))
IFS=$'\n' keywords_all=($(cat "impacts_table.txt"))


header="transcript_ablation"
for keyword in "${keywords[@]}"
do
    header="$header,$keyword"
done
echo "$header" >> "$output_file"


# Loop through each file
for file in Table4_repeatmasker_byChr/Table4_repeatmasker_novel/Table4_repeatmasker_chr*_not_RS.vcf
do
    file_chr="${file:49}"
    echo "start ${file_chr}"
    start_time=$(TZ="Asia/Dubai" date +"%H:%M:%S")
    
    echo "$file_chr start time: $start_time" >> "$status"
    
    # Initialize an array to store counts
    counts=("$file_chr")  
    counts1=()
    # Loop through the keywords and count occurrences in the file
    for keyword in "${keywords_all[@]}"
    do
        count=("$(grep -o -i -w -c "$keyword" "$file")")
        counts+=("$count")
        counts1=($(echo "${counts[*]}" | tr '\n' ','))
    done
    counts1=${counts1::-1}
    # Append the results to the CSV file
    #echo "${counts1[*]}" >> "$output_file"
    echo "${counts1[*]}" >> "$output_file"
    end_time=$(TZ="Asia/Dubai" date +"%H:%M:%S")
    echo "$file_chr end time: $end_time" >> "$status"
done
```


#### Step 7:

```python
#!/bin/bash

# Create a CSV file to store the results
output_file="../output_files/step7_inRepeat_Known.csv"
status="status/status_step7.txt"

# Array of keywords
IFS=$'\n' keywords=($(cat "impacts_table_1less.txt"))
IFS=$'\n' keywords_all=($(cat "impacts_table.txt"))


header="transcript_ablation"
for keyword in "${keywords[@]}"
do
    header="$header,$keyword"
done
echo "$header" >> "$output_file"

val=("20" "21" "22" "X" "Y" "M")
# Loop through each file
for file in Table4_repeatmasker_byChr/Table4_repeatmasker_novel/Table4_repeatmasker_chr*_In_RS.vcf
do
    file_chr="${file:49}"
    echo "start ${file_chr}"
    start_time=$(TZ="Asia/Dubai" date +"%H:%M:%S")
    
    echo "$file_chr start time: $start_time" >> "$status"
    
    # Initialize an array to store counts
    counts=("$file_chr")  
    counts1=()
    # Loop through the keywords and count occurrences in the file
    for keyword in "${keywords_all[@]}"
    do
        count=("$(grep -o -i -w -c "$keyword" "$file")")
        counts+=("$count")
        counts1=($(echo "${counts[*]}" | tr '\n' ','))
    done
    counts1=${counts1::-1}
    # Append the results to the CSV file
    echo "${counts1[*]}" >> "$output_file"
    end_time=$(TZ="Asia/Dubai" date +"%H:%M:%S")
    echo "$file_chr end time: $end_time" >> "$status"
done
```

