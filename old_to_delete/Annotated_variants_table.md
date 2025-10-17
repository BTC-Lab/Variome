#### This code is dependent on Annotated_variants_shell_script.md

```python
import os
import glob
import pandas as pd
import numpy as np
pd.set_option('display.max_rows', 200)
```


```python
# function to capitalize first letter only
def cap_first(input):
    return input[0].upper() + input[1:]
```


```python
# step 2 All variants
all_varinats = pd.read_csv("../output_files/step1_allVariants_noFiler.csv", delimiter=',')
all_varinats = all_varinats.drop(['chr12.txt_warnings.txt', 'chr17.txt_warnings.txt', 'transcript_ablation', 'chrrange21.txt', 'chrrange22.txt'])
all_varinats = all_varinats.astype(int)
all_varinats.head(2)
```






```python
#step 3 Novel not in repeat masker region
novel_notRepeatmasker = pd.read_csv("../output_files/step3_novelVariants_notRepeat.csv", delimiter=',', index_col=0)
novel_notRepeatmasker
novel_notRepeatmasker.head(2)
```





```python
#step 4 Novel in repeat masker region
novel_inRepeatmasker = pd.read_csv("../output_files/step4_novelVariants_inRepeat.csv", delimiter=',')
novel_inRepeatmasker
novel_inRepeatmasker.head(2)
```



```python
# all variants in repeat
all_inRepeatmasker = pd.read_csv("../output_files/step5_inRepeat.txt", delimiter=',')
all_inRepeatmasker.head(2)
```



```python
# all variants in repeat
notRepeatmasker_sub = pd.read_csv("../output_files/step5_notRepeat.csv", delimiter=',', names = list(novel_inRepeatmasker.columns))
notRepeatmasker_sub.head(1)
```



```python
#step 5 All variants in repeat masker region
inRepeatmasker = all_inRepeatmasker[~all_inRepeatmasker.index.str.endswith("_out.vcf")]
inRepeatmasker.head(1)#.sum(axis=0)
```



```python
# All variants not in repeat masker region
notRepeatmasker = all_inRepeatmasker[all_inRepeatmasker.index.str.endswith("_out.vcf")]
notRepeatmasker.head(4)#.sum(axis=0)
```




```python
# Known variants  in repeat masker region
inRepeatmasker_Known = pd.read_csv("../output_files/step6_notRepeat_Known.csv", delimiter=',')
inRepeatmasker_Known = inRepeatmasker_Known.drop(['transcript_ablation'])
inRepeatmasker_Known.head(4)
```




```python
# Known variants  in repeat masker region
inRepeatmasker_Known = pd.read_csv("../output_files/step7_inRepeat_Known.csv", delimiter=',')
inRepeatmasker_Known.head(4)
```





```python
imapact_sev= pd.read_csv("Impact_sev.txt", delimiter='\t', header=None, index_col=0)
imapact_sev =imapact_sev.rename_axis('impact',index=True)
imapact_sev= imapact_sev.sort_index(ascending=True)
imapact_sev = imapact_sev.rename(columns={1:'severity'})
imapact_sev.index = imapact_sev.index.to_series().apply(cap_first).str.replace('_', ' ' )
```


```python
# table of combined sums of the number of variants
combinedsums = pd.DataFrame({'All': all_varinats.sum(axis=0),'In repeat regions':inRepeatmasker.sum(axis=0), 'In repeat regions "Novel"': novel_inRepeatmasker.sum(axis=0), 'Not in repeat regions "Novel"':novel_notRepeatmasker.sum(axis=0)}) 
combinedsums.index = combinedsums.index.to_series().apply(cap_first).str.replace('_', ' ' )
combinedsums['Not in repeat regions "Novel"'] = combinedsums['Not in repeat regions "Novel"'].astype("int")
combinedsums['Not in repeat regions'] = combinedsums['All'] - combinedsums['In repeat regions']
combinedsums =combinedsums.rename_axis('impact',index=True)
combinedsums = combinedsums.sort_index(ascending=True)
```


```python
Table4 = combinedsums.join(imapact_sev)
Table4['Not in repeat regions'] = Table4['All'] - Table4['In repeat regions']
impact_order =  ["HIGH", "MODERATE", "LOW", "MODIFIER"]
Table4['severity'] = pd.Categorical(Table4['severity'], categories=impact_order)
Table4=Table4.sort_values(by=['severity', 'All', 'In repeat regions','In repeat regions "Novel"','Not in repeat regions "Novel"'] , ascending=[True, False, False, False, False])
Table4
```


```python
Table4_cleaned = combinedsums[combinedsums.values.sum(axis=1)!=0]

Table4_cleaned = Table4_cleaned.join(imapact_sev)
impact_order =  ["HIGH", "MODERATE", "LOW", "MODIFIER"]
Table4_cleaned['severity'] = pd.Categorical(Table4_cleaned['severity'], categories=impact_order)
Table4_cleaned = Table4_cleaned.sort_values(by=['severity', 'All', 'In repeat regions','In repeat regions "Novel"','Not in repeat regions "Novel"'] , ascending=[True, False, False, False, False])
```


```python
z = Table4_cleaned['In repeat regions']/Table4_cleaned['All']*100
Table4_cleaned.insert(2, 'In RR', z.round(0))

z = Table4_cleaned['In repeat regions "Novel"']/Table4_cleaned['All']*100
Table4_cleaned.insert(4, 'In RR Novel', z.round(0))

z = Table4_cleaned['Not in repeat regions "Novel"']/Table4_cleaned['All']*100
Table4_cleaned.insert(6, 'Not RR Novel', z.round(0))

z = Table4_cleaned['Not in repeat regions']/Table4_cleaned['All']*100
Table4_cleaned.insert(8, 'Not RR', z.round(0))
```


```python
Table4_cleaned
```


