

### Datasets in experiments:

```
[public real-world data]
	Voltage dataset: https://userscloud.com/rq69d3s865qc
	Price dataset: https://userscloud.com/zhhm350e5qn8
	Custom dataset: https://userscloud.com/0vtf66hjbhkx
[synthetic data]
	Zipf dataset: built with MainForZipfDataTxt.java
	Non-duplicate DHard dataset: built with MainForNonDupliDataTxt.java
```



### Out-of-database experiments:

Fig. 11(a):

```
main.testError() in MainForKLL.java and MainForKLLDupliPair.java
```

Fig. 11(b):

```
main.testM() in MainForKLL.java and MainForKLLDupliPair.java
```



Average Error in Fig. 12,13,14,15:

```
DD: main.testError() in MainForDD.java
DSS: main.testError() in MainForDSS.java
GK: main.testError() in MainForGK.java
TD: main.testError() in MainForTD.java
KLL: main.testError() in MainForKLL.java
KLL-Dupli: main.testError() in MainForKLLDupliPair.java

Select code in main() to vary parameters.
```



### In-database experiments:

Time cost in Fig. 12,13,14,15:

varyingN(),varyingMemory(),varyingQueriedQuantile(),varyingZipf() in DupliQuantile_CompareBaselinesIT.java in the branch of IoTDB:

```
https://github.com/apache/iotdb/tree/research/dupli-quantile
```

