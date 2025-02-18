

### Datasets in experiments:

```
[public real-world data]
	Voltage dataset: https://userscloud.com/rq69d3s865qc
	Price dataset: https://userscloud.com/zhhm350e5qn8
	Custom dataset: https://userscloud.com/0vtf66hjbhkx
[synthetic data]
	Zipf dataset: built with MainForZipfDataTxt.java
	Lognormal dataset: built with MainForLognormalDataTxt.java
	Continuous dataset: built with prepareUniform() in MainForKLL.java
```



### Out-of-database experiments:

Fig. 10 about verifying analyses on continuous duplicates

Fig. 10(a) about Average Error:

```
main.testError() in MainForKLL.java and MainForKLLDupliPair.java
```

Fig. 10(b) about Space Cost:

```
main.testM() in MainForKLL.java and MainForKLLDupliPair.java
```



Fig. 11,12,19,20 about verifying frequent items:

```
MainForFreqCountOfRealFrequency.java
	To obtain results on different datasets, change the parameter "dataType" in Line 236
```



Average Error in Fig. 13,14,15,16,17:

```
DD: main.testError() in MainForDD.java
DSS: main.testError() in MainForDSS.java
GK: main.testError() in MainForGK.java
TD: main.testError() in MainForTD.java
KLL: main.testError() in MainForKLL.java
Req: main.testError() in MainForReqSketch.java
DupliSketch: main.testError() in MainForKLLDupliPairFT.java

Select code in main() to change the parameters.
```



### In-database experiments:

Time cost in Fig. 13,14,15,16,17:

varyingN(), varyingMemory(), varyingQueriedQuantile(), varyingZipf(), varyingLognormal() in DupliQuantile_CompareBaselinesIT.java in the branch of IoTDB:

```
https://github.com/apache/iotdb/tree/research/dupli-quantile
```

