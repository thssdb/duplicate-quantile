import it.unimi.dsi.fastutil.HashCommon;
import it.unimi.dsi.fastutil.ints.IntArrayList;
import it.unimi.dsi.fastutil.longs.Long2LongMap;
import it.unimi.dsi.fastutil.longs.Long2LongOpenHashMap;
import it.unimi.dsi.fastutil.objects.ObjectArrayList;
import it.unimi.dsi.util.XoRoShiRo128PlusRandom;
import org.apache.commons.math3.distribution.NormalDistribution;

import java.util.Arrays;
import java.util.Comparator;

public class KLLDupliFT extends KLLSketchForQuantile {
  IntArrayList compactionNumInLevel;
  public XoRoShiRo128PlusRandom randomForReserve = new XoRoShiRo128PlusRandom();
  long MIN_V=Long.MAX_VALUE,MAX_V=Long.MIN_VALUE;
  int[] levelPosT; // range in level L:   [  levelPos[L]  ,  lvT[L]  ),  [  lvT[L]  ,  levelPos[L+1]  )
  Long2LongOpenHashMap freqValueCount=null;
  boolean freqValueCountFrozen=false,frozenCheckingHashMapForNonFreqValue=false;
  long[] freqValue,freqCountSum;
  int maxMemoryByte,maxMemoryNumForSketch,lastUpdateFreqLevel=0;
  long lastCheckFreqValueN=0,newSketchNAfterLastCheck=0,lastThreshold=0;
  final int BitsPerItem=64;
  int AmortizedRatioForBuffer=50,MaxFreqThresholdRatio=10;

  public KLLDupliFT(int maxMemoryByte) {
    N = 0;
    this.maxMemoryByte=maxMemoryByte;
    calcParameters(maxMemoryByte);
    calcLevelMaxSize(1);
  }

  private void calcParameters(int maxMemoryByte) {
    maxMemoryNumForSketch=calcMaxMemoryNum(maxMemoryByte);
    num = new long[maxMemoryNumForSketch];
    level0Sorted = false;
    cntLevel = 0;
    compactionNumInLevel = new IntArrayList();
  }

  @Override
  protected int calcMaxMemoryNum(int maxMemoryByte) {
    return Math.min(1 << 20, maxMemoryByte*8 / BitsPerItem);
  }

  @Override
  protected void calcLevelMaxSize(int setLevel) { // set cntLevel.  make sure cntLevel won't decrease
    int[] tmpArr = new int[setLevel + 1];
    int[] tmpArrT = new int[setLevel + 1];
    int maxPos = cntLevel > 0 ? Math.max(maxMemoryNumForSketch, levelPos[cntLevel]) : maxMemoryNumForSketch;
    for (int i = 0; i < setLevel + 1; i++) {
      tmpArr[i] = i < cntLevel ? levelPos[i] : maxPos;
      tmpArrT[i] = i < cntLevel ? levelPosT[i] : maxPos;
    }
    levelPos = tmpArr;
    levelPosT = tmpArrT;
    for(int i=cntLevel;i<setLevel;i++) {
      compactionNumInLevel.add(0);
    }
    cntLevel = setLevel;
    levelMaxSize=calcLevelMaxSizeByLevel(maxMemoryNumForSketch,cntLevel);
//    System.out.println("\t\tmaxMemoryNumForSketch:"+maxMemoryNumForSketch+"\tmaxPos:"+maxPos+"\tnum.len:"+num.length+"\tN:"+N+"\tsetLV:"+setLevel+"\tlvmaxsize:"+ Arrays.toString(levelMaxSize));
  }



  public String toString() {
    sortLV0();
    final StringBuilder sb = new StringBuilder();
    sb.append(N);
    sb.append(levelMaxSize[cntLevel - 1]);
    sb.append(cntLevel);
    sb.append((levelPos[cntLevel] - levelPos[0]));
    for (int i = 0; i < cntLevel; i++) sb.append(levelPos[i]);
    return sb.toString();
  }


  public void showLevelStat() {
    int numLEN = levelPos[cntLevel] - levelPos[0];
    System.out.println("\t\t//maxMemoryNumForSketch:" + maxMemoryNumForSketch + "\t//N:" + N);
    for (int i = 0; i < cntLevel; i++)
      System.out.print("\t\t" + levelMaxSize[i] + "\t");
    System.out.println();
    System.out.println("-------------------------------------------------------");
  }
  public void show(){
    sortLV0();
    for(int i=0;i<cntLevel;i++){
      System.out.print("\t");
      System.out.print("["+(levelPos[i+1]-levelPos[i])+"]"+"("+(levelPos[i+1]-levelPosT[i])+")");
      System.out.print("\t");
    }System.out.println("\tmaxLV="+(cntLevel-1));
  }
  public void showDoubleNum(){
    sortLV0();
    for(int i=0;i<cntLevel;i++){
      System.out.print("\t|");
      for(int j=levelPos[i];j<levelPos[i+1];j++)
//        System.out.print(num[j]+",");
        System.out.print(longToResult(num[j])+"("+(j>=levelPosT[i])+")"+", ");
      System.out.print("|\t");
    }System.out.println();
  }
  public void showNum(){
    sortLV0();
    for(int i=0;i<cntLevel;i++){
      System.out.print("\t|");
      for(int j=levelPos[i];j<levelPos[i+1];j++)
//        System.out.print(num[j]+",");
        System.out.print((num[j])+"("+(j>=levelPosT[i]?"T":"F")+")"+", ");
      System.out.print("|\t");
    }System.out.println();
  }

  @Override
  public void update(long x) { // signed long
    freqValueCountFrozen=false;
//    checkN();

    if(freqValueCount!=null&&freqValueCount.containsKey(x)) {
      freqValueCount.addTo(x, 1);
      N++;
    }
    else {
      while (levelPos[0] <= 0) {
        while(levelPos[0]<=Math.max(MinimumLevelMaxSize,getLevelSize(0)/AmortizedRatioForBuffer)) {
          compact();
          sumForAvgDupli += checkDupliElementsInSketch();
          countForAvgDupli += 1;
        }
        if(needToCheckFreqValue())
          checkFreqValue();
      }
      --levelPos[0];
      --levelPosT[0];
      num[levelPos[0]] = x;
      level0Sorted = false;
      N++;
      newSketchNAfterLastCheck++;
    }
//    showNum();
  }

  protected void randomlyHalveDownToLeft(int len,long[] tmpNum,boolean[] tmpIsT){
    int delta = getNextRand01();
    int mid = len>>>1;
    for(int i=0,j=0;i<mid;i++,j+=2) {
      tmpIsT[i] = tmpNum[j]==tmpNum[j+1]&&tmpIsT[j]&&tmpIsT[j+1];
      tmpNum[i] = tmpNum[j + delta];
    }
  }

  public int PairT0T0 =0, PairTF =0, PairFF =0, PairT1T2 =0;
  protected void mergeSortWithoutSpace(int level, int L1,int mid,long[] tmpNum, boolean[] tmpIsT,int L2,int PT2, int R2){
    int p1 = L1, p2=L2, cntPos=mid;
    int dupliPair=0;
    while(true){
      while(p1<mid&&tmpIsT[p1-L1])p1++;
      if(p1>=mid&&p2>=PT2)break;
      if(p1<mid&&(p2==PT2||tmpNum[p1-L1]<num[p2])){
        num[cntPos]=tmpNum[p1-L1];
        p1++;
      }
      else {
        num[cntPos]=num[p2++];
      }
      cntPos++;
    }
    levelPosT[level+1] = cntPos;
    p1=L1;
    p2=PT2;
    while(true){
      while(p1<mid&&!tmpIsT[p1-L1])p1++;
      if(p1>=mid&&p2>=PT2)break;
      if(p1<mid&&(p2==R2||tmpNum[p1-L1]<num[p2])){
        num[cntPos]=tmpNum[p1-L1];
        p1++;
      }
      else {
        num[cntPos]=num[p2++];
      }
      cntPos++;
    }
    for(int i=levelPosT[level+1]+1;i<R2;i++)
      if(num[i]==num[i-1]){
        dupliPair++;
        i++;
      }
  }
  private long getMinVInSortedLevel(int level){
    long minV=num[levelPos[level]];
    if(levelPosT[level]<levelPos[level+1])minV=Math.min(minV,num[levelPosT[level]]);
    return minV;
  }
  private long getMaxVInSortedLevel(int level){
    long maxV=num[levelPos[level+1]-1];
    if(levelPosT[level]>levelPos[level])maxV=Math.max(maxV,num[levelPosT[level]-1]);
    return maxV;
  }

  private void compactOneLevel(int level) {
    if (level == cntLevel - 1)
      calcLevelMaxSize(cntLevel + 1);
    int L1 = levelPos[level], PT1=levelPosT[level], R1 = levelPos[level + 1]; // [L,R)
    if (level == 0 && !level0Sorted)
      sortLV0();
    if (R1-L1<=1) return;
    int reservedF=0;
    if(((R1-L1)&1)==1){
      reservedF=PT1>L1?1:0;
      L1++;
      PT1=Math.max(PT1,L1);
    }
//    System.out.println("\t\t\t???L1,PT1,R1:\t"+L1+","+PT1+","+R1+"\t\tnum.len:"+num.length+"\t\tN:"+N+"\tlevel:"+level);
//    checkN();
    boolean[] tmpIsT = new boolean[R1-L1];
    long[] tmpNum=new long[R1-L1];
    int tmpCnt=0,tmpP1=L1,tmpP2=PT1;
    for(int i=0;i<R1-L1;i++){
      if(tmpP1<PT1&&(tmpP2>=levelPos[level+1]||num[tmpP1]<=num[tmpP2])){
        tmpNum[tmpCnt]=num[tmpP1++];
        tmpIsT[tmpCnt]=false;
      }else{
        tmpNum[tmpCnt]=num[tmpP2++];
        tmpIsT[tmpCnt]=true;
      }
      tmpCnt++;
    }
    int tmpPairT0T0=0,tmpPairFF=0,tmpPairTF=0,tmpPairT1T2=0;
    for(int i=1;i<tmpCnt;i+=2){
      if(tmpIsT[i-1]&&tmpIsT[i]){
        if(tmpNum[i-1]==tmpNum[i]) tmpPairT0T0++;
        else tmpPairT1T2++;
      }else if(!tmpIsT[i-1]&&!tmpIsT[i]) tmpPairFF++;
      else tmpPairTF++;
    }
    PairT0T0+=tmpPairT0T0;
    PairFF+=tmpPairFF;
    PairTF+=tmpPairTF;
    PairT1T2+=tmpPairT1T2;
    int toMergeT0T0Pairs =
//        1<<20;
        2;
//        tmpCnt/2/10+1; // TODO DEBUG    选择无损压缩的阈值
    int toMergePairs =
//        tmpCnt/2;
        (levelMaxSize[level]+1)/2; // TODO DEBUG    有损时是否全部压缩掉
    boolean
        toMergeFF=tmpPairT0T0<toMergeT0T0Pairs,
        toMergeTF=toMergeFF&&tmpPairT0T0+tmpPairFF<toMergePairs,
        toMergeT1T2=toMergeFF&&tmpPairT0T0+tmpPairFF+tmpPairTF<toMergePairs;
//    if(!toMergeFF)
//      System.out.println("\t\tisToMergeTypes: "+toMergeFF+"\t"+toMergeTF+"\t"+toMergeT1T2);
    boolean[] isToMergeNum = new boolean[tmpCnt];
    if(toMergeFF)
      addRecord(false, getMinVInSortedLevel(level), getMaxVInSortedLevel(level), level);
    for(int i=1;i<tmpCnt;i+=2){
      if(tmpIsT[i-1]&&tmpIsT[i]){
        if(tmpNum[i-1]==tmpNum[i])isToMergeNum[i-1]=isToMergeNum[i]=true;
        else isToMergeNum[i-1]=isToMergeNum[i]=toMergeT1T2;
      }else if(!tmpIsT[i-1]&&!tmpIsT[i]) isToMergeNum[i-1]=isToMergeNum[i]=toMergeFF;
      else  isToMergeNum[i-1]=isToMergeNum[i]=toMergeTF;
    }
    for(int i=0;i<tmpCnt;i++)if(!tmpIsT[i]&&!isToMergeNum[i]){
      num[L1]=tmpNum[i];
      reservedF++;
      L1++;
    }
    for(int i=0;i<tmpCnt;i++)if(tmpIsT[i]&&!isToMergeNum[i]){
      num[L1]=tmpNum[i];
      L1++;
    }
    int toMergeTmpCnt=0;
    for(int i=0;i<tmpCnt;i++)if(isToMergeNum[i]){
      tmpNum[toMergeTmpCnt]=tmpNum[i];
      tmpIsT[toMergeTmpCnt]=tmpIsT[i];
      toMergeTmpCnt++;
    }
//    if(toMergeTmpCnt!=tmpCnt){
//      System.out.println("\t\t\tDEBUG\ttoMergeTmpCnt!=tmpCnt\t"+toMergeTmpCnt+"\t"+tmpCnt);
//    }

//    System.out.print("compactOneFullLevel:LV="+level+"\t F:");
//    for(int i=L1;i<PT1;i++)System.out.print(" "+num[i]);
//    System.out.print("\t T:");
//    for(int i=PT1;i<levelPos[level+1];i++)System.out.print(" "+num[i]);
//    System.out.print("\n\t\t\ttmpNum:\t");
//    for(int i=0;i<R1-L1;i++)System.out.print(" "+tmpNum[i]+"("+(tmpIsT[i]?"T":"F")+")");
//    System.out.print("\n");
    randomlyHalveDownToLeft(toMergeTmpCnt,tmpNum,tmpIsT);
    mergeSortWithoutSpace(level, L1,L1+toMergeTmpCnt/2, tmpNum, tmpIsT, levelPos[level + 1], levelPosT[level+1], levelPos[level + 2]);
    levelPos[level + 1] = L1+toMergeTmpCnt/2;
    int newP = levelPos[level + 1] - 1, oldP = L1 - 1;
    for (int i = oldP; i >= levelPos[0]; i--) {
      num[newP] = num[oldP];
      newP--;oldP--;
    }
    levelPos[level] = levelPos[level + 1] - (L1 - levelPos[level]);
    levelPosT[level] = levelPos[level]+reservedF;
//    if(reservedF)System.out.println("\t\t\treservedF:\t"+num[levelPos[level]]+"\t");
    int numReduced = (R1 - L1) >>> 1;
    for (int i = level - 1; i >= 0; i--) {
      levelPos[i] += numReduced;
      levelPosT[i] +=numReduced;
    }
//    this.checkN();
  }

  private long getFreqMinimalThreshold(){
//    return (8L<<(cntLevel-1));
    return (4L<<(cntLevel-1));
  }
  private boolean needToCheckFreqValue(){
    return lastUpdateFreqLevel<=cntLevel&&newSketchNAfterLastCheck>=getFreqMinimalThreshold()/2 && !frozenCheckingHashMapForNonFreqValue;
  }
  private void checkFreqValue() {
    lastUpdateFreqLevel = cntLevel;
    lastCheckFreqValueN = N;
    newSketchNAfterLastCheck=0;
    Long2LongOpenHashMap fullFreqValueCount = freqValueCount!=null?new Long2LongOpenHashMap(freqValueCount):new Long2LongOpenHashMap();
    long cntValue, cntCount;
    for (int lv = 0; lv < cntLevel; lv++)
      for (int i = levelPosT[lv], pair = 0; i < levelPos[lv+1]; i++) {
        cntValue = num[i];
        cntCount = 1;
        fullFreqValueCount.addTo(cntValue, cntCount << lv);
//        System.out.println("\t\t\t\tV:"+longToResult(cntValue)+"\tC:"+(cntCount << lv));
      }

    int freqThresholdRatio = 1;
    Long2LongOpenHashMap newFreqValueCount = null;
    for (; freqThresholdRatio <= MaxFreqThresholdRatio; freqThresholdRatio++) {
      if(freqThresholdRatio*getFreqMinimalThreshold()<lastThreshold)continue;
      Long2LongOpenHashMap tmpFreqValueCount = new Long2LongOpenHashMap();
      for (Long2LongMap.Entry entry : fullFreqValueCount.long2LongEntrySet())
        if (entry.getLongValue() >= freqThresholdRatio*getFreqMinimalThreshold())
          tmpFreqValueCount.put(entry.getLongKey(), entry.getLongValue());
      if (tmpFreqValueCount.isEmpty()) {
//        compact();
//        System.out.println("\t[checkFreqValue] Empty NewFreq. no-op. exit.");
        return;
        // TODO: if old values do not meet the threshold, should clear freqValueCount.
      }
      int tmpMaxMemoryNumForSketch = (maxMemoryByte - 2 * 8 * HashCommon.arraySize(tmpFreqValueCount.size(), Long2LongOpenHashMap.DEFAULT_LOAD_FACTOR)) * 8 / BitsPerItem;
      int tmpSketchNumLen = 0;
      for (int lv = 0; lv < cntLevel; lv++) {
        tmpSketchNumLen += levelPosT[lv] - levelPos[lv];
        for (int i = levelPosT[lv]; i < levelPos[lv+1]; i++) {
          cntValue = num[i];
          if (!tmpFreqValueCount.containsKey(cntValue)) {
            tmpSketchNumLen++;
          }
        }
      }
      // todo: add up tmpSketchNumLen for some non-freq in hashmap (will be moved to lv0 pairs).
//      System.out.println("\t\ttry threshold:"+freqThresholdRatio*getFreqMinimalThreshold()+"\tnewNumLen:"+tmpSketchNumLen+"\tshould not exceed "+tmpMaxMemoryNumForSketch);
      if (tmpSketchNumLen < tmpMaxMemoryNumForSketch) {
        newFreqValueCount = tmpFreqValueCount;
        break;
      }
    }
    if(freqThresholdRatio>MaxFreqThresholdRatio){
//      System.out.println("\t[checkFreqValue] NewFreq makes newNum too large. no-op. exit.");
//      System.out.print("\t[checkFreqValue] Sketch:");showNum();
      return;
    }
    maxMemoryNumForSketch = (maxMemoryByte - 2 * 8 * HashCommon.arraySize(newFreqValueCount.size(), Long2LongOpenHashMap.DEFAULT_LOAD_FACTOR)) * 8 / BitsPerItem;
    levelMaxSize = calcLevelMaxSizeByLevel(maxMemoryNumForSketch, cntLevel);
    lastThreshold=freqThresholdRatio*getFreqMinimalThreshold();
//    System.out.println("\t[checkFreqValue]\tthreshold:\t"+freqThresholdRatio*getFreqMinimalThreshold()+"\tmaxMemoryNumForSketch:"+maxMemoryNumForSketch+"\tN:"+N);
//    System.out.print("\t[checkFreqValue] Sketch:\t");showNum();
    long[] newNum = new long[maxMemoryNumForSketch];
    int[] newLevelPos = new int[cntLevel + 1], newLvPosT = new int[cntLevel + 1];


    int cntPos = newLevelPos[0] = 0;
    for (int lv = 0; lv < cntLevel; lv++) {
      for (int i = levelPos[lv]; i < levelPosT[lv]; i++) newNum[cntPos++] = num[i];
      newLvPosT[lv] = cntPos;
      for (int i = levelPosT[lv]; i < levelPos[lv+1]; i++) {
        cntValue = num[i];
        if (!newFreqValueCount.containsKey(cntValue)) {
//          System.out.println("\t\t\t!nonFreq cntV,C:"+longToResult(cntValue)+","+cntCount);
          newNum[cntPos] = cntValue;
          cntPos++;
        }
      }
//      System.out.println("\t\t\t\t???space for newPairT:"+((tmpNewPair + 1) / 2));
      newLevelPos[lv + 1] = cntPos;
    }
    int delta = maxMemoryNumForSketch - cntPos;
    for (int i = maxMemoryNumForSketch-1; i >= delta; i--) newNum[i] = newNum[i - delta];
    for (int i = 0; i <= cntLevel; i++) {
      newLevelPos[i] += delta;
      newLvPosT[i] += delta;
    }
    num = newNum;
    levelPos = newLevelPos;
    levelPosT = newLvPosT;
    Long2LongOpenHashMap oldMap = freqValueCount; // todo: simply move to lv0 pairs.
    freqValueCount = newFreqValueCount;
    if (oldMap != null) {
      Long2LongOpenHashMap toRemoveMap = new Long2LongOpenHashMap();
      long mmp=0,mmpDistinct=0;
//      System.out.print("PRE checkN.\t");
//      checkN();
      for (Long2LongMap.Entry entry : oldMap.long2LongEntrySet())
        if (!freqValueCount.containsKey(entry.getLongKey())) {
          N -= entry.getLongValue();
          mmp+=entry.getLongValue();
          mmpDistinct++;
          toRemoveMap.put(entry.getLongKey(), entry.getLongValue());
        }
//      System.out.print("MID checkN.\t");
//      checkN();
//      if(mmp>0)System.out.println("MID toRemoveN:"+mmp+"\t\tavgCount:"+mmp/mmpDistinct);
      frozenCheckingHashMapForNonFreqValue=true;
      for (Long2LongMap.Entry entry : toRemoveMap.long2LongEntrySet()) {
//        System.out.println("\t\t!!!WARN\t delete non-freq item:\t" + longToResult(entry.getLongKey()) + "," + entry.getLongValue());
        for (long j = entry.getLongValue(); j > 0; j--)
          update(entry.getLongKey());
        mmp-=entry.getLongValue();
      }
      frozenCheckingHashMapForNonFreqValue=false;
//      System.out.print("FINAL checkN.\t");
//      checkN();
    }
//        checkN();
  }


  @Override
  public void compact() {
    sortLV0();
    for (int i = 0; i < cntLevel; i++)
      if (getLevelSize(i) > levelMaxSize[i]||i==cntLevel-1) {
        compactOneLevel(i);
        break;
      }
  }

//  @Override
//  public void compact() {
//    sortLV0();
//    for (int i = 0; i < cntLevel-1; i++)
//      if (getLevelSize(i) > levelMaxSize[i]) {
////        showNum();
////        System.out.println("\t\tcompact:\tlv:"+i+"\tdupli:"+dupliPairInLevel.getInt(i) * 2+"/"+getLevelSize(i)+"="+(1.0*dupliPairInLevel.getInt(i) * 2 / getLevelSize(i)));
////        if (dupliPairInLevel.getInt(i) * 2 >= getLevelSize(i) / 2)
//////          compactOneFullLevel(i);
////          compactDupliInLevel(i);
////        else
//          compactOneLevel(i);
////          compactNonDupliInLevel(i);
////        System.out.println("\t\tafterLVSize:"+getLevelSize(i));
////        showNum();
////        System.out.println("\n");
//        return;
//      }
//
//
//    if(needToCheckFreqValue()) {
//      lastUpdateFreqLevel = cntLevel;
//      lastCheckFreqValueN = N;
////      for (int i = 0; i < cntLevel - 1; i++)
////        if (dupliPairInLevel.getInt(i) >= 1)
////          compactDupliInLevel(i);
//      for (int threshold = 2; threshold <= 4; threshold++) {
//        LongArrayList freqValue = new LongArrayList(), freqCount = new LongArrayList();
//        for (int i = levelPosT[cntLevel - 1]; i < levelPos[cntLevel]; i++){
//            int j = i + 1;
//            while (j < levelPos[cntLevel] && num[j] == num[i]) j++;
//            if (j - i >= threshold) {
//              freqValue.add(num[i]);
//              freqCount.add(0);
//            }
//            i = j - 1;
//          }
////      for(int i=1;i<freqValue.size();i++)if(freqValue.getLong(i)<=freqValue.getLong(i-1))System.out.println("\t\t????freqValue??\t"+freqValue.getLong(i-1)+" "+freqValue.getLong(i));
////      System.out.println("\t\tfreqValuesInTopLevel:\t"+freqValue);
//        for (int lv = 0; lv <= cntLevel - 1; lv++) {
//          int tmpPos = 0;
//          for (int i = levelPosT[lv]; i < levelPos[lv + 1]; i++) {
//            while (tmpPos < freqValue.size() && freqValue.getLong(tmpPos) < num[i]) tmpPos++;
//            if (tmpPos >= freqValue.size()) break;
//            if (freqValue.getLong(tmpPos) != num[i]) continue;
//            freqCount.set(tmpPos, freqCount.getLong(tmpPos) + (1L << lv));
//          }
//        }
//        long tmpFreqN = 0;
//        for (long count : freqCount) tmpFreqN += count;
////      System.out.println("\t\t\t\t\ttmpFreqN:"+tmpFreqN);
//
//        Long2LongOpenHashMap newFreqValueCount = new Long2LongOpenHashMap();
//        for (int i = 0; i < freqValue.size(); i++)
//          newFreqValueCount.put(freqValue.getLong(i), freqCount.getLong(i));
////      System.out.println("\t\tcheck freq value distinct:\t\t"+freqValue.size()+"\t"+newFreqValueCount.size());
//        if (freqValueCount != null) {
//          for (Long2LongMap.Entry entry : freqValueCount.long2LongEntrySet())
//            if (!newFreqValueCount.containsKey(entry.getLongKey())) {
//              if (entry.getLongValue() >= ((long) threshold << (cntLevel - 1)))
//                newFreqValueCount.put(entry.getLongKey(), entry.getLongValue());
//            } else
//              newFreqValueCount.addTo(entry.getLongKey(), entry.getLongValue());
//        }
////      System.out.println("\t\t\thashMap:" + 2 * 8 * HashCommon.arraySize(newFreqValueCount.size(), Long2LongOpenHashMap.DEFAULT_LOAD_FACTOR) + "Byte");
//
//        if (newFreqValueCount.isEmpty()) {
//          compact();
//          return;
//          // TODO  check non-freq items in hashmap
//        }
//
//        int tmpRestNum = (maxMemoryByte - 2 * 8 * HashCommon.arraySize(newFreqValueCount.size(), Long2LongOpenHashMap.DEFAULT_LOAD_FACTOR)) * 8 / BitsPerItem;
//
//
//        int cntPos = tmpRestNum - 1;
//        for (int lv = cntLevel - 1; lv >= 0; lv--)
//          for (int i = levelPos[lv + 1] - 1; i >= levelPos[lv]; i--)
//            if (i<levelPosT[lv] || !newFreqValueCount.containsKey(num[i])) {
//              cntPos--;
//            }
//        if(cntPos<0){
//          continue;
//        }
//
//        maxMemoryNumForSketch=tmpRestNum;
//        levelMaxSize = calcLevelMaxSizeByLevel(maxMemoryNumForSketch, cntLevel);
////        System.out.println("\t\t\tthreshold:\t"+threshold+"\tmaxMemoryNumForSketch:"+maxMemoryNumForSketch+"\tN:"+N+"\tleP:"+ Arrays.toString(levelPos));
//        long[] newNum = new long[maxMemoryNumForSketch];
//        int[] newLevelPos = new int[cntLevel + 1];
//        int[] newLevelPosT = new int[cntLevel + 1];
//
//        cntPos = maxMemoryNumForSketch - 1;
//        newLevelPos[cntLevel] = maxMemoryNumForSketch;
//        for (int lv = cntLevel - 1; lv >= 0; lv--) {
//          for (int i = levelPos[lv + 1] - 1; i >= levelPosT[lv]; i--)
//            if (!newFreqValueCount.containsKey(num[i])) {
//              newNum[cntPos] = num[i];
//              cntPos--;
//            }
//          newLevelPosT[lv] = cntPos+1;
//          for (int i = levelPosT[lv] - 1; i >= levelPos[lv]; i--){
//              newNum[cntPos] = num[i];
//              cntPos--;
//            }
//          newLevelPos[lv] = cntPos + 1;
//        }
//        levelPos = newLevelPos;
//        levelPosT=newLevelPosT;
//        num = newNum;
//        Long2LongOpenHashMap oldMap = freqValueCount;
//        freqValueCount = newFreqValueCount;
//        if (oldMap != null)
//          for (Long2LongMap.Entry entry : oldMap.long2LongEntrySet())
//            if (!freqValueCount.containsKey(entry.getLongKey()))
//              N -= entry.getLongValue();
////        checkN();
//        if (oldMap != null)
//          for (Long2LongMap.Entry entry : oldMap.long2LongEntrySet())
//            if (!freqValueCount.containsKey(entry.getLongKey())) {
////          System.out.println("\t\t!!!WARN\t delete non-freq item:\t");
//              for (long j = entry.getLongValue(); j > 0; j--)
//                update(entry.getLongKey());
//            }
////        checkN();
//        if (levelPos[0] > 0) return;
//        else break;
//      }
//    }
//    calcLevelMaxSize(cntLevel+1);
//    compactOneLevel(cntLevel - 2);
////    this.showLevelStat();
////    this.showNum();
//  }

  public int getMaxErr(){
    int totERR=0;
    for(int i=0;i<cntLevel;i++)
      totERR+=compactionNumInLevel.getInt(i)<<i;
    return totERR;
  }

  public void showCompact(){
    int totERR=0;
    long totSIGMA=0;
    for(int i=0;i<cntLevel;i++){
      System.out.print("\t");
      System.out.print("["+compactionNumInLevel.getInt(i)+"]");
      System.out.print("\t");
      totERR+=compactionNumInLevel.getInt(i)<<i;
      totSIGMA+= (long) compactionNumInLevel.getInt(i) <<(i*2L);
    }System.out.println("\tmaxLV="+(cntLevel-1)+"\ttotERR="+totERR+"\ttotSIGMA="+totSIGMA);
  }

  private void addRecordInLevel(int level) {
    int hisNum = compactionNumInLevel.getInt(level);
      hisNum++;
    compactionNumInLevel.set(level, hisNum);
  }

  public void addRecord(boolean isSum, long minV, long maxV, int level) {
    MIN_V=Math.min(MIN_V,minV);
    MAX_V=Math.max(MAX_V,maxV);
    if (!isSum) {
      addRecordInLevel(level);
    }else{
      for(int i=0;i<=level;i++)
        addRecordInLevel(i);
    }
  }
  static long lastBound;
  static double lastSig2,lastPr=-1;
  static NormalDistribution lastNormalDis;

  public int[] getRelatedCompactNum() {
    int[] relatedCompactNum = new int[cntLevel - 1];
    for (int i = 0; i < cntLevel - 1; i++)
        relatedCompactNum[i] = compactionNumInLevel.getInt(i);
    return relatedCompactNum;
  }

  public void sortLV0(){// dupli=true in lv0, just sort num[]
    if(level0Sorted)return;
    Arrays.sort(num, levelPos[0], levelPos[1]);
    level0Sorted=true;

    boolean hasLastNum=false;
    int dupliPair=0;
    for(int i=levelPos[0];i<levelPos[1];i++)
      if(hasLastNum&&num[i-1]==num[i]) {
        hasLastNum = false;
        dupliPair++;
      }
    else hasLastNum=true;
  }


  public static long queryRankErrBoundGivenParameter(double sig2, long maxERR, double Pr) {
    if(maxERR==0)return 0;
    NormalDistribution dis = sig2==lastSig2&&lastNormalDis!=null?lastNormalDis:new NormalDistribution(0, Math.sqrt(sig2));
    if(sig2==lastSig2&&Pr==lastPr)return lastBound;
    double tmpBound = dis.inverseCumulativeProbability(0.5*(1-Pr));
    tmpBound=Math.min(-Math.floor(tmpBound),maxERR);
//    System.out.println("\t\t\t\t..\t\ttmp\t\t"+sig2+"\t"+maxERR+"\t"+Pr+"\t\t\ttmpBound:\t"+tmpBound);
    lastBound=(long)tmpBound;
    lastPr=Pr;
    lastSig2=sig2;
    lastNormalDis=dis;
//    System.out.println("\tPr:"+Pr+"\t"+"sig2:\t"+sig2+"\t\tbyInvCum:\t"+(long)tmpBound);
    return (long)tmpBound;
  }

  public static long queryRankErrBound(int[] relatedCompactNum, double Pr) {
    double sig2=0;
    long maxERR=0;
    for(int i=0;i<relatedCompactNum.length;i++)sig2+=0.5*relatedCompactNum[i]*Math.pow(2,i*2);
    for(int i=0;i<relatedCompactNum.length;i++)maxERR+=(long)relatedCompactNum[i]<<i;
//    System.out.println("\t\t\t\t\t\tqueryRankErrBound\tmaxERR:\t"+maxERR+"\t\tsig2:\t"+ sig2+"\t\t\n"+ Arrays.toString(relatedCompactNum));
    return queryRankErrBoundGivenParameter(sig2,maxERR,Pr);
  }

  public long queryRankErrBound(double Pr){
    return queryRankErrBound(compactionNumInLevel.toIntArray(),Pr);
  }

  // 返回误差err，该数值在sketch里估计排名的偏差绝对值有Pr的概率<err
  public long queryRankErrBound(long result, double Pr) {
    int[] relatedCompactNum = getRelatedCompactNum();
    return queryRankErrBound(relatedCompactNum,Pr);
  }


  protected int findRankInLevel(int level,long v){
    if(levelPos[level]>=levelPos[level+1])return 0;
    if(level==0&&!level0Sorted)
      sortLV0();
    long ansF=0,ansT=0;
    if(levelPos[level]>=levelPosT[level]||num[levelPos[level]]>v)ansF=0;
    else {
      int L = levelPos[level],R = levelPosT[level]-1;
      while (L < R) {
        int mid = (L + R + 1) >> 1;
        if (num[mid] <= v) L = mid;
        else R = mid - 1;
      }
      ansF = (L - levelPos[level] + 1) * (1 << level);
    }
    if(levelPosT[level]>=levelPos[level+1]||num[levelPosT[level]]>v)ansT=0;
    else {
      int L = levelPosT[level],R = levelPos[level+1]-1;
      while (L < R) {
        int mid = (L + R + 1) >> 1;
        if (num[mid] <= v) L = mid;
        else R = mid - 1;
      }
      ansT = (L - levelPosT[level] + 1) * (1 << level);
    }
    return (int)(ansF+ansT);
  }
  @Override
  public int getApproxRank(long v){
    int approxRank = 0;
    for(int i=0;i<cntLevel;i++) if(levelPos[i]<levelPos[i+1]){
      approxRank+=findRankInLevel(i,v);
//      for (int j = levelPos[i]; j < levelPos[i + 1]; j++)
//        if (num[j] < v) approxRank += 1 << i;
    }
    if(freqValueCount!=null) {
      if(!freqValueCountFrozen){
        freqValue=new long[freqValueCount.size()];
        freqCountSum=new long[freqValueCount.size()];
        ObjectArrayList<Long2LongMap.Entry> entries = new ObjectArrayList<>(freqValueCount.size());
        entries.addAll(freqValueCount.long2LongEntrySet());
        entries.sort(Comparator.comparingLong(Long2LongMap.Entry::getLongKey));
        int i=0;long tmpSum=0;
        for (Long2LongMap.Entry entry : entries){
          freqValue[i]=entry.getLongKey();
          tmpSum+=entry.getLongValue();
          freqCountSum[i]=tmpSum;
          i++;
        }
        freqValueCountFrozen=true;
      }
      int pos=-1;
      for(int k=Integer.highestOneBit(freqValue.length);k>0;k/=2){
        if(pos+k<freqValue.length&&freqValue[pos+k]<=v)pos+=k;
      }
      if(pos>=0)approxRank+=freqCountSum[pos];
//      for (Long2LongMap.Entry entry : freqValueCount.long2LongEntrySet())
//        if (entry.getLongKey() <= v)
//          approxRank += entry.getLongValue();
//      for (Long2LongMap.Entry entry : freqValueCount.long2LongEntrySet())
//        if (entry.getLongKey() <= v)
//          approxRank += entry.getLongValue();
    }
    return approxRank;
  }
  public long findMinValueWithRank(long K){
    long L=Long.MIN_VALUE,R=Long.MAX_VALUE,mid;
    while(L<R){
      mid = L + ((R - L) >>>1);
//      System.out.println("2fen\t\t"+L+"..."+R+"\t\tmid="+mid+"\t\t"+(getApproxRank(mid)>=K));
      if(getApproxRank(mid)>K)R=mid;
      else L=mid+1;
    }
//    System.out.println("FT K:"+K+"\tN:"+getN()+" rank(L):"+getApproxRank(L));
    return L;
  }

  public void checkN(){
    long numN=0;
    for(int i=0;i<cntLevel;i++)numN+=getLevelSize(i)*(1L<<i);//System.out.println("\t\t\tN_in_sketch:"+numN);
    if(freqValueCount!=null)
    for (Long2LongMap.Entry entry : freqValueCount.long2LongEntrySet())
      numN += entry.getLongValue();
    if(numN!=N)
      System.out.println("\t\t[ERR]checkN="+numN+"\tN="+N);
  }

  public double sumForAvgDupli=0,countForAvgDupli=0;
  public double checkDupliElementsInSketch(){
    if(levelPos[0]==levelPos[cntLevel])return 1.0;
    sortLV0();
    int tmpSize=0,tmpDistinct=0;
    for(int lv=0;lv<cntLevel;lv++) {
      int lvSize = getLevelSize(lv), lvDistinct=0;
      for (int i = levelPos[lv]; i < levelPosT[lv]; i++) {
        if(i==levelPos[lv]||num[i]!=num[i-1])lvDistinct++;
      }
      for (int i = levelPosT[lv]; i < levelPos[lv+1]; i++) {
        if(i==levelPosT[lv]||num[i]!=num[i-1])lvDistinct++;
      }
      tmpSize+=lvSize;
      tmpDistinct+=lvDistinct;
    }
//    System.out.println("\t\t\t[DEBUG]\tSize/Distinct:\t"+1.0*tmpSize/tmpDistinct);
    return 1.0*tmpSize/tmpDistinct;
  }
  public double getAvgDupliInSketch(){
    if(countForAvgDupli>0)return sumForAvgDupli/countForAvgDupli;
    else return checkDupliElementsInSketch();
  }


}