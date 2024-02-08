import it.unimi.dsi.fastutil.HashCommon;
import it.unimi.dsi.fastutil.ints.IntArrayList;
import it.unimi.dsi.fastutil.longs.Long2LongMap;
import it.unimi.dsi.fastutil.longs.Long2LongOpenHashMap;
import it.unimi.dsi.fastutil.longs.LongArrayList;
import it.unimi.dsi.fastutil.objects.ObjectArrayList;
import it.unimi.dsi.util.XoRoShiRo128PlusRandom;
import org.apache.commons.math3.distribution.NormalDistribution;

import java.util.Arrays;
import java.util.Comparator;

public class KLLDupliWBitsSwitch extends KLLSketchForQuantile {
  static public boolean onlyDealFullyDupli=true;
  IntArrayList compactionNumInLevel,dupliPairInLevel;
  public XoRoShiRo128PlusRandom randomForReserve = new XoRoShiRo128PlusRandom();
  long MIN_V=Long.MAX_VALUE,MAX_V=Long.MIN_VALUE;
  boolean[] dupli;
  boolean[] dupliNum;
  Long2LongOpenHashMap freqValueCount=null;
  int maxMemoryByte,maxMemoryNumForSketch,lastUpdateFreqLevel=0;
  int wBits;
  int BitsPerItem;

  public KLLDupliWBitsSwitch(int maxMemoryByte, int wBits) {
    N = 0;
    this.maxMemoryByte=maxMemoryByte;
    this.wBits=wBits;
    BitsPerItem=64;//+wBits/*+1*/;
    calcParameters(maxMemoryByte);
    calcLevelMaxSize(1);
  }

  private void calcParameters(int maxMemoryByte) {
    maxMemoryNumForSketch=maxMemoryNum = calcMaxMemoryNum(maxMemoryByte);
    num = new long[maxMemoryNum];
    dupli = new boolean[maxMemoryNum];
    dupliNum = new boolean[maxMemoryNum*wBits];
    level0Sorted = false;
    cntLevel = 0;
    compactionNumInLevel = new IntArrayList();
    dupliPairInLevel = new IntArrayList();
  }

  private long getDupliNum(int level,int p){
    long dup=0;
    for(int i=p*wBits;i<p*wBits+wBits;i++)dup=dup<<1|(dupliNum[i]?1:0);
    if(level>=wBits)dup<<=(level+1-wBits);
    return dup+1;
  }
  private long mergeDupNum(long dupNum1,long dupNum2){
    return (dupNum1+dupNum2)>>>1;
  }
  private void setDupliNum(int level, int p,long dupNum){
    dupNum-=1;
    dupNum>>>=Math.max(0,level+1-wBits);
    for(int i=0;i<wBits;i++)dupliNum[p*wBits+i] = (dupNum&(1L<<(wBits-i-1)))>0;
  }

  @Override
  protected int calcMaxMemoryNum(int maxMemoryByte) {
    return Math.min(1 << 20, maxMemoryByte*8 / BitsPerItem);
  }

  @Override
  protected void calcLevelMaxSize(int setLevel) { // set cntLevel.  make sure cntLevel won't decrease
    int[] tmpArr = new int[setLevel + 1];
    int maxPos = cntLevel > 0 ? Math.max(maxMemoryNumForSketch, levelPos[cntLevel]) : maxMemoryNumForSketch;
    for (int i = 0; i < setLevel + 1; i++) tmpArr[i] = i < cntLevel ? levelPos[i] : maxPos;
    levelPos = tmpArr;
    for(int i=cntLevel;i<setLevel;i++) {
      compactionNumInLevel.add(0);
      dupliPairInLevel.add(0);
    }
    cntLevel = setLevel;
    levelMaxSize=calcLevelMaxSizeByLevel(maxMemoryNumForSketch,cntLevel);
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
    System.out.println("\t\t//maxMemNum:" + maxMemoryNum + "\t//N:" + N);
    for (int i = 0; i < cntLevel; i++)
      System.out.print("\t\t" + levelMaxSize[i] + "\t");
    System.out.println();
    System.out.println("-------------------------------------------------------");
  }
  public void show(){
    sortLV0();
    for(int i=0;i<cntLevel;i++){
      System.out.print("\t");
      System.out.print("["+(levelPos[i+1]-levelPos[i])+"]"+"("+dupliPairInLevel.getInt(i)+")");
      System.out.print("\t");
    }System.out.println("\tmaxLV="+(cntLevel-1));
  }
  public void showNum(){
    sortLV0();
    for(int i=0;i<cntLevel;i++){
      System.out.print("\t|");
      for(int j=levelPos[i];j<levelPos[i+1];j++)
//        System.out.print(num[j]+",");
        System.out.print(longToResult(num[j])+"("+dupli[j]+")"+"("+getDupliNum(i,j)+")"+", ");
      System.out.print("|\t");
    }System.out.println();
  }

  @Override
  public void update(long x) { // signed long
    if (levelPos[0] == 0) compact();
    N++;
    if(freqValueCount!=null&&freqValueCount.containsKey(x))
      freqValueCount.addTo(x,1);
    else {
      num[--levelPos[0]] = x;
      dupli[levelPos[0]] = true;
      setDupliNum(0,levelPos[0], 1);
      level0Sorted = false;
    }
  }


  protected void randomlyHalveDownToLeft(int level, int L,int R){
    int delta = getNextRand01();
    int mid = (L+R)>>>1;
    for(int i=L,j=L;i<mid;i++,j+=2) {
      num[i] = num[j + delta];
      if(num[j]==num[j+1]){
        dupli[i] = dupli[j]&&dupli[j+1];
        setDupliNum(level+1,i,mergeDupNum(getDupliNum(level,j),getDupliNum(level,j+1)));
      }else{
        dupli[i]=false;
        setDupliNum(level+1,i,getDupliNum(level,j+delta));
      }


    }
    for(int i=L;i<mid;){
      int j=i,count1=0,pos1=i;
      while(j<mid&&num[j]==num[i]){
        if(dupli[j]){
          dupli[j]=dupli[pos1];
          dupli[i]=true;
          long tmp=getDupliNum(level+1,j);
          setDupliNum(level+1,j,getDupliNum(level+1,i));
          setDupliNum(level+1,i,tmp);
          pos1++;
        }
        j++;
      }
//      for(int k=i;k<i+count1;k++)dupli[k]=true;
//      for(int k=i+count1;k<j;k++)dupli[k]=false;
      i=j;
    }
  }
  protected void mergeSortWithoutSpace(int level, int L1,int mid,int L2,int R2){
    int p1 = L1, p2=L2, cntPos=mid;
    boolean hasLastNum=false;
    int dupliPair=0;
    while(p1<mid||p2<R2){
      if(p1<mid&&(p2==R2||num[p1]<num[p2]||(num[p1]==num[p2]&&dupli[p1]))){
        num[cntPos]=num[p1];
        dupli[cntPos]=dupli[p1];
        for(int j=0;j<wBits;j++)dupliNum[cntPos*wBits+j]=dupliNum[p1*wBits+j];
        p1++;
      }
      else {
        num[cntPos]=num[p2];
        dupli[cntPos]=dupli[p2];
        for(int j=0;j<wBits;j++)dupliNum[cntPos*wBits+j]=dupliNum[p2*wBits+j];
        p2++;
      }
      if(hasLastNum&&num[cntPos]==num[cntPos-1]&&dupli[cntPos]&&dupli[cntPos-1]){
        hasLastNum=false;
        dupliPair++;
      }else hasLastNum=true;
      cntPos++;
    }
    dupliPairInLevel.set(level,dupliPair);
  }
  private void compactOneLevel(int level) { // compact half of data when numToReduce is small
    if (level == cntLevel - 1)
      calcLevelMaxSize(cntLevel + 1);
    int L1 = levelPos[level], R1 = levelPos[level + 1]; // [L,R)
//    System.out.println("T_T\t"+(R1-L1));
    if (level == 0 && !level0Sorted)
      sortLV0();
    if (R1-L1<=1) return;
    addRecord(false, num[L1], num[R1 - 1], level);
    if(((R1-L1)&1)==1){
      int reserve_P = randomForReserve.nextBoolean()?L1:(R1-1);
      long res_num = num[reserve_P];
      boolean res_dupli=dupli[reserve_P];
      long res_dupliNum = getDupliNum(level,reserve_P);
      for(int i=reserve_P;i>L1;i--){
        num[i]=num[i-1];
        dupli[i]=dupli[i-1];
        for(int j=0;j<wBits;j++)dupliNum[i*wBits+j]=dupliNum[(i-1)*wBits+j];
      }
      num[L1] = res_num;
      dupli[L1]=res_dupli;
      setDupliNum(level,L1,res_dupliNum);
      L1++;
    }
    randomlyHalveDownToLeft(level,L1, R1);
    dupliPairInLevel.set(level,0);
    int mid = (L1 + R1) >>> 1;
    mergeSortWithoutSpace(level+1, L1, mid, levelPos[level + 1], levelPos[level + 2]);
    levelPos[level + 1] = mid;
    int newP = levelPos[level + 1] - 1, oldP = L1 - 1;
    for (int i = oldP; i >= levelPos[0]; i--) {
      num[newP] = num[oldP];
      dupli[newP] = dupli[oldP];
      for(int j=0;j<wBits;j++)dupliNum[newP*wBits+j]=dupliNum[oldP*wBits+j];
      newP--;
      oldP--;
    }
    levelPos[level] = levelPos[level + 1] - (L1 - levelPos[level]);
    int numReduced = (R1 - L1) >>> 1;
    for (int i = level - 1; i >= 0; i--) levelPos[i] += numReduced;
  }

  private void quicksort(int level, int L,int R){
    if(L>=R)return;
    int mid=(L+R)>>1;
    long pNum=num[mid],tmpNum,tmpDupliNum;
    boolean pDupli=dupli[mid],tmpDupli;
    int tmpPos=L;
    for(int i=L;i<=R;i++)if(num[i]<pNum||(num[i]==pNum&&(dupli[i]||!pDupli))){
      tmpNum=num[tmpPos];
      num[tmpPos]=num[i];
      num[i]=tmpNum;
      tmpDupli=dupli[tmpPos];
      dupli[tmpPos]=dupli[i];
      dupli[i]=tmpDupli;
      tmpDupliNum=getDupliNum(level,tmpPos);
      setDupliNum(level,tmpPos,getDupliNum(level,i));
      setDupliNum(level,i,tmpDupliNum);
      tmpPos++;
    }
//    for(int i=L;i<=tmpPos-1;i++) if(num[i]>pNum)System.out.println("\t\t??!!! ERR_1 quicksort OUT OF ORDER");
//    for(int i=tmpPos;i<=R;i++) if(num[i]<pNum)System.out.println("\t\t??!!! ERR_1 quicksort OUT OF ORDER");
    quicksort(level,tmpPos,R);
//    tmpPos--;
//    while(tmpPos>L&&num[tmpPos]==num[tmpPos-1]&&dupli[tmpPos]==dupli[tmpPos-1])tmpPos--;
//    if(L<tmpPos)quicksort(L,tmpPos-1);
//    int dePos=L;
//    for(int i=L;i<tmpPos;i++)
//      if(num[i]!=pNum||dupli[i]!=pDupli){
//        tmpNum=num[dePos];
//        num[dePos]=num[i];
//        num[i]=tmpNum;
//        tmpDupli=dupli[dePos];
//        dupli[dePos]=dupli[i];
//        dupli[i]=tmpDupli;
//        dePos++;
//      }
//    quicksort(L,dePos-1);
    int samePos=tmpPos-1;
    for(int i=tmpPos-1;i>=L;i--)
      if(num[i]==pNum&&dupli[i]==pDupli){
        num[i]=num[samePos];
        dupli[i]=dupli[samePos];
        boolean[] tmpDupliNumArr = Arrays.copyOfRange(dupliNum,i*wBits,i*wBits+wBits);
        for(int j=0;j<wBits;j++)dupliNum[i*wBits+j]=dupliNum[samePos*wBits+j];
        for(int j=0;j<wBits;j++)dupliNum[samePos*wBits+j]=tmpDupliNumArr[j];
        num[samePos]=pNum;
        dupli[samePos]=pDupli;
        samePos--;
      }
    quicksort(level,L,samePos);
//    for(int i=samePos+1;i<=tmpPos-1;i++) if(num[i]!=pNum||dupli[i]!=pDupli)System.out.println("\t\t??!!! ERR_2 quicksort OUT OF ORDER");
//    for(int i=L+1;i<=R;i++) if(num[i]<num[i-1])System.out.println("\t\t??!!! quicksort OUT OF ORDER");
  }

  private void compactDupliInLevel(int level) {
    if (level == cntLevel - 1)
      calcLevelMaxSize(cntLevel + 1);
    int L1 = levelPos[level], R1 = levelPos[level + 1]; // [L,R)
    if (level == 0 && !level0Sorted)
      sortLV0();
    if (R1-L1<=1) return;
    boolean hasLastNum=false;
    long tmpNum,tmpDupliNum;
    int dupliPos=R1-1;
    for(int i=R1-1;i>=L1;i--) // put duplicate pairs at the right end
      if(hasLastNum&&num[i]==num[i+1]&&dupli[i]&&dupli[i+1]){
        hasLastNum=false;

        dupli[i+1]=dupli[dupliPos];
        dupli[dupliPos]=true;
        tmpNum=num[i+1];
        num[i+1]=num[dupliPos];
        num[dupliPos]=tmpNum;
        tmpDupliNum=getDupliNum(level,i+1);
        setDupliNum(level,i+1,getDupliNum(level,dupliPos));
        setDupliNum(level,dupliPos,tmpDupliNum);

        dupli[i]=dupli[dupliPos-1];
        dupli[dupliPos-1]=true;
        tmpNum=num[i];
        num[i]=num[dupliPos-1];
        num[dupliPos-1]=tmpNum;
        tmpDupliNum=getDupliNum(level,i);
        setDupliNum(level,i,getDupliNum(level,dupliPos-1));
        setDupliNum(level,dupliPos-1,tmpDupliNum);

        dupliPos-=2;
      }else{
        hasLastNum=true;
      }
    int disorderL=L1,disorderR=dupliPos;
    if(disorderL<disorderR){
      quicksort(level,disorderL,disorderR);
    }
    L1=dupliPos+1; // only compact duplicates (at the right end)
    randomlyHalveDownToLeft(level,L1, R1);
//    System.out.println("\t\t\t\tdupli L1,R1:\t"+L1+","+R1);
    dupliPairInLevel.set(level,0);
    int mid = (L1 + R1) >>> 1;
    mergeSortWithoutSpace(level+1, L1, mid, levelPos[level + 1], levelPos[level + 2]);
    levelPos[level + 1] = mid;
    int newP = levelPos[level + 1] - 1, oldP = L1 - 1;
    for (int i = oldP; i >= levelPos[0]; i--) {
      num[newP] = num[oldP];
      dupli[newP] = dupli[oldP];
      for(int j=0;j<wBits;j++)dupliNum[newP*wBits+j]=dupliNum[oldP*wBits+j];
      newP--;oldP--;
    }
    levelPos[level] = levelPos[level + 1] - (L1 - levelPos[level]);
    int numReduced = (R1 - L1) >>> 1;
    for (int i = level - 1; i >= 0; i--) levelPos[i] += numReduced;
  }

  boolean redoingNonFreqItems=false;

  @Override
  public void compact() {
    sortLV0();
    for (int i = 0; i < cntLevel-1; i++)
      if (getLevelSize(i) > levelMaxSize[i]) {
//        compactOneLevel(i);
        if (dupliPairInLevel.getInt(i) * 2 >= getLevelSize(i) / 2)
          compactDupliInLevel(i);
        else {
          compactOneLevel(i);
//          if(i==cntLevel-2&&!redoingNonFreqItems){
//            checkFreqItemInTopLevel();
//          }
        }
        return;
      }
//    for(int i=0;i<cntLevel;i++){
//      for(int j=levelPos[i]+1;j<levelPos[i+1];j++)
//        if(num[j]<num[j-1])System.out.println("\t\t??!!! OUT OF ORDER");
//    }
    /** overfull, need to increase height **/
//    System.out.println("\t\tcntLevel:"+cntLevel);
    if(cntLevel==1){
      calcLevelMaxSize(2);
      compactOneLevel(0);
      return;
    }
    if(lastUpdateFreqLevel!=cntLevel) {
      lastUpdateFreqLevel = cntLevel;
      checkFreqItemInTopLevel();
    }
    calcLevelMaxSize(cntLevel+1);
    compactOneLevel(cntLevel - 2);
//    this.showLevelStat();
//    this.showNum();
  }

  private void checkFreqItemInTopLevel(){
//    System.out.println("\t\tchecking checkFreqItemInTopLevel\tN:"+N);
    for (int i = 0; i < cntLevel - 1; i++)
      if (dupliPairInLevel.getInt(i) >= 1)
        compactDupliInLevel(i);
//      int freqItemCount = 0, dupliCount = 0;
//      long lastFreqItem = 0;
//      for (int i = levelPos[cntLevel - 1]; i < levelPos[cntLevel]; i++)
//        if (dupli[i]) {
//          dupliCount++;
//          if (freqItemCount == 0 || num[i] != lastFreqItem) {
//            freqItemCount++;
//            lastFreqItem = num[i];
//          }
//        }
//      System.out.println("\t\t\tfreqItemCount:" + freqItemCount + "\t" + "dupliRatioInTopLevel:" + 1.0 * dupliCount / getLevelSize(cntLevel - 1) + "\t\t");
//    for(int i=levelPos[cntLevel-2];i<levelPos[cntLevel-1];i++)
//      if(dupli[i])System.out.print("\t"+num[i]);
    for (int threshold = 1; threshold <= 4; threshold++) {
      LongArrayList freqValue = new LongArrayList(), freqCount = new LongArrayList();
      for (int i = levelPos[cntLevel - 1]; i < levelPos[cntLevel]; i++){
        int j=i;
        long dupliSum=0;
        while (j < levelPos[cntLevel] && num[j] == num[i]) {
          if(dupli[j])dupliSum+=(1L<<(cntLevel-1));
          else if(!onlyDealFullyDupli)
            dupliSum+=getDupliNum(cntLevel-1,j);
          j++;
        }
        if(dupliSum>=((long) threshold <<(cntLevel-1))){
          freqValue.add(num[i]);
          freqCount.add(0);
        }
        i=j-1;
      }
//        for (int i = levelPos[cntLevel - 1]; i < levelPos[cntLevel]; i++)
//          if (dupli[i]) {
//            int j = i + 1;
//            while (j < levelPos[cntLevel] && num[j] == num[i] && dupli[j]) j++;
//            if (j - i >= threshold) {
//              freqValue.add(num[i]);
//              freqCount.add(0);
//            }
//            i = j - 1;
//          }
//      for(int i=1;i<freqValue.size();i++)if(freqValue.getLong(i)<=freqValue.getLong(i-1))System.out.println("\t\t????freqValue??\t"+freqValue.getLong(i-1)+" "+freqValue.getLong(i));
//      System.out.println("\t\tfreqValuesInTopLevel:\t"+freqValue);
      for (int lv = 0; lv <= cntLevel - 1; lv++) {
        int tmpPos = 0;
        for (int i = levelPos[lv]; i < levelPos[lv + 1]; i++) {
          if (onlyDealFullyDupli&&!dupli[i]) continue;
          while (tmpPos < freqValue.size() && freqValue.getLong(tmpPos) < num[i]) tmpPos++;
          if (tmpPos >= freqValue.size()) break;
          if (freqValue.getLong(tmpPos) != num[i]) continue;
          freqCount.set(tmpPos, freqCount.getLong(tmpPos) + (1L << lv));
        }
      }
      long tmpFreqN = 0;
      for (long count : freqCount) tmpFreqN += count;
//      System.out.println("\t\t\t\t\ttmpFreqN:"+tmpFreqN);

      Long2LongOpenHashMap newFreqValueCount = new Long2LongOpenHashMap();
      for (int i = 0; i < freqValue.size(); i++)
        newFreqValueCount.put(freqValue.getLong(i), freqCount.getLong(i));
//      System.out.println("\t\tcheck freq value distinct:\t\t"+freqValue.size()+"\t"+newFreqValueCount.size());
      if (freqValueCount != null) {
        for (Long2LongMap.Entry entry : freqValueCount.long2LongEntrySet())
          if (!newFreqValueCount.containsKey(entry.getLongKey())) {
            if (entry.getLongValue() >= ((long) threshold << (cntLevel - 1)))
              newFreqValueCount.put(entry.getLongKey(), entry.getLongValue());
          } else
            newFreqValueCount.addTo(entry.getLongKey(), entry.getLongValue());
      }
//      System.out.println("\t\t\thashMap:" + 2 * 8 * HashCommon.arraySize(newFreqValueCount.size(), Long2LongOpenHashMap.DEFAULT_LOAD_FACTOR) + "Byte");

      if (newFreqValueCount.isEmpty()) {
        compact();
        return;
        // TODO  check non-freq items in hashmap
      }

      int tmpRestNum = (maxMemoryByte - 2 * 8 * HashCommon.arraySize(newFreqValueCount.size(), Long2LongOpenHashMap.DEFAULT_LOAD_FACTOR)) * 8 / BitsPerItem;


      int cntPos = tmpRestNum - 1;
      for (int lv = cntLevel - 1; lv >= 0; lv--)
        for (int i = levelPos[lv + 1] - 1; i >= levelPos[lv]; i--)
          if ((onlyDealFullyDupli&&!dupli[i]) || !newFreqValueCount.containsKey(num[i])) {
            cntPos--;
          }
      if(cntPos<0){
        continue;
      }
//        System.out.println("\t\t\tthreshold:\t"+threshold);

      maxMemoryNumForSketch=tmpRestNum;
      levelMaxSize = calcLevelMaxSizeByLevel(maxMemoryNumForSketch, cntLevel);
      long[] newNum = new long[maxMemoryNumForSketch];
      boolean[] newDupli = new boolean[maxMemoryNumForSketch];
      boolean[] newDupliNum = new boolean[maxMemoryNumForSketch*wBits];
      int[] newLevelPos = new int[cntLevel + 1];
//        System.out.println("\t\t???\t\t"+newDupliNum.length+"  oldDNLen:"+dupliNum.length);

      cntPos = maxMemoryNumForSketch - 1;
      for (int lv = cntLevel - 1; lv >= 0; lv--) {
        newLevelPos[lv + 1] = cntPos + 1;
        for (int i = levelPos[lv + 1] - 1; i >= levelPos[lv]; i--)
          if ((onlyDealFullyDupli&&!dupli[i]) || !newFreqValueCount.containsKey(num[i])) {
            newNum[cntPos] = num[i];
            newDupli[cntPos] = dupli[i];
            for(int j=0;j<wBits;j++)newDupliNum[cntPos*wBits+j]=dupliNum[i*wBits+j];
            cntPos--;
          }
      }
      newLevelPos[0] = cntPos + 1;
      levelPos = newLevelPos;
      num = newNum;
      dupli = newDupli;
      dupliNum=newDupliNum;
      Long2LongOpenHashMap oldMap = freqValueCount;
      freqValueCount = newFreqValueCount;
      if (oldMap != null) {
        redoingNonFreqItems=true;
        for (Long2LongMap.Entry entry : oldMap.long2LongEntrySet())
          if (!freqValueCount.containsKey(entry.getLongKey())) {
            N -= entry.getLongValue();
//          System.out.println("\t\t!!!WARN\t delete non-freq item:\t");
            for (long j = entry.getLongValue(); j > 0; j--)
              update(entry.getLongKey());
          }
        redoingNonFreqItems=false;
      }
      long checkN = 0;
      for (int i = 0; i < cntLevel; i++) checkN += (1L << i) * getLevelSize(i);
      for (Long2LongMap.Entry entry : freqValueCount.long2LongEntrySet())
        checkN += entry.getLongValue();
      if (checkN != N)
        System.out.println("\t\t\t\t[ERROR!!]\tchecking N:\t N:" + N + "\t\tcheckN:" + checkN);
      if (levelPos[0] > 0) return;
      else break;
    }
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

  private void addRecordInLevel(long minV, long maxV, int level) {
    int hisNum = compactionNumInLevel.getInt(level);
      hisNum++;
    compactionNumInLevel.set(level, hisNum);
  }

  public void addRecord(boolean isSum, long minV, long maxV, int level) {
    MIN_V=Math.min(MIN_V,minV);
    MAX_V=Math.max(MAX_V,maxV);
    if (!isSum) {
      addRecordInLevel(minV, maxV, level);
    }else{
      for(int i=0;i<=level;i++)
        addRecordInLevel(minV, maxV, i);
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

  public void sortLV0(){// dupli=true, dupliNum=0 in lv0, just sort num[]
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
    dupliPairInLevel.set(0,dupliPair);
  }

  public long getMin() {
    long mn = MIN_V;
    if(levelPos[0]<levelPos[1]) {
      if (!level0Sorted)
        sortLV0();
      mn = Math.min(mn, num[levelPos[0]]);
    }
    return mn;
  }

  public long getMax() {
    long mx = MAX_V;
    if(levelPos[0]<levelPos[1]) {
      if (!level0Sorted)
        sortLV0();
      mx = Math.max(mx, num[levelPos[1]-1]);
    }
    return mx;
  }


  public static long queryRankErrBoundGivenParameter(double sig2, long maxERR, double Pr) {
    if(maxERR==0)return 0;
    NormalDistribution dis = sig2==lastSig2&&lastNormalDis!=null?lastNormalDis:new NormalDistribution(0, Math.sqrt(sig2));
    if(sig2==lastSig2&&Pr==lastPr)return lastBound;
    double tmpBound = dis.inverseCumulativeProbability(0.5*(1-Pr));
    tmpBound=Math.min(-Math.floor(tmpBound),maxERR);
    lastBound=(long)tmpBound;
    lastPr=Pr;
    lastSig2=sig2;
    lastNormalDis=dis;
    return (long)tmpBound;
  }

  public static long queryRankErrBound(int[] relatedCompactNum, double Pr) {
    double sig2=0;
    long maxERR=0;
    for(int i=0;i<relatedCompactNum.length;i++)sig2+=0.5*relatedCompactNum[i]*Math.pow(2,i*2);
    for(int i=0;i<relatedCompactNum.length;i++)maxERR+=(long)relatedCompactNum[i]<<i;
    return queryRankErrBoundGivenParameter(sig2,maxERR,Pr);
  }

  public long queryRankErrBound(double Pr){
    return queryRankErrBound(compactionNumInLevel.toIntArray(),Pr);
  }


  boolean freqValueCountFrozen=false;
  long[] freqValue,freqCountSum;
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
        ObjectArrayList<Long2LongMap.Entry>entries = new ObjectArrayList<>(freqValueCount.size());
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
}