import it.unimi.dsi.fastutil.doubles.DoubleArrayList;
import it.unimi.dsi.fastutil.ints.IntArrayList;
import it.unimi.dsi.fastutil.longs.LongArrayList;
import it.unimi.dsi.util.XoRoShiRo128PlusRandom;
import org.apache.commons.math3.distribution.NormalDistribution;

import java.util.Arrays;
import java.util.Collections;
import java.util.List;

public class KLLDupli01 extends KLLSketchForQuantile {
  IntArrayList compactionNumInLevel,dupliPairInLevel;
  public XoRoShiRo128PlusRandom randomForReserve = new XoRoShiRo128PlusRandom();
  long MIN_V=Long.MAX_VALUE,MAX_V=Long.MIN_VALUE;
  boolean[] dupli;

  public KLLDupli01(int maxMemoryByte) {
    N = 0;
    calcParameters(maxMemoryByte);
    calcLevelMaxSize(1);
  }

  private void calcParameters(int maxMemoryByte) {
    maxMemoryNum = calcMaxMemoryNum(maxMemoryByte);
    num = new long[maxMemoryNum];
    dupli = new boolean[maxMemoryNum];
    level0Sorted = false;
    cntLevel = 0;
    compactionNumInLevel = new IntArrayList();
    dupliPairInLevel = new IntArrayList();
  }

  @Override
  protected int calcMaxMemoryNum(int maxMemoryByte) {
    return Math.min(1 << 20, maxMemoryByte*8 / (64+1));
  }

  @Override
  protected void calcLevelMaxSize(int setLevel) { // set cntLevel.  make sure cntLevel won't decrease
    int[] tmpArr = new int[setLevel + 1];
    int maxPos = cntLevel > 0 ? Math.max(maxMemoryNum, levelPos[cntLevel]) : maxMemoryNum;
    for (int i = 0; i < setLevel + 1; i++) tmpArr[i] = i < cntLevel ? levelPos[i] : maxPos;
    levelPos = tmpArr;
    for(int i=cntLevel;i<setLevel;i++) {
      compactionNumInLevel.add(0);
      dupliPairInLevel.add(0);
    }
    cntLevel = setLevel;
    levelMaxSize=calcLevelMaxSizeByLevel(maxMemoryNum,cntLevel);
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
        System.out.print(longToResult(num[j])+"("+dupli[j]+")"+", ");
      System.out.print("|\t");
    }System.out.println();
  }

  @Override
  public void update(long x) { // signed long
    if (levelPos[0] == 0) compact();
    num[--levelPos[0]] = x;
    dupli[levelPos[0]]=true;
    N++;
    level0Sorted = false;
  }


  @Override
  protected void randomlyHalveDownToLeft(int L,int R){
    int delta = getNextRand01();
    int mid = (L+R)>>>1;
    for(int i=L,j=L;i<mid;i++,j+=2) {
      num[i] = num[j + delta];
      dupli[i] = dupli[j + delta]&&dupli[j+(1-delta)]&&(num[j+delta]==num[j+(1-delta)]);
    }
    for(int i=L;i<mid;){
      int j=i,count1=0;
      while(j<mid&&num[j]==num[i]){
        count1+=dupli[j]?1:0;
        j++;
      }
      for(int k=i;k<i+count1;k++)dupli[k]=true;
      for(int k=i+count1;k<j;k++)dupli[k]=false;
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
        dupli[cntPos]=dupli[p1++];
      }
      else {
        num[cntPos]=num[p2];
        dupli[cntPos]=dupli[p2++];
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
      for(int i=reserve_P;i>L1;i--){
        num[i]=num[i-1];
        dupli[i]=dupli[i-1];
      }
      num[L1] = res_num;
      dupli[L1]=res_dupli;
      L1++;
    }
    randomlyHalveDownToLeft(L1, R1);
    dupliPairInLevel.set(level,0);
    int mid = (L1 + R1) >>> 1;
    mergeSortWithoutSpace(level+1, L1, mid, levelPos[level + 1], levelPos[level + 2]);
    levelPos[level + 1] = mid;
    int newP = levelPos[level + 1] - 1, oldP = L1 - 1;
    for (int i = oldP; i >= levelPos[0]; i--) {
      num[newP] = num[oldP];
      dupli[newP--] = dupli[oldP--];
    }
    levelPos[level] = levelPos[level + 1] - (L1 - levelPos[level]);
    int numReduced = (R1 - L1) >>> 1;
    for (int i = level - 1; i >= 0; i--) levelPos[i] += numReduced;
  }

  private void quicksort(int L,int R){
    if(L>=R)return;
    int mid=(L+R)>>1;
    long pNum=num[mid],tmpNum;
    boolean pDupli=dupli[mid],tmpDupli;
    int tmpPos=L;
    for(int i=L;i<=R;i++)if(num[i]<pNum||(num[i]==pNum&&(dupli[i]||!pDupli))){
      tmpNum=num[tmpPos];
      num[tmpPos]=num[i];
      num[i]=tmpNum;
      tmpDupli=dupli[tmpPos];
      dupli[tmpPos]=dupli[i];
      dupli[i]=tmpDupli;
      tmpPos++;
    }
    quicksort(tmpPos,R);
//    tmpPos--;
//    while(tmpPos>L&&num[tmpPos]==num[tmpPos-1]&&dupli[tmpPos]==dupli[tmpPos-1])tmpPos--;
//    if(L<tmpPos)quicksort(L,tmpPos-1);
    int dePos=L;
    for(int i=L;i<tmpPos;i++)
      if(num[i]!=pNum||dupli[i]!=pDupli){
        tmpNum=num[dePos];
        num[dePos]=num[i];
        num[i]=tmpNum;
        tmpDupli=dupli[dePos];
        dupli[dePos]=dupli[i];
        dupli[i]=tmpDupli;
        dePos++;
      }
    quicksort(L,dePos-1);
  }

  private void compactDupliInLevel(int level) {
    if (level == cntLevel - 1)
      calcLevelMaxSize(cntLevel + 1);
    int L1 = levelPos[level], R1 = levelPos[level + 1]; // [L,R)
    if (level == 0 && !level0Sorted)
      sortLV0();
    if (R1-L1<=1) return;
    boolean hasLastNum=false;
    long tmpNum;
    int dupliPos=R1-1;
    for(int i=R1-1;i>=L1;i--) // put duplicate pairs at the right end
      if(hasLastNum&&num[i]==num[i+1]&&dupli[i]&&dupli[i+1]){
        hasLastNum=false;

        dupli[i+1]=dupli[dupliPos];
        dupli[dupliPos]=true;
        tmpNum=num[i+1];
        num[i+1]=num[dupliPos];
        num[dupliPos]=tmpNum;

        dupli[i]=dupli[dupliPos-1];
        dupli[dupliPos-1]=true;
        tmpNum=num[i];
        num[i]=num[dupliPos-1];
        num[dupliPos-1]=tmpNum;

        dupliPos-=2;
      }else{
        hasLastNum=true;
      }
    int disorderL=L1,disorderR=dupliPos;
    if(disorderL<disorderR){
      //TODO quicksort.
      quicksort(disorderL,disorderR);
    }
    L1=dupliPos+1; // only compact duplicates (at the right end)
    randomlyHalveDownToLeft(L1, R1);
//    System.out.println("\t\t\t\tdupli L1,R1:\t"+L1+","+R1);
    dupliPairInLevel.set(level,0);
    int mid = (L1 + R1) >>> 1;
    mergeSortWithoutSpace(level+1, L1, mid, levelPos[level + 1], levelPos[level + 2]);
    levelPos[level + 1] = mid;
    int newP = levelPos[level + 1] - 1, oldP = L1 - 1;
    for (int i = oldP; i >= levelPos[0]; i--) {
      num[newP] = num[oldP];
      dupli[newP--] = dupli[oldP--];
    }
    levelPos[level] = levelPos[level + 1] - (L1 - levelPos[level]);
    int numReduced = (R1 - L1) >>> 1;
    for (int i = level - 1; i >= 0; i--) levelPos[i] += numReduced;
  }


  private void compactNonDupliInLevel(int level) {
    if (level == cntLevel - 1)
      calcLevelMaxSize(cntLevel + 1);
    int L1 = levelPos[level], R1 = levelPos[level + 1]; // [L,R)
    if (level == 0 && !level0Sorted)
      sortLV0();
    if (R1-L1<=1) return;
    boolean hasLastNum=false,tmpDupli;
    int lastP=-1;
    long tmpNum;
    int nonDupliPos=R1-1;
    for(int i=R1-1;i>=L1;i--) // put non-duplicate items at the right end
      if(hasLastNum&&num[i]==num[i+1]&&dupli[i]&&dupli[i+1]){
        hasLastNum=false;
      }else{
        if(hasLastNum){
          tmpDupli=dupli[lastP];
          dupli[lastP]=dupli[nonDupliPos];
          dupli[nonDupliPos]=tmpDupli;
          tmpNum=num[lastP];
          num[lastP]=num[nonDupliPos];
          num[nonDupliPos]=tmpNum;
          nonDupliPos--;
        }
        hasLastNum=true;
        lastP=i;
      }
    if(hasLastNum){
      tmpDupli=dupli[lastP];
      dupli[lastP]=dupli[nonDupliPos];
      dupli[nonDupliPos]=tmpDupli;
      tmpNum=num[lastP];
      num[lastP]=num[nonDupliPos];
      num[nonDupliPos]=tmpNum;
      nonDupliPos--;
    }
    int disorderL=L1,disorderR;
    L1=nonDupliPos+1;// not compact duplicates (at the left end)
    addRecord(false, num[L1], num[R1 - 1], level);
    if(((R1-L1)&1)==1)
      L1++;
    disorderR=L1-1;
    if(disorderL<disorderR)
      quicksort(disorderL,disorderR);
    randomlyHalveDownToLeft(L1, R1);
    int mid = (L1 + R1) >>> 1;
    mergeSortWithoutSpace(level+1, L1, mid, levelPos[level + 1], levelPos[level + 2]);
    levelPos[level + 1] = mid;
    int newP = levelPos[level + 1] - 1, oldP = L1 - 1;
    for (int i = oldP; i >= levelPos[0]; i--) {
      num[newP] = num[oldP];
      dupli[newP--] = dupli[oldP--];
    }
    levelPos[level] = levelPos[level + 1] - (L1 - levelPos[level]);
    int numReduced = (R1 - L1) >>> 1;
    for (int i = level - 1; i >= 0; i--) levelPos[i] += numReduced;
  }

  @Override
  public void compact() {
    sortLV0();
    for (int i = 0; i < cntLevel; i++)
      if (getLevelSize(i) > levelMaxSize[i]) {
//        showNum();
//        System.out.println("\t\tcompact:\tlv:"+i+"\tdupli:"+dupliPairInLevel.getInt(i) * 2+"/"+getLevelSize(i)+"="+(1.0*dupliPairInLevel.getInt(i) * 2 / getLevelSize(i)));
        if (dupliPairInLevel.getInt(i) * 2 >= getLevelSize(i) / 2)
//          compactOneLevel(i);
          compactDupliInLevel(i);
        else
          compactOneLevel(i);
//          compactNonDupliInLevel(i);
//        System.out.println("\t\tafterLVSize:"+getLevelSize(i));
//        showNum();
//        System.out.println("\n");
        return;
      }
    for (int i = 0; i < cntLevel; i++)
      if (getLevelSize(i) > levelMaxSize[i]) {
        compactOneLevel(i);
        return;
      }
    calcLevelMaxSize(cntLevel + 1);
    compactOneLevel(cntLevel - 2);
//    this.showLevelStat();
//    this.showNum();
  }

  public void merge(KLLSketchForQuantile another) {
    // todo: dupli
    if (another.cntLevel > cntLevel)
      calcLevelMaxSize(another.cntLevel);
    for (int i = 0; i < another.cntLevel; i++) {
      int numToMerge = another.levelPos[i + 1] - another.levelPos[i];
      if (numToMerge == 0) continue;
      int mergingL = another.levelPos[i];
      while (numToMerge > 0) {
        if (levelPos[0] == 0) compact();
        int delta = Math.min(numToMerge, levelPos[0]);
        if (i > 0) { // move to give space for level i
          for (int j = 0; j < i; j++) levelPos[j] -= delta;
          System.arraycopy(num, delta, num, 0, levelPos[i] - delta);
        }
        System.arraycopy(another.num, mergingL, num, levelPos[i] - delta, delta);
        levelPos[i] -= delta;
        numToMerge -= delta;
        mergingL += delta;
      }
    }
    this.N += another.N;
  }

  public void mergeWithTempSpace(KLLSketchForQuantile another,long minV,long maxV) {
    LongArrayList mn = new LongArrayList(1);mn.add(minV);
    LongArrayList mx = new LongArrayList(1);mx.add(maxV);
    mergeWithTempSpace(Collections.singletonList(another),mn,mx);
  }

  public void mergeWithTempSpace(List<KLLSketchForQuantile> otherList,LongArrayList minV,LongArrayList maxV) {
    // todo: dupli
    int[] oldLevelPos = Arrays.copyOf(levelPos, cntLevel + 1);
    int oldCntLevel = cntLevel;
    int otherNumLen = 0;
    long otherN = 0;
//    System.out.print("\t\t\t\t[mergeWithTempSpace] others:");
    for (int i=0;i<otherList.size();i++) {
      KLLSketchForQuantile another = otherList.get(i);
      if (another != null) {
//      System.out.print("\t"+another.getN());
        if (another.cntLevel > cntLevel)
          calcLevelMaxSize(another.cntLevel);
        if (another.cntLevel >= 2)
          addRecord(true, minV.getLong(i), maxV.getLong(i), another.cntLevel - 2);
        otherNumLen += another.getNumLen();
        otherN += another.getN();
      }
    }
//    System.out.println();
//    System.out.println("[mergeWithTempSpace]\totherNumLen:"+otherNumLen);
    if (getNumLen() + otherNumLen <= maxMemoryNum) {
      int cntPos = oldLevelPos[0] - otherNumLen;
      for (int i = 0; i < cntLevel; i++) {
        levelPos[i] = cntPos;
        if (i < oldCntLevel) {
          System.arraycopy(num, oldLevelPos[i], num, cntPos,
              oldLevelPos[i + 1] - oldLevelPos[i]);
          cntPos += oldLevelPos[i + 1] - oldLevelPos[i];
        }
        for (KLLSketchForQuantile another : otherList)
          if (another != null && i < another.cntLevel) {
            System.arraycopy(another.num, another.levelPos[i], num, cntPos, another.getLevelSize(i));
            cntPos += another.getLevelSize(i);
          }
        Arrays.sort(num, levelPos[i], cntPos);
//        System.out.println("\t\t!!\t"+cntPos);
      }
      levelPos[cntLevel] = cntPos;
      this.N += otherN;
    } else {
      long[] oldNum = num;
      num = new long[getNumLen() + otherNumLen];
//      System.out.println("\t\t\t\ttmp_num:"+num.length+"  old_num:"+levelPos[0]+"..."+levelPos[oldCntLevel]);
      int numLen = 0;
      for (int i = 0; i < cntLevel; i++) {
        levelPos[i] = numLen;
        if (i < oldCntLevel) {
//          System.out.println("\t\t\tlv"+i+"\toldPos:"+oldLevelPos[i]+"\t"+numLen+" this_level_old_len:"+(oldLevelPos[i + 1] - oldLevelPos[i]));
//          System.out.println("\t\t\t"+oldNum[oldLevelPos[i + 1]-1]);
          System.arraycopy(oldNum, oldLevelPos[i], num, numLen,
              oldLevelPos[i + 1] - oldLevelPos[i]);
          numLen += oldLevelPos[i + 1] - oldLevelPos[i];
        }
        for (KLLSketchForQuantile another : otherList)
          if (another != null && i < another.cntLevel) {
            System.arraycopy(another.num, another.levelPos[i], num, numLen, another.getLevelSize(i));
            numLen += another.getLevelSize(i);
          }
        Arrays.sort(num, levelPos[i], numLen);
      }
      levelPos[cntLevel] = numLen;
      this.N += otherN;
//    System.out.println("-------------------------------.............---------");
//      show();System.out.println("\t?\t"+levelPos[0]);
      while (getNumLen() > maxMemoryNum) compact();
//      show();System.out.println("\t?\t"+levelPos[0]);
//    System.out.println("\t\t??\t\t"+Arrays.toString(num));
      int newPos0 = maxMemoryNum - getNumLen();
      System.arraycopy(num, levelPos[0], oldNum, newPos0, getNumLen());
      for (int i = cntLevel; i >= 0; i--) levelPos[i] += newPos0 - levelPos[0];
      num = oldNum;
    }
//    System.out.println("\t\t??\t\t"+Arrays.toString(num));
//    System.out.println("\t\t??\t\t"+Arrays.toString(levelPos));
//    System.out.println("-------------------------------.............---------");
//    System.out.println("[MERGE result]");
//    show();
//    System.out.println();
  }

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
//    NormalDistribution dis = sig2==lastSig2&&lastNormalDis!=null?lastNormalDis:new NormalDistribution(0, Math.sqrt(sig2));
//    if(sig2==lastSig2&&Pr==lastPr)return lastBound;
//    double tmpBound = dis.inverseCumulativeProbability(0.5*(1-Pr));
//    tmpBound=Math.min(-Math.floor(tmpBound),maxERR);
//    lastBound=(int)tmpBound;
//    lastPr=Pr;
//    lastSig2=sig2;
//    lastNormalDis=dis;
////    System.out.println("\tPr:"+Pr+"\t"+"sig2:\t"+sig2+"\t\tbyInvCum:\t"+(long)tmpBound);
//    return (int)tmpBound;
  }

  public long queryRankErrBound(double Pr){
    return queryRankErrBound(compactionNumInLevel.toIntArray(),Pr);
  }

  // 返回误差err，该数值在sketch里估计排名的偏差绝对值有Pr的概率<err
  public long queryRankErrBound(long result, double Pr) {
    int[] relatedCompactNum = getRelatedCompactNum();
    return queryRankErrBound(relatedCompactNum,Pr);
  }


  private long getLowerBound(long queryRank, double Pr) {
    long L = getMin()-1, R = getMax()-1, mid;
    int approxRankL = getApproxRank(L);
    while (L < R) {
      mid = L + ((R - L + 1) >>> 1);
      assert L <= mid && mid <= R;
      int approxRankMid = getApproxRank(mid);
      if(approxRankMid==approxRankL){L=mid;continue;}

      long rankErrBound = queryRankErrBound(mid, Pr);
//      System.out.println("\t\t\t\t\t"+longToResult(mid)+"\t\trank:"+approxRankMid+"\t\terr:"+rankErrBound+"\t\t\tL,R,mid:"+L+" "+R+" "+mid);
      if (approxRankMid + rankErrBound < /*<*/ queryRank) {
        L = mid;
        approxRankL = approxRankMid;
      }

      else R = mid - 1;
    }
    L++;
//    if(getApproxRank(L)+queryBound(L,Pr)<queryRank)L++;
//    System.out.println("\t\t[]exactKLL lowerBound.\t\tPr:"+Pr+"\trank:"+queryRank+"\t\t\t\tlowerBoundV:"+L+"(longToResult:"+longToResult(L)+")"+ "\t\tL_rank:"+getApproxRank(L)+"\t\tPrErr:"+queryRankErrBound(L,Pr)+"\t1.0Err:"+queryRankErrBound(L,1.0));
    return L;
  }

  private long getUpperBound(long queryRank, double Pr) {
    long L = getMin(), R = getMax(), mid;
    int approxRankR = getApproxRank(R);
    while (L < R) {
      mid = L + ((R - L) >>> 1);
      int approxRankMid = getApproxRank(mid);
      if(approxRankMid==approxRankR){R=mid;continue;}
      long rankErrBound = queryRankErrBound(mid, Pr);
      if (approxRankMid - rankErrBound >= queryRank){
        R = mid;
        approxRankR = approxRankMid;
      }
      else L = mid + 1;
    }
//    L--;
//    System.out.println("\t\t[]exactKLL upperBound.\t\tPr:"+Pr+"\trank:"+queryRank+"\t\t\t\tupperBoundV:"+L+"(longToResult:"+longToResult(L)+")"+"\t\tR_rank:"+getApproxRank(L)+"\t\tePrErr:"+queryRankErrBound(L,Pr)+"\t1.0Err:"+queryRankErrBound(L,1.0));
//    System.out.println("\t\t\t\trank()");
    return L;
  }

  public double[] findResultRange(long K1, long K2, double Pr) {
    DoubleArrayList result = new DoubleArrayList();
    long valL, valR;
    if (exactResult()) {
      valL = getExactResult((int) K1 - 1);
      valR = getExactResult((int) K2 - 1);
//      System.out.println("\t\tEXACT RESULT!\tn:"+getN()+"\tK1,2:"+K1+","+K2+"\t\tvalL,R:"+valL+","+valR);
      result.add(longToResult(valL));
      result.add(longToResult(valR));
      result.add(-233);
    } else
    if(getMin()==getMax()){
      valL = valR = getMin();
//      System.out.println("\t\tEXACT RESULT!\tn:"+getN()+"\tK1,2:"+K1+","+K2+"\t\tvalL,R:"+valL+","+valR);
      result.add(longToResult(valL));
      result.add(longToResult(valR));
      result.add(-233);
    }else{
//      queriedBound.clear();
      valL = getLowerBound(K1, Pr);
      valR = getUpperBound(K2, Pr);
      result.add(longToResult(valL));
      result.add(longToResult(valR));
      result.add(getApproxRank(valL));
      result.add(getApproxRank(valR));
      result.add(queryRankErrBound(valL,Pr));
      result.add(queryRankErrBound(valR,Pr));
    }
    return result.toDoubleArray();
  }

  public double[] getFilterL(long CountOfValL,long CountOfValR,double valL,double valR,long K,double Pr){
    if(K<=CountOfValL)return new double[]{valL,-233.0};
    if(K>CountOfValL+getN())return new double[]{valR,-233.0};
    K-=CountOfValL;
    if(exactResult()){
      return new double[]{longToResult(getExactResult((int) K - 1)),-233.0};
    }
    if(getMin()==getMax())return new double[]{longToResult(getMin()),-233.0};
    return new double[]{longToResult(getLowerBound(K, Pr))};
  }

  public double[] getFilterR(long CountOfValL,long CountOfValR,double valL,double valR,long K,double Pr){
    if(K<=CountOfValL)return new double[]{valL,-233.0};
    if(K>CountOfValL+getN())return new double[]{valR,-233.0};
    K-=CountOfValL;
    if(exactResult()){
      return new double[]{longToResult(getExactResult((int) K - 1)),-233.0};
    }
    if(getMin()==getMax())return new double[]{longToResult(getMin()),-233.0};
    return new double[]{longToResult(getUpperBound(K, Pr))};
  }

  public double[] getFilter(long CountOfValL,long CountOfValR,double valL,double valR,long K1,long K2,double Pr){
    double[] filterL =getFilterL(CountOfValL,CountOfValR,valL,valR,K1,Pr);
    double[] filterR =getFilterR(CountOfValL,CountOfValR,valL,valR,K2,Pr);
    if(filterL.length+filterR.length==4)return new double[]{filterL[0],filterR[0],-233};
    else return new double[]{filterL[0],filterR[0]};
  }

}