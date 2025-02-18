import it.unimi.dsi.fastutil.doubles.DoubleArrayList;
import it.unimi.dsi.fastutil.ints.Int2DoubleOpenHashMap;
import it.unimi.dsi.fastutil.ints.IntIntPair;
import it.unimi.dsi.fastutil.longs.LongArrayList;
import it.unimi.dsi.fastutil.objects.Object2DoubleOpenHashMap;
import it.unimi.dsi.fastutil.objects.ObjectArrayList;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.util.Arrays;
import java.util.Date;
import java.util.Random;


public class MainForCompareOurAndApprox {
    int dataType = -233;
    static int startType = 2, endType = 2;
    static int[] Ns=new int[]{(int)3e7,(int)3e7,(int)3e7,(int)1.1e8};
    public static int TEST_CASE = 1,QUERY_PER_TEST=1; // CHECK IT  1e8:430/case
    static double[] a;
    static double LINEAR=1.0;
    ObjectArrayList<StringBuilder>RESULT=new ObjectArrayList<>();
    int RESULT_LINE = 0;
    static Object2DoubleOpenHashMap<IntIntPair> typeMem2Alpha;
    static Int2DoubleOpenHashMap type2MinV;


    public void prepareA(int dataType) throws IOException {
        if (a == null || a.length!=Ns[dataType]) a = new double[Ns[dataType]];
        this.dataType = dataType;
        BufferedReader reader = null;
        if (dataType == 1) reader = new BufferedReader(new FileReader(new File("DupliTorqueVoltage.txt")));
        if (dataType == 2) reader = new BufferedReader(new FileReader(new File("DupliECommercePrice.txt")));
        assert reader != null;
        reader.readLine(); // ignore first line.
        String line;
        int cntN = 0;
        while ((line = reader.readLine()) != null) {
            a[cntN++] = Double.parseDouble(line);
            if (cntN == Ns[dataType]) break;
        }
    }

    public int getValueActualRank(double[] sortedA, int queryN, double v) { // number of elements <= v
        int L = 0, R = queryN - 1;
        while (L < R) {
            int mid = (L + R + 1) >>> 1;
            if (v < sortedA[mid]) R = mid - 1;
            else L = mid;
        }
        return L;
    }

    public int getValueLessThan(double[] sortedA, int queryN, double v) { // number of elements <= v
        int L = 0, R = queryN - 1;
        while (L < R) {
            int mid = (L + R + 1) >>> 1;
            if (sortedA[mid] < v) L = mid;
            else R = mid - 1;
        }
        return sortedA[L] < v ? L : L - 1;
    }

    public int getDeltaRank(double[] sortedA, int queryN, double v, int targetRank) {
        int rank_L = getValueLessThan(sortedA, queryN, v) + 1;
        int rank_R = getValueActualRank(sortedA, queryN, v);
//        System.out.println("\t\t\t"+targetRank+"\t\tresultLR:"+rank_L+"..."+rank_R+"\t\tresV:"+v);
        if (targetRank >= rank_L && targetRank <= rank_R) return 0;
        else return targetRank < rank_L ? (targetRank - rank_L) : (targetRank - rank_R);
    }

    private long dataToLong(double data) {
        long result = Double.doubleToLongBits((double) data);
        return data >= 0d ? result : result ^ Long.MAX_VALUE;
    }

    private double longToResult(long result) {
        result = (result >>> 63) == 0 ? result : result ^ Long.MAX_VALUE;
        return Double.longBitsToDouble(result);
    }
    public static LongArrayList getRankList(double[] sortedA,DoubleArrayList values){
        LongArrayList result=new LongArrayList();
        int pos=0,lastPos=0,idealBucketSize=(sortedA.length+values.size()-1)/values.size();
        double lastValue=Double.MAX_VALUE;
        for(double value:values){
//            while(pos<sortedA.length&&sortedA[pos]<=value)pos++;
            while(pos<sortedA.length&&sortedA[pos]<value)pos++;
            while(pos<sortedA.length&&sortedA[pos]==value&&pos-lastPos<idealBucketSize)pos++;
            result.add(pos);
            lastValue=value;lastPos=pos;
        }
        result.add(sortedA.length);
//        for(int i=0;i<values.size();i++){
//            int j=i+1;
//            while(j<values.size()&&values.getDouble(i)==values.getDouble(j))j++;
////            j--;
//            if(j>=i+3){
//                long avgDelta=(result.getLong(j)-result.getLong(i))/(j-i),maxDelta=0;
//                for(int k=i+1;k<=j;k++)maxDelta=Math.max(maxDelta,result.getLong(k)-result.getLong(k-1));
//                if(avgDelta!=idealBucketSize&&maxDelta>idealBucketSize)System.out.println("\t\tavgDelta:"+avgDelta+"\t\tideal:"+idealBucketSize+"\t\tmaxDelta:"+maxDelta);
//                for(int k=i+1;k<=j;k++)result.set(k,result.getLong(k-1)+avgDelta);
//            }
//            i=j-1;
//        }
        return result;
    }

    public void testMaxPartition(int totN, int sketchSizeByte, int bucketNum){
        int[] LL = new int[TEST_CASE];
        int[] RR = new int[TEST_CASE];
        Random random = new Random();
        for (int i = 0; i < TEST_CASE; i++) {
            LL[i] = random.nextInt(Ns[dataType] - totN + 1);
            RR[i] = LL[i] + totN;
        }

        double avgExactMaxPartition=0,avgKLLMaxPartition=0,avgKLLDupliMaxPartition=0,avgDDMaxPartition=0,avgDSSMaxPartition=0,avgTDMaxPartition=0,avgGKMaxPartition=0;

        for (int T = 0; T < TEST_CASE; T++) {
            int L = LL[T], R = RR[T];

            double DataMinV=Double.MAX_VALUE,DataMaxV=-DataMinV;
            for (int i = L; i < R; i++){
                DataMinV=Math.min(DataMinV,a[i]);
                DataMaxV=Math.max(DataMaxV,a[i]);
            }
            double dataset_V = DataMaxV-DataMinV+1;
            int DDLimit = sketchSizeByte/DDSketchPositiveForDupli.bucketNumPerByteInMemoryForApprox;
            double DDSketch_GAMMA = Math.pow(10, Math.log10(dataset_V) / (DDLimit - 1));
            double DDSketch_ALPHA = 1 - 2 / (DDSketch_GAMMA + 1);
//            System.out.println("\t\t\tDD gamma:"+DDSketch_GAMMA);

            KLLSketchLazyExactPriori KLLWorker = new KLLSketchLazyExactPriori(sketchSizeByte);
            KLLDupliPair KLLDupliWorker = new KLLDupliPair(sketchSizeByte);
            DDSketchPositiveForDupli DDWorker = new DDSketchPositiveForDupli(DDSketch_ALPHA,DDLimit);
            DyadicSpaceSavingForDupli DSSWorker = new DyadicSpaceSavingForDupli(sketchSizeByte);
            TDigestForDupli TDWorker = new TDigestForDupli(sketchSizeByte, 1);
            GKBandForDupli GKWorker = new GKBandForDupli(sketchSizeByte);
            DoubleArrayList divideQuantiles = new DoubleArrayList();


            for (int i = L; i < R; i++) {
                KLLWorker.update(dataToLong(a[i]));
                KLLDupliWorker.update(dataToLong(a[i]));
                DDWorker.update(a[i]-DataMinV+1);
                DSSWorker.update(a[i]);
                TDWorker.update(a[i]);
                GKWorker.update(a[i]);
            }
            //System.out.println("--- Compacts:");KLLWorker.showCompact();KLLDupliWorker.showCompact();

            DDSketchPositiveForDupli tmpDD=DDWorker;
            while(tmpDD.sketch_size()<=DDLimit){
                DDWorker=tmpDD;
                DDSketch_ALPHA*=0.9;
                tmpDD=new DDSketchPositiveForDupli(DDSketch_ALPHA,DDLimit);
                for(int i=L;i<R;i++)tmpDD.update(a[i]-DataMinV+1);
            }

            LongArrayList divideRanks=new LongArrayList();
            for(int i=0;i<bucketNum-1;i++){
                divideQuantiles.add(1.0/bucketNum*(i+1));
                divideRanks.add(totN/bucketNum*(i+1L)+(Math.min(i, totN % bucketNum)));
            }
            DoubleArrayList KLLDivideValues = new DoubleArrayList(),KLLDupliDivideValues = new DoubleArrayList(),DDDivideValues=new DoubleArrayList(),TDDivideValues=new DoubleArrayList(),DSSDivideValues=new DoubleArrayList(),GKDivideValues=new DoubleArrayList();
            GKDivideValues=GKWorker.query(divideQuantiles);
            for(int i=0;i<bucketNum-1;i++) {
                KLLDivideValues.add(KLLWorker.longToResult(KLLWorker.findMinValueWithRank(divideRanks.getLong(i))));
                KLLDupliDivideValues.add(KLLDupliWorker.longToResult(KLLDupliWorker.findMinValueWithRank(divideRanks.getLong(i))));
                DDDivideValues.add(DDWorker.getQuantile(divideQuantiles.getDouble(i))+DataMinV-1);
                DSSDivideValues.add(DSSWorker.getQuantile(divideQuantiles.getDouble(i)));
                TDDivideValues.add(TDWorker.getQuantile(divideQuantiles.getDouble(i)));
//                GKDivideValues.add(GKWorker.query(divideQuantiles.getDouble(i)));

//                System.out.println("\t\t\tq:"+divideQuantiles.getDouble(i)+"\tKLL_Value:"+KLLDivideValues.getDouble(i)+"\tDD_Value:"+DDDivideValues.getDouble(i)+"\tTD_Value:"+TDDivideValues.getDouble(i)+"\tGK_Value:"+GKDivideValues.getDouble(i));
            }

            int ExactMaxPartition = (totN+bucketNum-1)/bucketNum;
            double[] query_a=new double[totN];
            System.arraycopy(a, L, query_a, 0, R - L);
            Arrays.sort(query_a);
//            System.out.println("\t\t___________----- KLL,KLLDupli:");
            LongArrayList KLLDivideRanks = getRankList(query_a,KLLDivideValues);
            LongArrayList KLLDupliDivideRanks = getRankList(query_a,KLLDupliDivideValues);
//            System.out.println("\t\t___________----- OVER");
            LongArrayList DDDivideRanks = getRankList(query_a,DDDivideValues);
            LongArrayList DSSDivideRanks = getRankList(query_a,DSSDivideValues);
            LongArrayList TDDivideRanks = getRankList(query_a,TDDivideValues);
            LongArrayList GKDivideRanks = getRankList(query_a,GKDivideValues);

            long KLLMaxPartition=KLLDivideRanks.getLong(0)-1;
            long KLLDupliMaxPartition=KLLDupliDivideRanks.getLong(0)-1;
            long DDMaxPartition=DDDivideRanks.getLong(0)-1;
            long DSSMaxPartition=DSSDivideRanks.getLong(0)-1;
            long TDMaxPartition=TDDivideRanks.getLong(0)-1;
            long GKMaxPartition=GKDivideRanks.getLong(0)-1;
            for(int i=0;i<bucketNum-1;i++){
                KLLMaxPartition=Math.max(KLLMaxPartition,KLLDivideRanks.getLong(i+1)-KLLDivideRanks.getLong(i));
                KLLDupliMaxPartition=Math.max(KLLDupliMaxPartition,KLLDupliDivideRanks.getLong(i+1)-KLLDupliDivideRanks.getLong(i));
                DDMaxPartition=Math.max(DDMaxPartition,DDDivideRanks.getLong(i+1)-DDDivideRanks.getLong(i));
                DSSMaxPartition=Math.max(DSSMaxPartition,DSSDivideRanks.getLong(i+1)-DSSDivideRanks.getLong(i));
                TDMaxPartition=Math.max(TDMaxPartition,TDDivideRanks.getLong(i+1)-TDDivideRanks.getLong(i));
                GKMaxPartition=Math.max(GKMaxPartition,GKDivideRanks.getLong(i+1)-GKDivideRanks.getLong(i));
//                if(GKDivideRanks.getLong(i+1)-GKDivideRanks.getLong(i)==GKMaxPartition)
//                    System.out.println("\t\tGK relaPos:"+1.0*(i+1)/(bucketNum-1));
            }
//            System.out.println("\t\tDD sketch size:"+ DDWorker.sketch_size());

            avgExactMaxPartition+=ExactMaxPartition;
            avgKLLMaxPartition+=KLLMaxPartition;
            avgKLLDupliMaxPartition+=KLLDupliMaxPartition;
            avgDDMaxPartition+=DDMaxPartition;
            avgDSSMaxPartition+=DSSMaxPartition;
            avgTDMaxPartition+=TDMaxPartition;
            avgGKMaxPartition+=GKMaxPartition;

        }

        String line = "";
        line += "\t\t"+bucketNum+"\t"+avgExactMaxPartition/TEST_CASE+"\t"+avgKLLMaxPartition/TEST_CASE+"\t"+avgKLLDupliMaxPartition/TEST_CASE+"\t"+avgDDMaxPartition/TEST_CASE+"\t"+avgDSSMaxPartition/TEST_CASE+"\t"+avgTDMaxPartition/TEST_CASE+"\t"+avgGKMaxPartition/TEST_CASE;
        RESULT.get(RESULT_LINE).append(line);
        System.out.println("\t\t"+line);
    }

    public void testSquarePartition(int totN){
//        int sketchSizeByte=totN/250*8/*totN/200*8*//*1024*256*/, bucketNum=totN/50; // totN/100*8: too large count
        int sketchSizeByte=1024*256, bucketNum=totN/10;
//        int sketchSizeByte=totN/10*8, bucketNum=totN/10; // OLD
        int[] LL = new int[TEST_CASE];
        int[] RR = new int[TEST_CASE];
        Random random = new Random();
        for (int i = 0; i < TEST_CASE; i++) {
            LL[i] = random.nextInt(Ns[dataType] - totN + 1);
            RR[i] = LL[i] + totN;
        }

        double avgExactSquareSum=0,avgKLLSquareSum=0,avgKLLDupliSquareSum=0,avgDDSquareSum=0,avgDSSSquareSum=0,avgTDSquareSum=0,avgGKSquareSum=0;
        long ST_TIME=new Date().getTime();
        for (int T = 0; T < TEST_CASE; T++) {
            int L = LL[T], R = RR[T];

            double MIN_V = type2MinV.getOrDefault(dataType, 0);
            double DataMinV=Double.MAX_VALUE,DataMaxV=-DataMinV;
            for (int i = L; i < R; i++){
                DataMinV=Math.min(DataMinV,a[i]);
                DataMaxV=Math.max(DataMaxV,a[i]);
            }
            double dataset_V = DataMaxV-DataMinV+1;
            int DDLimit = sketchSizeByte/DDSketchPositiveForDupli.bucketNumPerByteInMemoryForApprox;
            double DDSketch_GAMMA = Math.pow(10, Math.log10(dataset_V) / (DDLimit - 1));
            double DDSketch_ALPHA = 1 - 2 / (DDSketch_GAMMA + 1);
            System.out.println("\t\t\tDD gamma:"+DDSketch_GAMMA);
            double DD_ALPHA = DDSketch_ALPHA;//typeMem2Alpha.getDouble(IntIntPair.of(dataType, sketchSizeByte / 1024));

            KLLSketchLazyExactPriori KLLWorker = new KLLSketchLazyExactPriori(sketchSizeByte);
            KLLDupliPair KLLDupliWorker = new KLLDupliPair(sketchSizeByte);
            DDSketchPositiveForDupli DDWorker = new DDSketchPositiveForDupli(DD_ALPHA,DDLimit);
            DyadicSpaceSavingForDupli DSSWorker = new DyadicSpaceSavingForDupli(sketchSizeByte);
            TDigestForDupli TDWorker = new TDigestForDupli(sketchSizeByte, 1);
            GKBandForDupli GKWorker = new GKBandForDupli(sketchSizeByte);
            DoubleArrayList divideQuantiles = new DoubleArrayList();

            for (int i = L; i < R; i++) {
                KLLWorker.update(dataToLong(a[i]));
                KLLDupliWorker.update(dataToLong(a[i]));
                DDWorker.update(a[i]-MIN_V+1);
                DSSWorker.update(a[i]);
                TDWorker.update(a[i]);
                GKWorker.update(a[i]);
            }
            System.out.println("--- Compacts:");KLLWorker.showCompact();KLLDupliWorker.showCompact();

            DDSketchPositiveForDupli tmpDD=DDWorker;
            while(tmpDD.sketch_size()<=DDLimit){
                DDWorker=tmpDD;
                DD_ALPHA*=0.95;
                tmpDD=new DDSketchPositiveForDupli(DD_ALPHA,DDLimit);
                for(int i=L;i<R;i++)tmpDD.update(a[i]-MIN_V+1);
            }DD_ALPHA/=0.95;//System.out.println("\t\t\tDD ratio:"+tmpDD.sketch_size()/DDLimit*1.0+"\tpositive_collapse_bound:"+DDWorker.positive_collapse_bound);
            System.out.println("\t\t\tfinish update. time:"+(new Date().getTime()-ST_TIME)/1000.0);

            LongArrayList divideRanks=new LongArrayList();
            for(int i=0;i<bucketNum-1;i++){
                divideQuantiles.add(1.0/bucketNum*(i+1));
                divideRanks.add(totN/bucketNum*(i+1L)+(Math.min(i, totN % bucketNum)));
            }
            DoubleArrayList KLLDivideValues = new DoubleArrayList(),KLLDupliDivideValues = new DoubleArrayList(),DDDivideValues=new DoubleArrayList(),DSSDivideValues=new DoubleArrayList(),TDDivideValues=new DoubleArrayList(),GKDivideValues;

            DDDivideValues=DDWorker.getQuantiles(divideQuantiles);
            TDDivideValues=TDWorker.quantiles(divideQuantiles);
            GKDivideValues=GKWorker.query(divideQuantiles);
            LongArrayList dupliLongs = KLLDupliWorker.findMinValuesWithRanks(divideRanks);
            for(long l:dupliLongs)KLLDupliDivideValues.add(KLLDupliWorker.longToResult(l));
//            System.out.println("\t\t\tgot DD,TD divide values. time:"+(new Date().getTime()-ST_TIME)/1000.0);
            for(int i=0;i<bucketNum-1;i++) {
                KLLDivideValues.add(KLLWorker.longToResult(KLLWorker.findMinValueWithRank(divideRanks.getLong(i))));
                DDDivideValues.set(i,DDDivideValues.getDouble(i)+MIN_V-1);
                DSSDivideValues.add(DSSWorker.getQuantile(divideQuantiles.getDouble(i)));
                TDDivideValues.add(TDWorker.getQuantile(divideQuantiles.getDouble(i)));
//                System.out.println("\t\t\tq:"+divideQuantiles.getDouble(i)+"\tKLL_Value:"+KLLDivideValues.getDouble(i)+"\tDD_Value:"+DDDivideValues.getDouble(i)+"\tTD_Value:"+TDDivideValues.getDouble(i));
            }
            System.out.println("\t\t\tgot divide values. time:"+(new Date().getTime()-ST_TIME)/1000.0);
            double ExactSquareSum = Math.pow((totN+bucketNum-1)/bucketNum,2)*bucketNum;
            double[] query_a=new double[totN];
            System.arraycopy(a, L, query_a, 0, R - L);
            Arrays.sort(query_a);

            LongArrayList KLLDivideRanks = getRankList(query_a,KLLDivideValues);
            LongArrayList KLLDupliDivideRanks = getRankList(query_a,KLLDupliDivideValues);
            LongArrayList DDDivideRanks = getRankList(query_a,DDDivideValues);
            LongArrayList DSSDivideRanks = getRankList(query_a,DSSDivideValues);
            LongArrayList TDDivideRanks = getRankList(query_a,TDDivideValues);
            LongArrayList GKDivideRanks = getRankList(query_a,GKDivideValues);
            System.out.println("\t\t\tgot rankList. time:"+(new Date().getTime()-ST_TIME)/1000.0);

            double KLLSquareSum=Math.pow(KLLDivideRanks.getLong(0)-1,2);
            double KLLDupliSquareSum=Math.pow(KLLDupliDivideRanks.getLong(0)-1,2);
            double DSSSquareSum=Math.pow(DSSDivideRanks.getLong(0)-1,2);
            double DDSquareSum=Math.pow(DDDivideRanks.getLong(0)-1,2);
            double TDSquareSum=Math.pow(TDDivideRanks.getLong(0)-1,2);
            double GKSquareSum=Math.pow(GKDivideRanks.getLong(0)-1,2);
            for(int i=0;i<bucketNum-1;i++){
                KLLSquareSum+=Math.pow(KLLDivideRanks.getLong(i+1)-KLLDivideRanks.getLong(i),2);
                KLLDupliSquareSum+=Math.pow(KLLDupliDivideRanks.getLong(i+1)-KLLDupliDivideRanks.getLong(i),2);
                DSSSquareSum+=Math.pow(DSSDivideRanks.getLong(i+1)-DSSDivideRanks.getLong(i),2);
                DDSquareSum+=Math.pow(DDDivideRanks.getLong(i+1)-DDDivideRanks.getLong(i),2);
                TDSquareSum+=Math.pow(TDDivideRanks.getLong(i+1)-TDDivideRanks.getLong(i),2);
                GKSquareSum+=Math.pow(GKDivideRanks.getLong(i+1)-GKDivideRanks.getLong(i),2);
            }
//            System.out.println("\t\tDD sketch size:"+ DDWorker.sketch_size());

            avgExactSquareSum+=ExactSquareSum;
            avgKLLSquareSum+=KLLSquareSum;
            avgKLLDupliSquareSum+=KLLDupliSquareSum;
            avgDSSSquareSum+=DSSSquareSum;
            avgDDSquareSum+=DDSquareSum;
            avgTDSquareSum+=TDSquareSum;
            avgGKSquareSum+=GKSquareSum;

        }

        String line = "";
        line += "\t\tN,M(KB),buckets:\t"+totN+","+sketchSizeByte/1024.0+"KB,"+bucketNum+"\t"+avgExactSquareSum/TEST_CASE+"\t"+avgKLLSquareSum/TEST_CASE+"\t"+avgKLLDupliSquareSum/TEST_CASE+"\t"+avgDDSquareSum/TEST_CASE+"\t"+avgDSSSquareSum/TEST_CASE+"\t"+avgTDSquareSum/TEST_CASE+"\t"+avgGKSquareSum/TEST_CASE+"\t\t|";
        RESULT.get(RESULT_LINE).append(line);
        System.out.println("\t\t"+line);
    }

    public double getRangeError(double[] sortedA,int rankL,int rankR,LongArrayList rank,DoubleArrayList value,double minV,double maxV){
        long lastRank=1;
        double valL=sortedA[rankL],valR=sortedA[rankR];
        double ansRankL=0,ansRankR=0;
        for(int i=0;i<rank.size();i++){
            if(lastRank<=rankL&&rankL<=rank.getLong(i)){
                double v1=i==0?minV:value.getDouble(i-1),v2=i==rank.size()-1?maxV:value.getDouble(i);
                ansRankL=lastRank+(rank.getLong(i)-lastRank)*(valL-v1)/(v2-v1);
            }
            if(lastRank<=rankR&&rankR<=rank.getLong(i)){
                double v1=i==0?minV:value.getDouble(i-1),v2=i==rank.size()-1?maxV:value.getDouble(i);
                ansRankR=lastRank+(rank.getLong(i)-lastRank)*(valR-v1)/(v2-v1);
            }
            lastRank=rank.getLong(i);
        }
        double esti=ansRankR-ansRankL;
        return Math.abs(esti-(rankR-rankL))/(rankR-rankL);
    }

    public double getRangeError(double[] sortedA,double valL,double valR,LongArrayList rank,DoubleArrayList value,double minV,double maxV){
        double ansRankL=0,ansRankR=0;
        LongArrayList anoRank=new LongArrayList();anoRank.add(0);anoRank.addAll(rank);
        DoubleArrayList anoValue=new DoubleArrayList();anoValue.add(minV);anoValue.addAll(value);anoValue.add(maxV);
        for(int i=1;i<anoRank.size();i++){
            if(anoValue.getDouble(i-1)<=valL&&valL<=anoValue.getDouble(i)){
                double v1=anoValue.getDouble(i-1),v2=anoValue.getDouble(i);
                ansRankL=anoRank.getLong(i-1)+LINEAR*(anoRank.getLong(i)-anoRank.getLong(i-1))*(valL-v1)/(v2-v1);
            }
            if(anoValue.getDouble(i-1)<=valR&&valR<=anoValue.getDouble(i)){
                double v1=anoValue.getDouble(i-1),v2=anoValue.getDouble(i);
                ansRankR=anoRank.getLong(i-1)+LINEAR*(anoRank.getLong(i)-anoRank.getLong(i-1))*(valR-v1)/(v2-v1);
            }
        }
        double esti=ansRankR-ansRankL,ans=getValueActualRank(sortedA, sortedA.length, valR)-getValueActualRank(sortedA, sortedA.length, valL);
//        System.out.println("\t\testi:"+esti+"\t\tans:"+ans+"\t\tvalL,R:"+valL+","+valR+"\t\tdataMinMaxV:"+minV+","+maxV);
        return Math.abs(esti-ans);// /sortedA.length;
    }

    public void testCardinality(int totN,int sketchSizeByte, int bucketNum){
        int[] LL = new int[TEST_CASE];
        int[] RR = new int[TEST_CASE];
        Random random = new Random();
        for (int i = 0; i < TEST_CASE; i++) {
            LL[i] = random.nextInt(Ns[dataType] - totN + 1);
            RR[i] = LL[i] + totN;
        }

        double avgExact=0,avgKLL=0,avgDD=0,avgTD=0;
        DoubleArrayList exacts=new DoubleArrayList(),KLLs=new DoubleArrayList(),DDs=new DoubleArrayList(),TDs=new DoubleArrayList();

        for (int T = 0; T < TEST_CASE; T++) {
            int L = LL[T], R = RR[T];

            double DataMinV=Double.MAX_VALUE,DataMaxV=-DataMinV;
            for (int i = L; i < R; i++){
                DataMinV=Math.min(DataMinV,a[i]);
                DataMaxV=Math.max(DataMaxV,a[i]);
            }
            double dataset_V = DataMaxV-DataMinV+1;
            int DDLimit = sketchSizeByte/DDSketchPositiveForDupli.bucketNumPerByteInMemoryForExact;
            double DDSketch_GAMMA = Math.pow(10, Math.log10(dataset_V) / (DDLimit - 1));
            double DDSketch_ALPHA = 1 - 2 / (DDSketch_GAMMA + 1);
            KLLSketchLazyExactPriori KLLWorker = new KLLSketchLazyExactPriori(sketchSizeByte);


            DDSketchPositiveForDupli DDWorker = new DDSketchPositiveForDupli(DDSketch_ALPHA,DDLimit);
            TDigestForDupli TDWorker = new TDigestForDupli(sketchSizeByte, 1);
            DoubleArrayList divideQuantiles = new DoubleArrayList();

            for (int i = L; i < R; i++) {
                KLLWorker.update(dataToLong(a[i]));
                DDWorker.update(a[i]-DataMinV+1);
                TDWorker.update(a[i]);
            }

            LongArrayList divideRanks=new LongArrayList();
            for(int i=0;i<bucketNum-1;i++){
                divideQuantiles.add(1.0/bucketNum*(i+1));
                divideRanks.add(totN/bucketNum*(i+1L)+(Math.min(i, totN % bucketNum)));
            }
            DoubleArrayList ExactDividedValues=new DoubleArrayList(),KLLDivideValues = new DoubleArrayList(),DDDivideValues=new DoubleArrayList(),TDDivideValues=new DoubleArrayList();

            double[] query_a=new double[totN];
            System.arraycopy(a, L, query_a, 0, R - L);
            Arrays.sort(query_a);
            for(int i=0;i<bucketNum-1;i++) {
                ExactDividedValues.add(query_a[(int)divideRanks.getLong(i)]);
                KLLDivideValues.add(KLLWorker.longToResult(KLLWorker.findMinValueWithRank(divideRanks.getLong(i))));
                DDDivideValues.add(DDWorker.getQuantile(divideQuantiles.getDouble(i))+DataMinV-1);
                TDDivideValues.add(TDWorker.getQuantile(divideQuantiles.getDouble(i)));

//                System.out.println("\t\t\tq:"+divideQuantiles.getDouble(i)+"\tKLL_Value:"+KLLDivideValues.getDouble(i)+"\tDD_Value:"+DDDivideValues.getDouble(i)+"\tTD_Value:"+TDDivideValues.getDouble(i));
            }


            LongArrayList ExactDivideRanks = getRankList(query_a,ExactDividedValues);
            LongArrayList KLLDivideRanks = getRankList(query_a,KLLDivideValues);
            LongArrayList DDDivideRanks = getRankList(query_a,DDDivideValues);
            LongArrayList TDDivideRanks = getRankList(query_a,TDDivideValues);

//            for(int i=0;i<QUERY_PER_TEST;i++){
//                int queryL=random.nextInt(totN),queryR=random.nextInt(totN);
//                if(queryL>queryR){int tmp=queryL;queryL=queryR;queryR=tmp;}
//                avgExact+=getRangeError(query_a,queryL,queryR,ExactDivideRanks,ExactDividedValues,query_a[queryL],query_a[queryR]);
//                avgKLL+=getRangeError(query_a,queryL,queryR,KLLDivideRanks,KLLDivideValues,query_a[queryL],query_a[queryR]);
//                avgDD+=getRangeError(query_a,queryL,queryR,DDDivideRanks,DDDivideValues,query_a[queryL],query_a[queryR]);
//                avgTD+=getRangeError(query_a,queryL,queryR,TDDivideRanks,TDDivideValues,query_a[queryL],query_a[queryR]);
//            }

//            double QueryMinV=query_a[0],QueryMaxV=query_a[totN-1];
//            double QueryMinV=query_a[totN/200],QueryMaxV=query_a[totN/200*199];
            for(int i=0;i<QUERY_PER_TEST;i++){
                int per=bucketNum/4,percentile=random.nextInt(per);
                double QueryMinV=query_a[totN/per*percentile],QueryMaxV=query_a[totN/per*(percentile+1)-1];

                double valL=QueryMinV+random.nextDouble()*(QueryMaxV-QueryMinV),valR=QueryMinV+random.nextDouble()*(QueryMaxV-QueryMinV);
                if(valL>valR){double t=valL;valL=valR;valR=t;}
                double err;
                exacts.add(err = getRangeError(query_a,valL,valR,ExactDivideRanks,ExactDividedValues,DataMinV,DataMaxV)/QUERY_PER_TEST);avgExact+=err;
                KLLs.add(err=getRangeError(query_a,valL,valR,KLLDivideRanks,KLLDivideValues,DataMinV,DataMaxV)/QUERY_PER_TEST);avgKLL+=err;
                DDs.add(err=getRangeError(query_a,valL,valR,DDDivideRanks,DDDivideValues,DataMinV,DataMaxV)/QUERY_PER_TEST);avgDD+=err;
                TDs.add(err=getRangeError(query_a,valL,valR,TDDivideRanks,TDDivideValues,DataMinV,DataMaxV)/QUERY_PER_TEST);avgTD+=err;
            }
        }
        for(DoubleArrayList s:new DoubleArrayList[]{exacts,KLLs,DDs,TDs})s.sort(Double::compare);

        String line = "";
        line += "\t\t"+bucketNum+"\t"+avgExact/TEST_CASE+"\t"+avgKLL/TEST_CASE+"\t"+avgDD/TEST_CASE+"\t"+avgTD/TEST_CASE;
        RESULT.get(RESULT_LINE).append(line);
        System.out.println("\t\t"+line);
        for(DoubleArrayList s:new DoubleArrayList[]{exacts,KLLs,DDs,TDs})System.out.println(s);
        for(DoubleArrayList s:new DoubleArrayList[]{exacts,KLLs,DDs,TDs})System.out.println("midERR:\t"+s.getDouble(s.size()/2)+"\t\t99%ERR:\t"+s.getDouble(s.size()*99/100));
//        System.out.println(exacts);
//        System.out.println(KLLs);
//        System.out.println(DDs);
//        System.out.println(TDs);
    }



    public static void main(String[] args) throws IOException {
        long START_T = new Date().getTime();
        MainForCompareOurAndApprox main = new MainForCompareOurAndApprox();
        typeMem2Alpha = new Object2DoubleOpenHashMap<>();
        int[] VoltageDDAlphaPow95For64KBTo256KB=new int[]{135,140,144,147,150,153,155,157,159,161,163,164,166};
        for(int i=64,j=0;i<=256;i+=16,j++)typeMem2Alpha.put(IntIntPair.of(1,i),Math.pow(0.95,VoltageDDAlphaPow95For64KBTo256KB[j]));
        int[] PriceDDAlphaPow95For64KBTo256KB=new int[]{129,134,138,141,144,147,150,152,154,156,157,159,161};
        for(int i=64,j=0;i<=256;i+=16,j++)typeMem2Alpha.put(IntIntPair.of(2,i),Math.pow(0.95,PriceDDAlphaPow95For64KBTo256KB[j]));
        int[] PriceDDAlphaPow95For256KBTo1024KB=new int[]{146,150,154,158,161,163,166,168,170,171,173,175,176};
        for(int i=256,j=0;i<=1024;i+=64,j++)typeMem2Alpha.put(IntIntPair.of(3,i),Math.pow(0.95,PriceDDAlphaPow95For256KBTo1024KB[j]));
        type2MinV=new Int2DoubleOpenHashMap();
        type2MinV.put(0,1.0);// Pareto data
        type2MinV.put(1,-205.25);type2MinV.put(2,0);type2MinV.put(3,0);

//        prList.add(1.0);
//
//        for (int dataType = startType; dataType <= endType; dataType++) {
////            for(int i=0;i<10;i++)main.RESULT.add(new StringBuilder());
//            main.prepareA(dataType);
//            for (int bucket : new int[]{64, 128, 256, 512, 1024, 2048, 4096, 8192/**/}) {
//                main.RESULT.add(new StringBuilder());
//                main.testMaxPartition((int)2e7,1024*256,bucket); // 3e7_512kb:tooLargeCount
//                main.RESULT_LINE++;
//            }
//        }
//        String LineColumn="\t\tx\tavgExact\tavgKLL\tavgKLLDupli\tavgDD\tavgDSS\tavgTD\tavgGK"+"\t\t|";
//        System.out.println("\t"+LineColumn);
//        for(StringBuilder sb:main.RESULT)System.out.println("\t"+sb.toString());

//        int[] buckets = new int[]{64,/*128,256,*/512,/*1024,2048,*/4096,/**//*16384/**/};
//        int tmpLine=0,ALL_LINE=buckets.length+4;
//        for (int dataType = startType; dataType <= endType; dataType++) { // CHECK IT
//
//            for(int i=0;i<ALL_LINE;i++)main.RESULT.add(new StringBuilder());
//
//            main.prepareA(dataType);
//
//
//            int[] nn=new int[]{10000000,1000000,5000000,10000000};
//            int[][] mm=new int[][]{new int[]{1024*256},new int[]{1024*256},new int[]{1024*256},new int[]{1024*256},new int[]{1024*256},new int[]{1024*256}};
//            int[] ttcc=new int[]{1,1,1,1,1,1,1};
//            main.RESULT.get(tmpLine+0).append("\n------------------------------------------------\n\t|||DATASET:"+"\t"+dataType);
//            for(int i=1;i<ALL_LINE;i++)main.RESULT.get(tmpLine+i).append("\t|||\t");
//
//            for(int ii=0;ii<0;ii++)
//                for(int jj=0;jj<1;jj++){
//                    int queryN=nn[ii],queryMem=mm[ii][jj];
//                    TEST_CASE*=ttcc[ii];
//
//                    main.RESULT.get(tmpLine+0).append("\t\t\t\t\t\t\t\t\t\t\t");
//                    String LineCondition="\t\tN:\t"+queryN+"\tMemory:\t"+queryMem/1024+"KB\t";
//                    String LineColumn="\t\tx\tavgExact\tavgKLL\tavgDD\tavgTD"+"\t\t|";
//                    main.RESULT.get(tmpLine+1).append(LineCondition);
//                    main.RESULT.get(tmpLine+2).append(LineColumn);
//                    System.out.println("\t\t"+LineCondition);
//                    System.out.println("\t\t"+LineColumn);
////                    System.out.println("show Var Of FilterSize!\tTEST_CASE=" + TEST_CASE + "\tDATASET:" + dataType + "\tqueryN:\t" + queryN + "\tmemory:\t" + queryMem);
////                    System.out.println("\t\tPr\tavgFilterSize\tvarFilterSize\testiFS\tsum_mw\tsum_0.5mw^2\t\tdetail");
//
//                    main.RESULT_LINE=tmpLine+3;
////                        main.testMapReduce(queryN, queryMem, bucketNum);
////                    for(int n:new int[]{100000,500000,1000000,5000000,10000000}) {
//                        main.testSquarePartition(queryN);
//                        main.RESULT_LINE++;//break;
////                    }
//                    TEST_CASE/=ttcc[ii];
//                }
//
//            tmpLine+=ALL_LINE;
////            System.out.println("\n-------------------------\n");
//        }
//
//        for(StringBuilder sb:main.RESULT)System.out.println("\t"+sb.toString());


        main.prepareA(2);
        for(int queryN:new int[]{1000000,2000000,3000000,4000000,5000000,6000000,7000000,8000000,9000000,10000000}) {
            main.RESULT.add(new StringBuilder());
            main.testSquarePartition(queryN);
            main.RESULT_LINE++;
        }
        String LineColumn="\t\t\tx\tavgExact\tavgKLL\tavgKLLDupli\tavgDD\tavgDSS\tavgTD\tavgGK"+"\t\t|";
        System.out.println("\t"+LineColumn);
        for(StringBuilder sb:main.RESULT)System.out.println("\t"+sb.toString());

        System.out.println("\t\t\tALL_TIME:" + (new Date().getTime() - START_T));
    }
}