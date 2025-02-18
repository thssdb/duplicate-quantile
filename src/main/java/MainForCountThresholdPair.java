import it.unimi.dsi.fastutil.doubles.Double2IntMap;
import it.unimi.dsi.fastutil.doubles.Double2IntOpenHashMap;
import it.unimi.dsi.fastutil.doubles.DoubleArrayList;
import it.unimi.dsi.fastutil.doubles.DoubleIntPair;
import it.unimi.dsi.fastutil.longs.LongArrayList;
import it.unimi.dsi.fastutil.objects.ObjectArrayList;
import it.unimi.dsi.util.XoRoShiRo128PlusPlusRandomGenerator;
import it.unimi.dsi.util.XoRoShiRo128PlusRandom;
import it.unimi.dsi.util.XorShift1024StarPhiRandomGenerator;
import org.apache.commons.math3.distribution.ZipfDistribution;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.text.DecimalFormat;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Comparator;
import java.util.Date;
import java.util.Random;

public class MainForCountThresholdPair {
    int dataType;
    static int startType = 4, endType = 4;
    static int[] Ns=new int[]{(int)6e7,(int)3e7,(int)3e7,(int)1.1e8,(int)3e7};
    static int N; // CHECK IT
    public static int TEST_CASE = 4*8*9; // CHECK IT
    static double[] a;
    static ArrayList<String> result_strings = new ArrayList<>();
    int RESULT_LINE = 0;

    public void prepareA(int dataType) throws IOException {
        N=Ns[dataType];
        if (a == null||a.length<N) a = new double[N];
        this.dataType = dataType;
        if (dataType == 4) {
//            LogNormalDistribution log21=new LogNormalDistribution(new XoRoShiRo128PlusPlusRandomGenerator(233),1,2);
//            for (int i = 0; i < Ns[dataType]; i++) a[i] = log21.sample();
            XoRoShiRo128PlusPlusRandomGenerator random=new XoRoShiRo128PlusPlusRandomGenerator(233);
            for (int i = 0; i < Ns[dataType]; i++) a[i] =
                Math.pow(-1, random.nextInt(2)) * Math.pow(10.0, (2 * Math.pow(random.nextDouble(), 8) - 1) * 300); // D_hard
            return;
        }
        BufferedReader reader = null;
        if (dataType == 0)reader = new BufferedReader(new FileReader(new File("Zipf3E7Alpha10.txt")));
        if (dataType == 1)reader = new BufferedReader(new FileReader(new File("DupliTorqueVoltage.txt")));
        if (dataType == 2)reader = new BufferedReader(new FileReader(new File("DupliECommercePrice.txt")));
        if (dataType == 3)reader = new BufferedReader(new FileReader(new File("DupliCustom.txt")));
        reader.readLine(); // ignore first line.
        String line;
        int cntN = 0;
        while ((line = reader.readLine()) != null) {
            a[cntN++] = Double.parseDouble(line);
            if (cntN == N) break;
        }
    }

    public void prepareZipf(int dataType,double alpha) throws IOException {
        N=Ns[dataType];
        if (a == null) a = new double[N];
        this.dataType = dataType;
        ZipfDistribution dis = new ZipfDistribution(new XorShift1024StarPhiRandomGenerator(233), N/100, alpha);
        int[] v2v=new int[N+1];
        XoRoShiRo128PlusRandom random=new XoRoShiRo128PlusRandom(233);
        for(int i=1;i<=N;i++){
            v2v[i]=i;
            if(i>1){
                int p=random.nextInt(i);
                v2v[i]=v2v[p];
                v2v[p]=i;
            }
        }
        for (int i = 0; i < N; i++)
            a[i] = v2v[dis.sample()];
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

    static private long dataToLong(double data) {
        long result = Double.doubleToLongBits((double) data);
        return data >= 0d ? result : result ^ Long.MAX_VALUE;
    }

    private double longToResult(long result) {
        result = (result >>> 63) == 0 ? result : result ^ Long.MAX_VALUE;
        return Double.longBitsToDouble(result);
    }

    DecimalFormat fnum = new DecimalFormat("#0.00");

    public int getPos(LongArrayList list,long v){
        int L=0,R=list.size()-1;
        while(L<R){
            int mid=(L+R)/2;
            if(list.getLong(mid)>=v)R=mid;
            else L=mid+1;
        }return L;
    }

    public double getApproxRank(long dataL,KLLDupliPair sketch,ObjectArrayList<DoubleIntPair> orderedHeavyHitters){
        double dataV=longToResult(dataL);
        long rank=sketch.getApproxRank(dataL);
        for(DoubleIntPair x:orderedHeavyHitters)if(x.leftDouble()<=dataV)rank+=x.rightInt();
        return rank;
    }
    public double findMinValueWithRank(long K,KLLDupliPair sketch,ObjectArrayList<DoubleIntPair> orderedHeavyHitters){
        long L=Long.MIN_VALUE,R=Long.MAX_VALUE,mid;
        while(L<R){
            mid = L + ((R - L) >>>1);
//      System.out.println("2fen\t\t"+L+"..."+R+"\t\tmid="+mid+"\t\t"+(getApproxRank(mid)>=K));
            if(getApproxRank(mid,sketch,orderedHeavyHitters)>K)R=mid;
            else L=mid+1;
        }
//    System.out.println("FT K:"+K+"\tN:"+getN()+" rank(L):"+getApproxRank(L));
        return longToResult(L);
    }
    public LongArrayList getApproxRanks(DoubleArrayList valueList,KLLDupliPair sketch,ObjectArrayList<DoubleIntPair> heavyHitters){
        int cntRankInHH=0,nextPosInHH=0;
        LongArrayList result=new LongArrayList();
        LongArrayList sketchRankList=sketch.getApproxRanks(valueList);
        for(int i=0;i<valueList.size();i++){
            double value=valueList.getDouble(i);
            while(nextPosInHH<heavyHitters.size()&&value>=heavyHitters.get(nextPosInHH).firstDouble()){
                cntRankInHH+=heavyHitters.get(nextPosInHH).secondInt();
                nextPosInHH++;
            }
            result.add(cntRankInHH+sketchRankList.getLong(i));
        }
//        System.out.println("\n\t\tvalueList:"+valueList+"\n\t\tapproxRanks:"+result+"\n");
        return result;
    }
    public double getDuplicateEstimationError(DoubleArrayList valueList,LongArrayList approxRanks){
        double ans=0;
        int queryN=approxRanks.size(),cntValueRank=0;
        for(int i=0;i<queryN;i++){
            int j=i+1;
            while(j<queryN&&valueList.getDouble(j)==valueList.getDouble(i))j++;
            j-=1;
            ans+=1.0*(j-i+1)*Math.abs(j+1-approxRanks.getLong(i))/queryN;
            i=j;
        }
        return ans/queryN;
    }

    public void testErrorWithThreshold(int queryN, int queryByte,double threshold) throws IOException {
        double[] query_a = new double[queryN];

        int[] LL = new int[TEST_CASE];
        int[] RR = new int[TEST_CASE];
        Random random = new Random(233);
        for (int i = 0; i < TEST_CASE; i++) {
            LL[i] = random.nextInt(N - queryN + 1);
            RR[i] = LL[i] + queryN;
        }

        double avgErr=0,avgSketchNRate=0;
        final double inf=1e6;
        Double2IntOpenHashMap v2c=new Double2IntOpenHashMap();
        for (int T = 0; T < TEST_CASE; T++) {
            int L = LL[T], R = RR[T];
            System.arraycopy(a, L, query_a, 0, R - L);
            Arrays.parallelSort(query_a);

            v2c.clear();
            int tmpC=0;
            for (int i = 0; i < R-L; i++){
                if(i>0&&query_a[i]!=query_a[i-1]){
                    v2c.put(query_a[i-1],tmpC);
                    tmpC=0;
                }
                else tmpC+=1;
            }
            v2c.put(query_a[R-L-1],tmpC);

            ObjectArrayList<DoubleIntPair> orderedHeavyHitters=new ObjectArrayList<>();
            int thresholdCount=(int)(queryN*threshold);
            for(Double2IntMap.Entry entry:v2c.double2IntEntrySet())
                if(entry.getIntValue()>thresholdCount){
                    orderedHeavyHitters.add(DoubleIntPair.of(entry.getDoubleKey(),entry.getIntValue()));
                }
            orderedHeavyHitters.sort(Comparator.comparingDouble(DoubleIntPair::keyDouble));
            long restByte=queryByte-orderedHeavyHitters.size()*8*2L;
            if(restByte<=Math.log(queryN)*8) {
                avgErr += inf / TEST_CASE;
                continue;
            }
            KLLDupliPair worker = new KLLDupliPair(queryByte-orderedHeavyHitters.size()*8*2);

            for (int i = L; i < R; i++)
                if(v2c.get(a[i])<=thresholdCount)
                    worker.update(dataToLong(a[i]));

            System.out.println("\t\torderedHeavyHitters.size():\t"+orderedHeavyHitters.size());
            worker.show();

            double err=0;

//            double q_add = 1e-4, q_start = q_add, q_end = 1-q_add, q_count = Math.floor((q_end - q_start - 1e-10) / q_add) + 1;
//            for (double q = q_start; q < q_end + 1e-10; q += q_add) {
//                int query_rank = (int) (q * queryN);
////                double full_v = longToResult(longVs.getLong(tmpQID++));
////                while(tmpHeavyHitterID<orderedHeavyHitters.size()&&orderedHeavyHitters.get(tmpHeavyHitterID).leftDouble()<=query_a[query_rank]) {
////                    tmpHeavyHitterRank += orderedHeavyHitters.get(tmpHeavyHitterID).rightInt();
////                    tmpHeavyHitterID++;
////                }
////                double full_v = longToResult(worker.findMinValueWithRank(query_rank-tmpHeavyHitterRank));
//                double full_v = findMinValueWithRank(query_rank,worker,orderedHeavyHitters);
//                int full_delta_rank = getDeltaRank(query_a, queryN, full_v, query_rank);
//                double full_relative_err = 1.0 * full_delta_rank / (queryN);
//                err += Math.abs(full_relative_err) / (q_count);
////                System.out.println("\t\t\tq:"+q+"\tquery_rank:"+query_rank+"\tanswer:"+query_a[query_rank]+"\ttmpHeavyHitterRank:"+tmpHeavyHitterRank+"\tfull_v:"+full_v+"\t\testi_rank_of_fullV:"+worker.getApproxRank(dataToLong(full_v)));
////                System.out.println("\t\t\tq:"+q+"\tquery_rank:"+query_rank+"\tanswer:"+query_a[query_rank]+"\tfull_v:"+full_v+"\tfull_relative_err:"+full_relative_err);
//            } // original error


            DoubleArrayList queriedDataV=new DoubleArrayList(query_a);
            err=getDuplicateEstimationError(queriedDataV,getApproxRanks(queriedDataV,worker,orderedHeavyHitters));
              //new error

            avgErr+=err/TEST_CASE;
            avgSketchNRate+=1.0*worker.getN()/queryN/TEST_CASE;
            System.out.println("\t\tthreshold:\t"+threshold+"\torderedHeavyHitters.size():\t"+orderedHeavyHitters.size()+"\tsketch.N:"+worker.getN()+"\terr:\t"+err);
        }
//        System.out.println("N:\t"+queryN+"\tM:\t"+queryByte+"\t\tavgDupli:\t"+avgDupliInSketch/TEST_CASE);
        String prefix= "N:\t"+queryN+"\tMemByte:\t"+queryByte+ "\tHeavyHitterThreshold:\t"+threshold+ "\tavgSketchNRate:\t"+avgSketchNRate+"\tavgErr:";
        String content="\t"+avgErr;
        if(RESULT_LINE>=result_strings.size())result_strings.add(prefix+content);
        else result_strings.set(RESULT_LINE,result_strings.get(RESULT_LINE)+content);
        System.out.println(prefix+content+"\n");
    }






    public static void main(String[] args) throws IOException {
        long START = new Date().getTime();
        System.out.println("MainForCountThresholdPair\nTEST_CASE=" + TEST_CASE);

        MainForCountThresholdPair main;
//
        result_strings=new ArrayList<>();
        DoubleArrayList HeavyHitterThreshold=new DoubleArrayList();
        for(double eps=1e-4,i=eps;i<1;i*=1.666)HeavyHitterThreshold.add(i);
        System.out.println("\t\t\tHeavyHitterThreshold:\t"+HeavyHitterThreshold);
        for (int dataType:new int[]{0}) { // CHECK IT
            System.out.println("DATASET:"+dataType);
            main = new MainForCountThresholdPair();
            for(double zipfAlpha:new double[]{0.8,1.1}) {

                main.prepareZipf(dataType, zipfAlpha);
                for (int queryN : new int[]{(int) 1e7})
                    for (int queryByte : new int[]{1024*8,1024 * 32}) {
                        result_strings.add("\nzipfAlpha:\t"+zipfAlpha+"\t\tMemByte:\t"+queryByte);
                        main.RESULT_LINE++;
                        for (double threshold : HeavyHitterThreshold) {
                            main.testErrorWithThreshold(queryN, queryByte, threshold);
                            main.RESULT_LINE++;
                        }
                    }
            }
            System.out.println();
        }
        for (String s : result_strings)System.out.println(s);
        System.out.println("\t\tALL_TIME:" + (new Date().getTime() - START));
    }
}
