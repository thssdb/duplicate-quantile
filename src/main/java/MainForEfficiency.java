import it.unimi.dsi.fastutil.ints.IntArrayList;
import it.unimi.dsi.fastutil.longs.LongArrayList;
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
import java.util.Date;
import java.util.Random;

public class MainForEfficiency {
    int dataType;
    static int startType = 4, endType = 4;
    static int[] Ns=new int[]{(int)6e7,(int)3e7,(int)3e7,(int)1.1e8,(int)3e7,(int)3e7};
    static int N; // CHECK IT
    public static int TEST_CASE = 1; // CHECK IT
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
        if (dataType == 5) {
            XoRoShiRo128PlusPlusRandomGenerator random=new XoRoShiRo128PlusPlusRandomGenerator(233);
            for (int i = 0; i < Ns[dataType]; i++)
                a[i] =-i;//random.nextDoubleFast();
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


    public void test(int queryN, int queryByte) throws IOException {
        double[] query_a = new double[queryN];

        int[] LL = new int[TEST_CASE];
        int[] RR = new int[TEST_CASE];
        Random random = new Random(233);
        for (int i = 0; i < TEST_CASE; i++) {
            LL[i] = random.nextInt(N - queryN + 1);
            RR[i] = LL[i] + queryN;
        }

        double pairTime=0,pairFastTime=0,pairFasterTime=0,KLLTime=0;
        double pairErr=0,fastErr=0,fasterErr=0,KLLErr=0;
        for (int T = 0; T < TEST_CASE; T++) {
            int L = LL[T], R = RR[T];
            System.arraycopy(a, L, query_a, 0, R - L);
            Arrays.parallelSort(query_a);

            KLLDupliPair pairWorker = new KLLDupliPair(queryByte);
            KLLDupliPairFast pairFastWorker = new KLLDupliPairFast(queryByte);
            KLLDupliPairFaster pairFasterWorker = new KLLDupliPairFaster(queryByte);
            KLLSketchLazyExactPriori KLLWorker = new KLLSketchLazyExactPriori(queryByte);
            pairTime-=new Date().getTime();
            for (int i = L; i < R; i++)
                pairWorker.update(dataToLong(a[i]));
            pairTime+=new Date().getTime();
            if(T==0)pairWorker.showCompact();

            pairFastTime-=new Date().getTime();
            for (int i = L; i < R; i++)
                pairFastWorker.update(dataToLong(a[i]));
            pairFastTime+=new Date().getTime();
            if(T==0)pairFastWorker.showCompact();
//            pairFastWorker.showNum();

            pairFasterTime-=new Date().getTime();
            for (int i = L; i < R; i++)
                pairFasterWorker.update(dataToLong(a[i]));
            pairFasterTime+=new Date().getTime();
            if(T==0)pairFasterWorker.showCompact();
//            pairFasterWorker.showNum();

            KLLTime-=new Date().getTime();
            for (int i = L; i < R; i++)
                KLLWorker.update(dataToLong(a[i]));
            KLLTime+=new Date().getTime();
            if(T==0)KLLWorker.showCompact();

            double q_add = 0.0001, q_start = q_add, q_end = 1-q_add, q_count = Math.floor((q_end - q_start - 1e-10) / q_add) + 1;
            LongArrayList Ks=new LongArrayList();
            for (double q = q_start; q < q_end + 1e-10; q += q_add)Ks.add((int) (q * queryN));
            LongArrayList pairVs=pairWorker.findMinValuesWithRanks(Ks),fastVs=pairFastWorker.findMinValuesWithRanks(Ks),fasterVs=pairFasterWorker.findMinValuesWithRanks(Ks);
            int tmpQID=0;
            for (double q = q_start; q < q_end + 1e-10; q += q_add) {
                int query_rank = (int) (q * queryN);
                pairErr += Math.abs(1.0*getDeltaRank(query_a, queryN, longToResult(pairVs.getLong(tmpQID)), query_rank))  / (queryN)/ (q_count * TEST_CASE);
                fastErr += Math.abs(1.0*getDeltaRank(query_a, queryN, longToResult(fastVs.getLong(tmpQID)), query_rank))  / (queryN)/ (q_count * TEST_CASE);
                fasterErr += Math.abs(1.0*getDeltaRank(query_a, queryN, longToResult(fasterVs.getLong(tmpQID)), query_rank))  / (queryN)/ (q_count * TEST_CASE);
                KLLErr += Math.abs(1.0*getDeltaRank(query_a, queryN, longToResult(KLLWorker.findMinValueWithRank(query_rank)), query_rank))  / (queryN)/ (q_count * TEST_CASE);
                tmpQID++;
            }
        }
//        System.out.println("N:\t"+queryN+"\tM:\t"+queryByte+"\t\tavgDupli:\t"+avgDupliInSketch/TEST_CASE);
        String prefix= "N:\t"+queryN+"\tMemByte:\t"+queryByte;
        String content="\t|||\tpairWorker_Err:\t"+pairErr+"\tpairFastWorker_Err:\t"+fastErr+"\tpairFasterWorker_Err:\t"+fasterErr+"\tKLL_Err:\t"+KLLErr+"\tKLLPair_Time:\t"+pairTime/TEST_CASE+"\tKLLPairFast_Time:\t"+pairFastTime/TEST_CASE+"\tKLLPairFaster_Time:\t"+pairFasterTime/TEST_CASE+"\tKLL_Time:\t"+KLLTime/TEST_CASE;
        if(RESULT_LINE>=result_strings.size())result_strings.add(prefix+content);
        else result_strings.set(RESULT_LINE,result_strings.get(RESULT_LINE)+content);
        System.out.println(prefix+content+"\n");
    }






    public static void main(String[] args) throws IOException {
        long START = new Date().getTime();
        System.out.println("MainForCountThresholdPair\nTEST_CASE=" + TEST_CASE);

        MainForEfficiency main;

        IntArrayList queryNList=new IntArrayList();
//        queryNList.add((int)2333);
        queryNList.add((int)2e7);
        result_strings=new ArrayList<>();
//        for (int dataType:new int[]{1}) { // CHECK IT
//            System.out.println("\n---------------------------\nDATASET:"+dataType);
//            main = new MainForEfficiency();
//            main.prepareA(dataType);
//            for (int queryN : queryNList)
//                for (int queryByte : new int[]{1024*64/*1024*16,1024*64,/*1024*256*/}){
//                    main.test(queryN, queryByte);
//                    main.RESULT_LINE++;
//                }
//            System.out.println();
//        }

        for (int dataType:new int[]{1}) { // CHECK IT
            System.out.println("\n---------------------------\nDATASET:"+dataType);
            main = new MainForEfficiency();
            for (double alpha:new double[]{0.01,0.3,0.6,0.8,1.0,1.1,1.2}){
                    main.prepareZipf(dataType,alpha);
                    main.test((int)2e7, 1024*128);
                    main.RESULT_LINE++;
                }
            System.out.println();
        }
        for (String s : result_strings)System.out.println(s);
        System.out.println("\t\tALL_TIME:" + (new Date().getTime() - START));
    }
}
