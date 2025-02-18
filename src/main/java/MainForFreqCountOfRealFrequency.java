import it.unimi.dsi.fastutil.longs.Long2LongOpenHashMap;
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

public class MainForFreqCountOfRealFrequency {
    int dataType;
    static int startType = 4, endType = 4;
    static int[] Ns=new int[]{(int)3.1e7,(int)3e7,(int)3e7,(int)1.1e8,(int)3e7,(int)3e7};
    static int N; // CHECK IT
    public static int TEST_CASE = 4000; // CHECK IT
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
        if (dataType == 0)reader = new BufferedReader(new FileReader(new File("Zipf3E7Alpha8.txt")));
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

    /**
     * compare the HeavyHitters' count in sketch and the real count
     */
    public void testFreqCount(int queryN, int queryByte) throws IOException {
        double[] query_a = new double[queryN];

        int[] LL = new int[TEST_CASE];
        int[] RR = new int[TEST_CASE];
        Random random = new Random(233);
        for (int i = 0; i < TEST_CASE; i++) {
            LL[i] = random.nextInt(N - queryN + 1);
            RR[i] = LL[i] + queryN;
        }

        double pairFTTime=0,pairFastTime=0,pairFasterTime=0,KLLTime=0;
        double pairFTErr=0,fastErr=0,fasterErr=0,KLLErr=0;
        Long2LongOpenHashMap realValueCount=new Long2LongOpenHashMap();
        ObjectArrayList<double[]> caseCountRelativeTrueCountRelativeFreqCountValue=new ObjectArrayList<>();
        double[] avgRelativeTrueCount=new double[10000],avgCount=new double[10000],avgRelativeHashMapCount=new double[10000],countCases=new double[10000];
        double[] hashMapNonZeroAvgRelaCount=new double[10000],hashMapNonZeroAvgCount=new double[10000],hashMapNonZeroCases=new double[10000];
        double avgFreqThreshold=0,avgRelativeFreqNum=0,avgRelativeFreqCountSum=0;
        int freqValueNumInAllCases=10000,freqValueNumInAnyCases=0;
        for (int T = 0; T < TEST_CASE; T++) {
            int L = LL[T], R = RR[T];
            System.arraycopy(a, L, query_a, 0, R - L);

            KLLDupliPairFT pairFTWorker = new KLLDupliPairFT(queryByte);
            pairFTTime-=new Date().getTime();
            for (int i = L; i < R; i++)
                pairFTWorker.update(dataToLong(a[i]));
            pairFTTime+=new Date().getTime();
            if(T==0){
                pairFTWorker.showCompact();
                System.out.println("Sketch.freqValueCount:"+pairFTWorker.freqValueCount);
            }

            realValueCount.clear();
            caseCountRelativeTrueCountRelativeFreqCountValue.clear();
            for(double v:query_a)
                realValueCount.put(dataToLong(v),realValueCount.getOrDefault(dataToLong(v),0)+1);
            for(long longV: realValueCount.keySet()){
                long realCount=realValueCount.get(longV),sketchTrueCount=pairFTWorker.getTrueCountInSketch(longV),sketchHashMapCount=pairFTWorker.getHashMapCountInSketch(longV);
                caseCountRelativeTrueCountRelativeFreqCountValue.add(new double[]{realCount,1.0*sketchTrueCount/realCount,1.0*sketchHashMapCount/realCount,longToResult(longV)});
            }
            caseCountRelativeTrueCountRelativeFreqCountValue.sort(Comparator.comparingDouble(x -> -x[0]));
            int caseFreqValueNum=0;
            for(int i=0;i<caseCountRelativeTrueCountRelativeFreqCountValue.size();i++)if(caseCountRelativeTrueCountRelativeFreqCountValue.get(i)[0]>=pairFTWorker.getFreqThreshold())
                caseFreqValueNum=i+1;
            caseFreqValueNum=Math.max(caseFreqValueNum,2);
//            caseCountRelativeTrueCountRelativeFreqCountValue.size(caseFreqValueNum);
            freqValueNumInAllCases=Math.min(freqValueNumInAllCases,caseFreqValueNum);
            freqValueNumInAnyCases=Math.max(freqValueNumInAnyCases,caseFreqValueNum);

            avgFreqThreshold+=1.0*pairFTWorker.getFreqThreshold()/TEST_CASE;

//            for(int i=0;i<10;i++)System.out.print(Arrays.toString(caseCountRelativeTrueCountRelativeFreqCountValue.get(i))+"\t\t");System.out.println();
//            for(double[] tmp:caseCountRelativeTrueCountRelativeFreqCountValue)System.out.println(Arrays.toString(tmp));
//            System.out.println("\t\t??caseFreqValueNum:"+caseFreqValueNum+"\tpairFTWorker.getFreqThreshold():"+pairFTWorker.getFreqThreshold());
//            pairFTWorker.showCompact();
            for(int i=0;i<caseFreqValueNum;i++){
                countCases[i]+=1;
                avgCount[i]+=caseCountRelativeTrueCountRelativeFreqCountValue.get(i)[0];
                avgRelativeTrueCount[i]+=caseCountRelativeTrueCountRelativeFreqCountValue.get(i)[1];
                avgRelativeHashMapCount[i]+=caseCountRelativeTrueCountRelativeFreqCountValue.get(i)[2];
                if(caseCountRelativeTrueCountRelativeFreqCountValue.get(i)[2]>0){
                    hashMapNonZeroCases[i]+=1;
                    hashMapNonZeroAvgCount[i]+=caseCountRelativeTrueCountRelativeFreqCountValue.get(i)[0];
                    hashMapNonZeroAvgRelaCount[i]+=caseCountRelativeTrueCountRelativeFreqCountValue.get(i)[2];
                }
            }
//            System.out.println("\tMostFreqValue:"+ Arrays.toString(caseCountRelativeTrueCountRelativeFreqCountValue.get(0)));
        }
        avgCount=Arrays.copyOf(avgCount,freqValueNumInAnyCases);
        avgRelativeTrueCount=Arrays.copyOf(avgRelativeTrueCount,freqValueNumInAnyCases);
        avgRelativeHashMapCount=Arrays.copyOf(avgRelativeHashMapCount,freqValueNumInAnyCases);
//        System.out.println("N:\t"+queryN+"\tM:\t"+queryByte+"\t\tavgDupli:\t"+avgDupliInSketch/TEST_CASE);
        String prefix= "N:\t"+queryN+"\tMemByte:\t"+queryByte+"\tTEST_CASE:\t"+TEST_CASE;
        int ValueNum=freqValueNumInAnyCases;//Math.min(avgRelativeTrueCount.size(),20);


        String content="\tdataset:\t"+dataType+ "\tMem:\t"+(queryByte/1024)+"KB\t\tTop\t"+ValueNum+"\tvalues\t\tavgFreqThreshold:\t"+avgFreqThreshold;
        while(ValueNum>10&&avgCount[ValueNum-1]<avgFreqThreshold)ValueNum--;
        content+="\n\tPointID:\t";
        for(int i=0;i<ValueNum;i++)content+=(i+1)+"\t";
        content+="\n\tavgCountOfMostFreq:\t";
        for(int i=0;i<ValueNum;i++)content+=fnum.format(avgCount[i]/countCases[i])+"\t";
        content+="\n\tavgRelativeTrueCount:\t";
        for(int i=0;i<ValueNum;i++)content+=fnum.format(avgRelativeTrueCount[i]/countCases[i])+"\t";
        content+="\n\tavgRelativeHashMapCount:\t";
        for(int i=0;i<ValueNum;i++)content+=fnum.format(avgRelativeHashMapCount[i]/countCases[i])+"\t";

        content+="\n\thashMapNonZeroRate:\t";
        for(int i=0;i<ValueNum;i++)content+=fnum.format(hashMapNonZeroCases[i]/countCases[i])+"\t";
        content+="\n\thashMapNonZeroAvgCount:\t";
        for(int i=0;i<ValueNum;i++)if(hashMapNonZeroCases[i]>0)content+=fnum.format(hashMapNonZeroAvgCount[i]/hashMapNonZeroCases[i])+"\t";else content+="\t";
        content+="\n\thashMapNonZeroAvgRelaCount:\t";
        for(int i=0;i<ValueNum;i++)if(hashMapNonZeroCases[i]>0)content+=(hashMapNonZeroAvgRelaCount[i]/hashMapNonZeroCases[i])+"\t";else content+="\t";

        if(RESULT_LINE>=result_strings.size())result_strings.add("\n\n"+prefix+"\n\t\t\t"+content);
        else result_strings.set(RESULT_LINE,result_strings.get(RESULT_LINE)+"\n\n\t\t\t"+content);
        System.out.println(prefix+content+"\n");
    }






    public static void main(String[] args) throws IOException {
        long START = new Date().getTime();
        System.out.println("MainForFreqCountOfRealFrequency\nTEST_CASE=" + TEST_CASE);

        MainForFreqCountOfRealFrequency main;

        result_strings=new ArrayList<>();
        for (int dataType:new int[]{1}) { // CHECK IT
            System.out.println("\n---------------------------\nDATASET:"+dataType);
            main = new MainForFreqCountOfRealFrequency();
            main.prepareA(dataType);
            for(int queryByte:new int[]{1024*64,1024*128}) {
                main.testFreqCount((int) 1e7, queryByte);
                main.RESULT_LINE++;
                System.out.println();
            }
        }
        for (String s : result_strings)System.out.println(s);

//        for (int dataType:new int[]{1}) { // CHECK IT
//            System.out.println("\n---------------------------\nDATASET:"+dataType);
//            main = new MainForFreqCountOfRealFrequency();
//            for (double alpha:new double[]{0.01,0.3,0.6,0.8,1.0,1.1,1.2}){
//                    main.prepareZipf(dataType,alpha);
//                    main.testFreqCount((int)2e7, 1024*128);
//                    main.RESULT_LINE++;
//                }
//            System.out.println();
//        }
//        for (String s : result_strings)System.out.println(s);


        System.out.println("\t\tALL_TIME:" + (new Date().getTime() - START));
    }
}
