import it.unimi.dsi.fastutil.longs.LongArrayList;
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

public class MainForKLLDupliFT {
    int dataType;
    static int startType = 1, endType = 2;
    static int pageN = 8192;
    static int N = (int)3e7, pageNum = N / pageN; // CHECK IT
    public static int TEST_CASE = 10*3*5; // CHECK IT
    static double[] a;
    static ArrayList<String> err_result = new ArrayList<>();
    static ArrayList<String> time_result = new ArrayList<>();
    int RESULT_LINE = 0;

    Random random = new Random(233);

    public void prepareA(int dataType) throws IOException {
        if (a == null) a = new double[N];
        this.dataType = dataType;
        BufferedReader reader = null;
        if (dataType == 1)reader = new BufferedReader(new FileReader(new File("DupliTorqueVoltage.txt")));
        if (dataType == 2)reader = new BufferedReader(new FileReader(new File("DupliECommercePrice.txt")));
        if (dataType == 3)reader = new BufferedReader(new FileReader(new File("DupliElectric.txt")));
        reader.readLine(); // ignore first line.
        String line;
        int cntN = 0;
        while ((line = reader.readLine()) != null) {
            a[cntN++] = Double.parseDouble(line);
            if (cntN == N) break;
        }
    }
    public void prepareZipf(int dataType,double alpha) throws IOException {
        if (a == null) a = new double[N];
        this.dataType = dataType;

        ZipfDistribution dis = new ZipfDistribution(new XorShift1024StarPhiRandomGenerator(233), N/233, alpha);
        for (int i = 0; i < N; i++)
            a[i] = dis.sample();
    }

    public void prepareUniform(int dataType,int repeat,boolean shuffle) throws IOException {
        if (a == null) a = new double[N];
        this.dataType = dataType;

//        UniformIntegerDistribution dis = new UniformIntegerDistribution(new XorShift1024StarPhiRandomGenerator(233),0,N/repeat);
//        for (int i = 0; i < N; i++)
//            a[i] = dis.sample();
        for(int i=0;i<N;i++)a[i]=(double)(i/repeat);
        if(shuffle) {
            for (int i = 1; i < N; i++) {
                int p = random.nextInt(i);
                double aa = a[p];
                a[p] = a[i];
                a[i] = aa;
            }
        }else{
            for (int i = repeat; i < N; i++) if(i%repeat==0){
                int p = random.nextInt(i/repeat);
                double aa = a[p*repeat];
                for(int j=p*repeat;j<(p+1)*repeat;j++)a[j]=a[i];
                for(int j=i;j<i+repeat&&j<N;j++)a[j]=aa;
//                a[p] = a[i];
//                a[i] = aa;
            }
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

    DecimalFormat fnum = new DecimalFormat("#0.00");

    public int getPos(LongArrayList list,long v){
        int L=0,R=list.size()-1;
        while(L<R){
            int mid=(L+R)/2;
            if(list.getLong(mid)>=v)R=mid;
            else L=mid+1;
        }return L;
    }


    public void testError(int queryN, int queryByte) throws IOException {
        int sketchM = Integer.highestOneBit(queryByte / (8 + 8 + 2));
        long full_time = 0, merge_time = 0;
        double err_full = 0, err_merge = 0;
        double MMP_full = 0;
        double[] query_a = new double[queryN];

        int[] LL = new int[TEST_CASE];
        int[] RR = new int[TEST_CASE];
        Random random = new Random(233);
        for (int i = 0; i < TEST_CASE; i++) {
            LL[i] = random.nextInt(N - queryN + 1);
            RR[i] = LL[i] + queryN;
        }
        long errBound99=0;
        long PairT0T0=0,PairFF=0,PairTF=0,PairT1T2=0;
        double avgDupliInSketch=0;
        for (int T = 0; T < TEST_CASE; T++) {
            int L = LL[T], R = RR[T];

            full_time -= new Date().getTime();
            KLLDupliFT worker = new KLLDupliFT(queryByte);
            for (int i = L; i < R; i++)
                worker.update(dataToLong(a[i]));
            PairT0T0+=worker.PairT0T0;
            PairFF+=worker.PairFF;
            PairTF+=worker.PairTF;
            PairT1T2+= worker.PairT1T2;
            errBound99=worker.queryRankErrBound(0.99);
            if(T==0)worker.showCompact();
//            if(T==0) {
//                worker.show();
////                worker.showNum();
//            }
//            System.out.println("\t\t|worker|:\t"+worker.getMaximumMapCapacity()+"\t\tsketchM:"+sketchM);
            full_time += new Date().getTime();
            avgDupliInSketch+=worker.getAvgDupliInSketch();

            if (R - L >= 0) System.arraycopy(a, L, query_a, 0, R - L);
            Arrays.sort(query_a);

            double q_add = 0.0001, q_start = q_add, q_end = 1-q_add, q_count = Math.floor((q_end - q_start - 1e-10) / q_add) + 1;
            for (double q = q_start; q < q_end + 1e-10; q += q_add) {
                int query_rank = (int) (q * queryN);
                double full_v = longToResult(worker.findMinValueWithRank(query_rank));
                int full_delta_rank = getDeltaRank(query_a, queryN, full_v, query_rank);
                double full_relative_err = 1.0 * full_delta_rank / (queryN);
                err_full += Math.abs(full_relative_err) / (q_count * TEST_CASE);
            }
        }
//        System.out.println("N:\t"+queryN+"\tM:\t"+queryByte+"\t\tavgDupli:\t"+avgDupliInSketch/TEST_CASE);
        System.out.println("\t\t" + queryN +/*"\t"+err_merge+*/"\t" + err_full+"\t99%bound:\t"+1.0*errBound99/queryN+"\tPairs(TT,FF,TF,T1T2):\t"+PairT0T0+"\t"+PairFF+"\t"+PairTF+"\t"+PairT1T2);

        err_result.set(RESULT_LINE, err_result.get(RESULT_LINE).concat("\t\t\t" + err_full+"\t"+1.0*errBound99/queryN));
        time_result.set(RESULT_LINE, time_result.get(RESULT_LINE).concat("\t\t\t" + full_time));

    }

    public void show_time_result() {
        System.out.println(time_result);
    }

    public static void setTestCase(int tc) {
        TEST_CASE = tc;
    }


    public static void main(String[] args) throws IOException {
        long START = new Date().getTime();
        System.out.println("KLLDupliFT\nTEST_CASE=" + TEST_CASE);
//        KLLDupliFT sketchFT = new KLLDupliFT(8*40);
//        KLLSketchLazyExactPriori sketch=new KLLSketchLazyExactPriori(8*40);
//        for(int i=0;i<120;i++){
//            sketchFT.update(i/15);
//            sketch.update(i/15);
//        }
//        sketchFT.showNum();
//        sketch.showNum();
//        for(int i=0;i<8;i++)System.out.println("\t\t"+sketchFT.findMinValueWithRank(i*15)+"\t\t"+sketch.findMinValueWithRank(i*15));
//        System.out.println("\t\t"+sketchFT.getApproxRank(6)+"\t"+sketch.getApproxRank(6));

        MainForKLLDupliFT main;

        for (int dataType = 0; dataType <= 0; dataType++) { // CHECK IT
            main = new MainForKLLDupliFT();
            for (int queryN : new int[]{5000000})
                for (int queryByte : new int[]{1024*64})
                    for (double alpha : new double[]{0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0,1.1,1.2/*,1.3,1.4,1.5,1.6,1.7,1.8,1.9,2.0,2.1,2.2,2.3,2.4*/}) {
                        if (dataType == 0) {
                            err_result.add("N:" + queryN + ", " + "queryByte:" + queryByte + ", " + "alpha:" + alpha + "\t");
                            time_result.add("N:" + queryN + ", " + "queryByte:" + queryByte + ", " + "alpha:" + alpha + "\t");
                        }
                        main.prepareZipf(dataType,alpha);
                        main.testError(queryN, queryByte);
                        main.RESULT_LINE++;
                    }
        }
        System.out.println("\nKLLDupliFT Zipf\tTEST_CASE=" + TEST_CASE);
        System.out.println("Error rate:");
        for (String s : err_result)
            System.out.println(s);
//
//
        err_result=new ArrayList<>();
        time_result=new ArrayList<>();
        for (int dataType = 0; dataType <= 0; dataType++) { // CHECK IT
            main = new MainForKLLDupliFT();
            for (int queryN : new int[]{5000000})
                for (int queryByte : new int[]{1024*64})
                    for (int repeat : new int[]{1,2,5,10,20,50,100,200,500,1000,2000,5000,10000}) {
                        if (dataType == 0) {
                            err_result.add("N:" + queryN + ", " + "queryByte:" + queryByte + ", " + "repeat:" + repeat + "\t");
                            time_result.add("N:" + queryN + ", " + "queryByte:" + queryByte + ", " + "repeat:" + repeat + "\t");
                        }
                        main.prepareUniform(dataType,repeat,false);
                        main.testError(queryN, queryByte);
                        main.RESULT_LINE++;
                    }
        }
        System.out.println("\nKLLDupliFT Uniform\tTEST_CASE=" + TEST_CASE);
        System.out.println("Error rate:");
        for (String s : err_result)
            System.out.println(s);


//        err_result=new ArrayList<>();
//        time_result=new ArrayList<>();
//        for (int dataType = startType; dataType <= endType; dataType++) { // CHECK IT
//            System.out.println("DATASET:"+dataType);
//            main = new MainForKLLDupliFT();
//            for (int queryN : new int[]{5000000})
//                for (int queryByte : new int[]{1024*16,1024*32,1024*64,1024*128,1024*256/*,1024*512,1024*1024/*,1024*2048*/}){
//                    if (dataType == startType) {
//                        err_result.add("N:" + queryN + ", " + "queryByte:\t" + queryByte+ "\t");
//                        time_result.add("N:" + queryN + ", " + "queryByte:\t" + queryByte+ "\t");
//                    }
//                    main.prepareA(dataType);
//                    main.testError(queryN, queryByte);
//                    main.RESULT_LINE++;
//                }
//            System.out.println();
//        }
//        System.out.println("\nError rate:");
//        for (String s : err_result)System.out.println(s);

        err_result=new ArrayList<>();
        time_result=new ArrayList<>();
        for (int dataType = startType; dataType <= endType; dataType++) { // CHECK IT
            System.out.println("DATASET:"+dataType);
            main = new MainForKLLDupliFT();
            for (int queryN : new int[]{2000000,4000000,6000000,8000000,10000000,15000000,20000000/**/})
                for (int queryByte : new int[]{1024*256}){
                    if (dataType == startType) {
                        err_result.add("N:" + queryN + ", " + "queryByte:" + queryByte+ "\t");
                        time_result.add("N:" + queryN + ", " + "queryByte:" + queryByte+ "\t");
                    }
                    main.prepareA(dataType);
                    main.testError(queryN, queryByte);
                    main.RESULT_LINE++;
                }
            System.out.println();
        }
        System.out.println("\nError rate:");
        for (String s : err_result)System.out.println(s);

        System.out.println("\t\tALL_TIME:" + (new Date().getTime() - START));
    }
}
