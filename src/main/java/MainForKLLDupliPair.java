import it.unimi.dsi.fastutil.ints.IntArrayList;
import it.unimi.dsi.fastutil.longs.LongArrayList;
import it.unimi.dsi.util.XoRoShiRo128PlusPlusRandomGenerator;
import it.unimi.dsi.util.XoRoShiRo128PlusRandom;
import it.unimi.dsi.util.XorShift1024StarPhiRandomGenerator;
import org.apache.commons.math3.distribution.ParetoDistribution;
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

public class MainForKLLDupliPair {
    int dataType;
    static int startType = 4, endType = 4;
    static int[] Ns=new int[]{(int)3e7,(int)3e7,(int)3e7,(int)1.1e8,(int)3e7,(int)3.1e7};
    static int[] Ms=new int[]{0,128,128,512};
    static int N = (int)3e7; // CHECK IT
    public static int TEST_CASE = 20,TEST_CASE_M = 1; // CHECK IT
    static double[] a;
    static ArrayList<String> err_result = new ArrayList<>();
    static ArrayList<String> time_result = new ArrayList<>();
    int RESULT_LINE = 0;

    Random random = new Random(233);

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
        if (a == null) a = new double[N];
        this.dataType = dataType;
        ZipfDistribution dis = new ZipfDistribution(new XorShift1024StarPhiRandomGenerator(233), N/233, alpha);
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
    public void preparePareto(int dataType,double alpha) throws IOException {
        if (a == null) a = new double[N];
        this.dataType = dataType;
        ParetoDistribution dis = new ParetoDistribution(new XorShift1024StarPhiRandomGenerator(233), 1, alpha);
        for (int i = 0; i < N; i++)
            a[i] = dis.sample();
        double Epsilon=1+5e-4;
        for (int i = 0; i < N; i++)a[i]=Math.pow(Epsilon,Math.ceil(Math.log(a[i])/Math.log(Epsilon)));
        double[] b=Arrays.copyOf(a,N);
        Arrays.sort(b);
        int count=0;
        for(int i=1;i<N;i++)if(b[i]!=b[i-1])count++;
        System.out.println("\t[Pareto]\talpha:\t"+alpha+"\tEpsilon:\t"+Epsilon+"\t\tcount:\t"+count+"\t\tN:\t"+N);
//        System.out.println("\t\t10%:\t"+a[N/10]);
//        System.out.println("\t\t50%:\t"+a[N/10*5]);
//        System.out.println("\t\t75%:\t"+a[N/100*75]);
//        System.out.println("\t\t99%:\t"+a[N/100*99]);
    }
    public void prepareLognormal(String muS,double sigma) throws IOException {
        N=Ns[dataType];
        if (a == null||a.length<N) a = new double[N];
        BufferedReader reader = new BufferedReader(new FileReader("Lognormal3E7Mu"+muS+"Sigma"+(int)(sigma*10+0.1)+".txt"));
        String line;
        int cntN = 0;
        while ((line = reader.readLine()) != null) {
            a[cntN++] = Double.parseDouble(line);
            if (cntN == N) break;
        }
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
            KLLDupliPairFaster worker = new KLLDupliPairFaster(queryByte);
            for (int i = L; i < R; i++)
                worker.update(dataToLong(a[i]));
            errBound99=worker.queryRankErrBound(0.99);
            if(T==0)worker.showCompact();
            if(T==0) {
                worker.show();
//                worker.showNum();
            }
//            System.out.println("\t\t|worker|:\t"+worker.getMaximumMapCapacity()+"\t\tsketchM:"+sketchM);
            full_time += new Date().getTime();

            if (R - L >= 0) System.arraycopy(a, L, query_a, 0, R - L);
            Arrays.parallelSort(query_a);

            double q_add = 0.0001, q_start = q_add, q_end = 1-q_add, q_count = Math.floor((q_end - q_start - 1e-10) / q_add) + 1;
            LongArrayList Ks=new LongArrayList();
            for (double q = q_start; q < q_end + 1e-10; q += q_add)Ks.add((int) (q * queryN));
            LongArrayList longVs=worker.findMinValuesWithRanks(Ks);
            int tmpQID=0,tmpMaxRankErr=0;
            for (double q = q_start; q < q_end + 1e-10; q += q_add) {
                int query_rank = (int) (q * queryN);
                double full_v = longToResult(longVs.getLong(tmpQID++));
                int full_delta_rank = getDeltaRank(query_a, queryN, full_v, query_rank);
                double full_relative_err = 1.0 * full_delta_rank / (queryN);
                err_full += Math.abs(full_relative_err) / (q_count * TEST_CASE);
                tmpMaxRankErr=Math.max(tmpMaxRankErr,full_delta_rank);
            }
//            System.out.println("\t\ttmpMaxRankErr:\t"+tmpMaxRankErr);
        }
//        System.out.println("N:\t"+queryN+"\tM:\t"+queryByte+"\t\tavgDupli:\t"+avgDupliInSketch/TEST_CASE);
        System.out.println("\t\t" + queryN +/*"\t"+err_merge+*/"\t" + err_full+"\t99%bound:\t"+1.0*errBound99/queryN);

        err_result.set(RESULT_LINE, err_result.get(RESULT_LINE).concat("\t\t\t" + err_full+"\t"+1.0*errBound99/queryN));
        time_result.set(RESULT_LINE, time_result.get(RESULT_LINE).concat("\t\t\t" + full_time));
    }





    public void testM(int queryN, double targetAE) throws IOException {
        int[] LL = new int[TEST_CASE_M];
        int[] RR = new int[TEST_CASE_M];
        Random random = new Random(233);
        for (int i = 0; i < TEST_CASE_M; i++) {
            LL[i] = random.nextInt(N - queryN + 1);
            RR[i] = LL[i] + queryN;
        }
        double avgM=0,avgH=0,avgTopCap=0,avgSig2=0;
        for (int T = 0; T < TEST_CASE_M; T++) {
            int ML = 1024 * 2, MR = 1024 * 1024;
            double tmpH=0,tmpTopCap=0,tmpSig2=0;
            while (ML * 1.01 < MR) {
                int queryByte = (ML + MR) / 2;
                double[] query_a = new double[queryN];
                double err_full = 0;
                int L = LL[T], R = RR[T];
                KLLDupliPair worker = new KLLDupliPair(queryByte);
                for (int i = L; i < R; i++)
                    worker.update(dataToLong(a[i]));
                if (R - L >= 0) System.arraycopy(a, L, query_a, 0, R - L);
                Arrays.parallelSort(query_a);

                double q_add = 1e-5, q_start = q_add, q_end = 1-q_add, q_count = Math.floor((q_end - q_start - 1e-10) / q_add) + 1;
                LongArrayList Ks=new LongArrayList();
                for (double q = q_start; q < q_end + 1e-10; q += q_add)Ks.add((int) (q * queryN));
                LongArrayList longVs=worker.findMinValuesWithRanks(Ks);
                int tmpQID=0,tmpMaxRankErr=0;
                for (double q = q_start; q < q_end + 1e-10; q += q_add) {
                    int query_rank = (int) (q * queryN);
                    double full_v = longToResult(longVs.getLong(tmpQID++));
                    int full_delta_rank = getDeltaRank(query_a, queryN, full_v, query_rank);
                    double full_relative_err = 1.0 * full_delta_rank / (queryN);
                    err_full += Math.abs(full_relative_err) / (q_count);
                    tmpMaxRankErr=Math.max(tmpMaxRankErr,full_delta_rank);
                }
//                System.out.println("\t\tqueryByte:" + queryByte + "\tavgERR:" + err_full);
                if (err_full <= targetAE) {
                    MR = queryByte;
                    tmpH=worker.cntLevel-1;
                    tmpTopCap=worker.getLevelSize(worker.cntLevel-1);
                    tmpSig2=worker.getSig2();
                }
                else ML = queryByte + 1;
            }
            avgM+=1.0*ML/TEST_CASE_M;
            avgH+=tmpH/TEST_CASE_M;
            avgTopCap+=tmpTopCap/TEST_CASE_M;
            avgSig2+=tmpSig2/TEST_CASE_M;
        }
        String result_s = "\t" + avgM+"\t" + avgH+"\t" + avgTopCap+"\t"+avgSig2;
        System.out.println("\t\ttargetAE:" + targetAE +result_s);
        err_result.set(RESULT_LINE, err_result.get(RESULT_LINE).concat("\t\t"+result_s));
    }


    public void testQuantiles(int queryN, int queryByte, double[] Qs){
        double[] Errs = new double[Qs.length];
        double[] query_a = new double[queryN];
        int[] LL = new int[TEST_CASE];
        int[] RR = new int[TEST_CASE];
        Random random = new Random(233);
        for (int i = 0; i < TEST_CASE; i++) {
            LL[i] = random.nextInt(N - queryN + 1);
            RR[i] = LL[i] + queryN;
        }

        for (int T = 0; T < TEST_CASE; T++) {
            int L = LL[T], R = RR[T];

            KLLDupliPair worker = new KLLDupliPair(queryByte);
            for (int i = L; i < R; i++)
                worker.update(dataToLong(a[i]));

            if (R - L >= 0) System.arraycopy(a, L, query_a, 0, R - L);
            Arrays.sort(query_a);

            for (int qid=0;qid<Qs.length;qid++) {
                double q=Qs[qid];
                int query_rank = (int) (q * queryN);
                double full_v = longToResult(worker.findMinValueWithRank(query_rank));
                int full_delta_rank = getDeltaRank(query_a, queryN, full_v, query_rank);
                double full_relative_err = 1.0 * full_delta_rank / (queryN);
                Errs[qid] += Math.abs(full_relative_err) / (TEST_CASE);
                System.out.println("\t\tq:"+ fnum.format(q)+"\testi_v:"+(full_v)+"\t\trel_err:"+Math.abs(full_relative_err)+"\tfull_delta_rank:"+full_delta_rank);
            }
        }
        for(int i=0;i<Qs.length;i++)
            System.out.println("\tQ:\t"+Qs[i]+"\t"+(i+1)+"\t"+Errs[i]);
    }

    public static void setTestCase(int tc) {
        TEST_CASE = tc;
    }


    public static void main(String[] args) throws IOException {
        long START = new Date().getTime();
        IntArrayList rep;
        System.out.println("KLLDupliPair\nTEST_CASE=" + TEST_CASE);
//        KLLDupliPair worker = new KLLDupliPair(8*16);
//        for(int i=0;i<32;i++){
//            worker.update(dataToLong((double)(i/3)));
//            worker.showNum();
//        }
//        for(int i=32;i<40;i++){
//            worker.update(dataToLong((double)(i)));
//            worker.showNum();
//        }
//        System.out.println("\t#20 Val:\t"+worker.longToResult(worker.findMinValueWithRank(20)));
//        sketchFT.showNum();
//        sketch.showNum();
//        for(int i=0;i<8;i++)System.out.println("\t\t"+sketchFT.findMinValueWithRank(i*15)+"\t\t"+sketch.findMinValueWithRank(i*15));
//        System.out.println("\t\t"+sketchFT.getApproxRank(6)+"\t"+sketch.getApproxRank(6));

        MainForKLLDupliPair main;
//
//        for (int dataType = 0; dataType <= 0; dataType++) { // CHECK IT
//            main = new MainForKLLDupliPair();
//            for (int queryN : new int[]{10000000})
//                for (int queryByte : new int[]{1024*128})
//                    for (double alpha : new double[]{0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0,1.1,1.2,1.3,1.4/*,1.5/*,1.6,1.7,1.8/*,1.9,2.0,2.1,2.2,2.3,2.4*/}) {
//                        if (dataType == 0) {
//                            err_result.add("N:\t" + queryN + ", " + "\tqueryByte:\t" + queryByte + ", " + "\talpha:\t" + alpha + "\t");
//                            time_result.add("N:\t" + queryN + ", " + "\tqueryByte:\t" + queryByte + ", " + "\talpha:\t" + alpha + "\t");
//                        }
//                        main.prepareZipf(dataType,alpha);
//                        main.testError(queryN, queryByte);
//                        main.RESULT_LINE++;
//                    }
//        }System.out.println("\nKLLDupliPair Zipf\tTEST_CASE=" + TEST_CASE);
//        System.out.println("Error rate:");
//        for (String s : err_result)
//            System.out.println(s);


        String muS="5";
        for (int dataType = 5; dataType <= 5; dataType++) { // CHECK IT
            main = new MainForKLLDupliPair();
            for (int queryN : new int[]{10000000})
                for (int queryByte : new int[]{1024*128})
                    for (double sigma : new double[]{0.2,0.4,0.6,0.8,1,1.2,1.4,1.6,1.8,2}) {
                        err_result.add("N:\t" + queryN + ", " + "\tqueryByte:\t" + queryByte + ", " + "\tsigma:\t" + sigma + "\t");
                        time_result.add("N:\t" + queryN + ", " + "\tqueryByte:\t" + queryByte + ", " + "\tsigma:\t" + sigma + "\t");
                        main.prepareLognormal(muS,sigma);
                        main.testError(queryN, queryByte);
                        main.RESULT_LINE++;
                    }
        }
        System.out.println("\nKLLDupliPair Lognormal mu="+muS+"\tTEST_CASE=" + TEST_CASE);
        System.out.println("Error rate:");
        for (String s : err_result)
            System.out.println(s);



//        for (int dataType = 0; dataType <= 0; dataType++) { // CHECK IT
//            main = new MainForKLLDupliPair();
//            for (int queryN : new int[]{10000000})
//                for (int queryByte : new int[]{1024*64})
//                    for (double alpha : new double[]{/*0.5,0.75,1,*/1.25,1.5,1.75,2,2.25,2.5,2.75,3,3.25,3.5}) {
//                        err_result.add("N:\t" + queryN + ", " + "\tqueryByte:\t" + queryByte + ", " + "\talpha:\t" + alpha + "\t");
//                        time_result.add("N:\t" + queryN + ", " + "\tqueryByte:\t" + queryByte + ", " + "\talpha:\t" + alpha + "\t");
//                        main.preparePareto(dataType,alpha);
//                        main.testError(queryN, queryByte);
//                        main.RESULT_LINE++;
//                    }
//        }
//        System.out.println("\nKLLDupliPair Pareto\tTEST_CASE=" + TEST_CASE);
//        System.out.println("Error rate:");
//        for (String s : err_result) System.out.println(s);
//
//
//        err_result=new ArrayList<>();
//        time_result=new ArrayList<>();
//        rep=new IntArrayList();//for(int i:new int[]{96})rep.add(i);
////        for(int i=4;i<=64;i+=4)rep.add(i);
////        for(int i=64+32;i<=1024;i+=32)rep.add(i);
//        rep.add(256);
//        for (int dataType = 0; dataType <= 0; dataType++) { // CHECK IT
//            main = new MainForKLLDupliPair();
//            for (int queryN : new int[]{10000000})
//                for (int queryByte : new int[]{1024*128})
////                    for (int repeat : new int[]{4,6,8,12,16,24,32,48,64,96,128,32*5,32*6,32*7,32*8,32*9,32*10,256,384,512,768,1024}) {
//                        for (int repeat : rep) {
//                        if (dataType == 0) {
//                            err_result.add("N:" + queryN + ", " + "queryByte:" + queryByte + ", " + "repeat:\t" + repeat + "\t");
//                            time_result.add("N:" + queryN + ", " + "queryByte:" + queryByte + ", " + "repeat:\t" + repeat + "\t");
//                        }
//                        main.prepareUniform(dataType,repeat,false);
//                        main.testError(queryN, queryByte);
//                        main.RESULT_LINE++;
//                    }
//        }
//        System.out.println("\nKLLDupliPair Uniform\tTEST_CASE=" + TEST_CASE);
//        System.out.println("Error rate:");
//        for (String s : err_result)
//            System.out.println(s);



//        err_result=new ArrayList<>();
//        time_result=new ArrayList<>();
//        rep=new IntArrayList();//for(int i:new int[]{20,24})rep.add(i);
////        for(int i=4;i<=64;i+=4)rep.add(i);
////        for(int i=64+16;i<=1024;i+=16)rep.add(i);
//        rep.add(1024);
//        for (int dataType = 0; dataType <= 0; dataType++) { // CHECK IT
//            main = new MainForKLLDupliPair();
//            for (int queryN : new int[]{10000000})
//                for (double targetAE: new double[]{1e-5})
//                    for (int repeat : rep) {
//                        if (dataType == 0) {
//                            err_result.add("N:" + queryN + ", " + "targetAE:" + targetAE + ", " + "repeat:\t" + repeat + "\t");
//                        }
//                        main.prepareUniform(dataType,repeat,false);
//                        main.testM(queryN, targetAE);
//                        main.RESULT_LINE++;
//                    }
//        }
//        System.out.println("\nTEST_CASE=" + TEST_CASE_M);
//        System.out.println("Space cost:");
//        for (String s : err_result)
//            System.out.println(s);


//        err_result=new ArrayList<>();
//        time_result=new ArrayList<>();
//        IntArrayList Ms=new IntArrayList();//for(int i=256;i<=1024;i+=64)Ms.add(i);
//        Ms.add(704);//Ms.add(832);Ms.add(896);
//        for (int dataType = startType; dataType <= endType; dataType++) { // CHECK IT
//            System.out.println("DATASET:"+dataType);
//            main = new MainForKLLDupliPair();
//            for (int queryN : new int[]{50000000})
//                for (int queryKB : Ms){
//                    if (dataType == startType) {
//                        err_result.add("N:\t" + queryN + "\tqueryKB:\t" + queryKB+ "\t");
//                        time_result.add("N:\t" + queryN + "\tqueryKB:\t" + queryKB+ "\t");
//                    }
//                    main.prepareA(dataType);
//                    main.testError(queryN, queryKB*1024);
//                    main.RESULT_LINE++;
//                }
//            System.out.println();
//        }
//        System.out.println("\nError rate:");for (String s : err_result) System.out.println(s);
//

//        err_result=new ArrayList<>();
//        time_result=new ArrayList<>();
//        IntArrayList Ns=new IntArrayList();//for(double i=1e6;i<1.01*2e7;i+=Math.min(i,2e6))Ns.add((int)i);//for(double i=1e7;i<1.1e8;i+=1e7)Ns.add((int)i);//
//        Ns.add((int)2e7);
//        for (int dataType = startType; dataType <= endType; dataType++) { // CHECK IT
//            System.out.println("DATASET:"+dataType);
//            main = new MainForKLLDupliPair();
//            for (int queryN : Ns)
//                for (int queryByte : new int[]{/*1024*512*/1024*64}){
//                    if (dataType == startType) {
//                        err_result.add("N:\t" + queryN + "\tqueryByte:\t" + queryByte+ "\t");
//                        time_result.add("N:\t" + queryN + "\tqueryByte:\t" + queryByte+ "\t");
//                    }
//                    main.prepareA(dataType);
//                    main.testError(queryN, queryByte);
//                    main.RESULT_LINE++;
//                }
//            System.out.println();
//        }
//        System.out.println("\nError rate:");
//        for (String s : err_result)System.out.println(s);

//        DoubleArrayList Qs=new DoubleArrayList();
////        Qs.add(0.75);
//        for(int i=1;i<=5;i++)Qs.add(i/100.0);for(int i=10;i<=90;i+=5)Qs.add(i/100.0);for(int i=95;i<=99;i++)Qs.add(i/100.0);
//        for (int dataType = startType; dataType <= endType; dataType++){
//            System.out.println("\n----dataType:"+dataType+"--------");
//            main = new MainForKLLDupliPair();
//            main.prepareA(dataType);
////            main.testQuantiles((int)5e7,512*1024,Qs.toDoubleArray());
//            main.testQuantiles((int)1e7,128*1024,Qs.toDoubleArray());
//        }

        System.out.println("\t\tALL_TIME:" + (new Date().getTime() - START));
    }
}
