import it.unimi.dsi.fastutil.longs.LongArrayList;
import it.unimi.dsi.util.XorShift1024StarPhiRandomGenerator;
import org.apache.commons.math3.distribution.ZipfDistribution;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.text.DecimalFormat;
import java.util.ArrayList;
import java.util.Date;
import java.util.Random;

public class MainForDomainChange {
    int dataType;
    static int zipfDataType=0,uniformDataType=0,synDataType = zipfDataType+uniformDataType;
    static int startType = 1, endType = 2;
    static int N = (int)3e7;
    public static int TEST_CASE = 15; // CHECK IT
    static double[] a;
    static ArrayList<String> err_result = new ArrayList<>();
    static ArrayList<String> time_result = new ArrayList<>();
    int RESULT_LINE = 0;

    Random random = new Random(233);

    public void prepareA(double[] a, int dataType, int AN) throws IOException {
//        if (a == null||a.length<AN) a = new double[AN];
//        System.out.println("\t\t\tan:"+a.length);
        this.dataType = dataType;
        if(dataType>=synDataType) {
//            System.out.println("\t\t??\t"+dataType);
            BufferedReader reader = null;
            if (dataType == synDataType+0) reader = new BufferedReader(new FileReader(new File("DupliTorqueVoltage.txt")));
            if (dataType == synDataType+1) reader = new BufferedReader(new FileReader(new File("DupliECommercePrice.txt")));
            if (dataType == synDataType+2) reader = new BufferedReader(new FileReader(new File("DupliElectric.txt")));
            reader.readLine(); // ignore first line.
            String line;
            int cntN = 0;
            while ((line = reader.readLine()) != null) {
                a[cntN++] = Double.parseDouble(line);
                if (cntN == AN) break;
            }
        }else{
            double[] zipfP = new double[]{0.8,0.8};
            double[] uniformP = new double[]{1000,10000};
            if(dataType<zipfDataType){
                ZipfDistribution dis = new ZipfDistribution(new XorShift1024StarPhiRandomGenerator(233), N/100, zipfP[dataType]);
                for (int i = 0; i < AN; i++) a[i] = 1.0*dataType*N+dis.sample();
            }else{
                for(int i=0;i<AN;i++)a[i]=1.0*dataType*N+i/uniformP[dataType-zipfDataType];
                for(int i=1;i<AN;i++){
                    int p=random.nextInt(i);double aa=a[p];a[p]=a[i];a[i]=aa;
                }
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
    public int getDeltaRank(SortedList<Double> sortedA, double v, int targetRank) {
        int rank_L = sortedA.getNodeNumLessThanValue(v);
        double v2=longToResult(dataToLong(v)+1);
        int rank_R = sortedA.getNodeNumLessThanValue(v2) - 1;
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

    DecimalFormat fnum = new DecimalFormat("#.##E0");

    public int getPos(LongArrayList list,long v){
        int L=0,R=list.size()-1;
        while(L<R){
            int mid=(L+R)/2;
            if(list.getLong(mid)>=v)R=mid;
            else L=mid+1;
        }return L;
    }


//    public void testError(int queryN, int queryByte) throws IOException {
//        int sketchM = Integer.highestOneBit(queryByte / (8 + 8 + 2));
//        long full_time = 0, merge_time = 0;
//        double err_full = 0, err_merge = 0;
//        double MMP_full = 0;
//        double[] query_a = new double[queryN];
//
//        int[] LL = new int[TEST_CASE];
//        int[] RR = new int[TEST_CASE];
//        Random random = new Random(233);
//        for (int i = 0; i < TEST_CASE; i++) {
//            LL[i] = random.nextInt(N - queryN + 1);
//            RR[i] = LL[i] + queryN;
//        }
//        long errBound99=0;
//        long PairT0T0=0,PairFF=0,PairTF=0,PairT1T2=0;
//        for (int T = 0; T < TEST_CASE; T++) {
//            int L = LL[T], R = RR[T];
//
//            full_time -= new Date().getTime();
//            KLLDupliFT worker = new KLLDupliFT(queryByte);
//            for (int i = L; i < R; i++)
//                worker.update(dataToLong(a[i]));
//            PairT0T0+=worker.PairT0T0;
//            PairFF+=worker.PairFF;
//            PairTF+=worker.PairTF;
//            PairT1T2+= worker.PairT1T2;
//            errBound99=worker.queryRankErrBound(0.99);
////            if(T==0) {
////                worker.show();
//////                worker.showNum();
////            }
////            System.out.println("\t\t|worker|:\t"+worker.getMaximumMapCapacity()+"\t\tsketchM:"+sketchM);
//            full_time += new Date().getTime();
//
//            if (R - L >= 0) System.arraycopy(a, L, query_a, 0, R - L);
//            Arrays.sort(query_a);
//
//            double q_add = 0.0001, q_start = q_add, q_end = 1-q_add, q_count = Math.floor((q_end - q_start - 1e-10) / q_add) + 1;
//            for (double q = q_start; q < q_end + 1e-10; q += q_add) {
//                int query_rank = (int) (q * queryN);
//                double full_v = longToResult(worker.findMinValueWithRank(query_rank));
//                int full_delta_rank = getDeltaRank(query_a, queryN, full_v, query_rank);
//                double full_relative_err = 1.0 * full_delta_rank / (queryN);
//                err_full += Math.abs(full_relative_err) / (q_count * TEST_CASE);
//            }
//        }
////        err_merge/=(queryPageNum*pageN);err_full/=(queryPageNum*pageN);
////        System.out.println("\t\t" + queryN +/*"\t"+err_merge+*/"\t" + err_full+"\t99%bound:\t"+1.0*errBound99/queryN+"\tPairs(TT,FF,TF,T1T2):\t"+PairT0T0+"\t"+PairFF+"\t"+PairTF+"\t"+PairT1T2);
//
//        err_result.set(RESULT_LINE, err_result.get(RESULT_LINE).concat("\t\t\t" + err_full+"\t"+1.0*errBound99/queryN));
//        time_result.set(RESULT_LINE, time_result.get(RESULT_LINE).concat("\t\t\t" + full_time));
//
//    }



    public void testDomainAB(int dataType1,int dataType2, int queryByte) throws IOException {
        double[] a1=new double[N/2],a2=new double[N/2];
        prepareA(a1,dataType1,N/2);
        prepareA(a2,dataType2,N/2);
        a = new double[N];
        for(int i=0;i<N/2;i++){a[i]=a1[i];a[i+N/2]=a2[i];}
//        prepareA(a,dataType2,N); // CHECK IT!

        int STEP_NUM = 100;
        double[] errs_ori=new double[STEP_NUM],errs_FT=new double[STEP_NUM];
        for(int T=0;T<TEST_CASE;T++) {
            System.out.println("TEST_CASE:\t"+T);
            KLLSketchLazyExactPriori workerOri = new KLLSketchLazyExactPriori(queryByte);
            KLLDupliFT workerFT = new KLLDupliFT(queryByte);
            SortedList<Double> sortedA = new SortedList<>(Double::compareTo);
            for (int step = 1; step <= STEP_NUM; step++) {
//            if(step>=3)break;
                for (int i = (step - 1) * (N / STEP_NUM); i < step * (N / STEP_NUM); i++) {
                    workerOri.update(dataToLong(a[i]));
                    workerFT.update(dataToLong(a[i]));
                    sortedA.add(a[i]);
                }
//            workerFT.checkN();
                double q_add = 0.0001, q_start = q_add, q_end = 1 - q_add, q_count = Math.floor((q_end - q_start - 1e-10) / q_add) + 1;
                int queryN = step * (N / STEP_NUM);
                double step_err_ori = 0, step_err_FT = 0;
                for (double q = q_start; q < q_end + 1e-10; q += q_add) {
                    int query_rank = (int) (q * queryN);
                    double ori_v = longToResult(workerOri.findMinValueWithRank(query_rank));
                    int ori_delta_rank = getDeltaRank(sortedA, ori_v, query_rank);
                    double ori_relative_err = 1.0 * ori_delta_rank / N;//(queryN);
                    step_err_ori += Math.abs(ori_relative_err) / (q_count);

                    double FT_v = longToResult(workerFT.findMinValueWithRank(query_rank));
                    int FT_delta_rank = getDeltaRank(sortedA, FT_v, query_rank);
                    double FT_relative_err = 1.0 * FT_delta_rank / N;//(queryN);
                    step_err_FT += Math.abs(FT_relative_err) / (q_count);
                }
                errs_ori[step-1]+=step_err_ori/TEST_CASE;
                errs_FT[step-1]+=step_err_FT/TEST_CASE;

                if(T==0)
                System.out.println("----[STEP]:\t" + step
                    + "\t\tBOUND99:\tOriKLL:\t" + workerOri.queryRankErrBound(0.99) + "\tFTKLL:\t" + workerFT.queryRankErrBound(0.99)
                    + "\t\tERR:\tOriKLL:\t" + fnum.format(step_err_ori) + "\tFTKLL:\t" + fnum.format(step_err_FT));
            }
        }
        System.out.println("\t\tstepID\t\toriKLL\tDupliKLL");
        for(int i=0;i<STEP_NUM;i++)System.out.println("\t----[STEP]:\t"+i+"\t\t"+fnum.format(errs_ori[i])+"\t"+fnum.format(errs_FT[i]));
//        long errBound99Ori,errBound99FT=0;
//        long PairT0T0=0,PairFF=0,PairTF=0,PairT1T2=0;
//        for (int T = 0; T < TEST_CASE; T++) {
//            int L = LL[T], R = RR[T];
//
//            full_time -= new Date().getTime();
//            KLLDupliFT worker = new KLLDupliFT(queryByte);
//            for (int i = L; i < R; i++)
//                worker.update(dataToLong(a[i]));
//            PairT0T0+=worker.PairT0T0;
//            PairFF+=worker.PairFF;
//            PairTF+=worker.PairTF;
//            PairT1T2+= worker.PairT1T2;
//            errBound99=worker.queryRankErrBound(0.99);
////            if(T==0) {
////                worker.show();
//////                worker.showNum();
////            }
////            System.out.println("\t\t|worker|:\t"+worker.getMaximumMapCapacity()+"\t\tsketchM:"+sketchM);
//            full_time += new Date().getTime();
//
//            if (R - L >= 0) System.arraycopy(a, L, query_a, 0, R - L);
//            Arrays.sort(query_a);
//
//            double q_add = 0.0001, q_start = q_add, q_end = 1-q_add, q_count = Math.floor((q_end - q_start - 1e-10) / q_add) + 1;
//            for (double q = q_start; q < q_end + 1e-10; q += q_add) {
//                int query_rank = (int) (q * queryN);
//                double full_v = longToResult(worker.findMinValueWithRank(query_rank));
//                int full_delta_rank = getDeltaRank(query_a, queryN, full_v, query_rank);
//                double full_relative_err = 1.0 * full_delta_rank / (queryN);
//                err_full += Math.abs(full_relative_err) / (q_count * TEST_CASE);
//            }
//        }
////        err_merge/=(queryPageNum*pageN);err_full/=(queryPageNum*pageN);
////        System.out.println("\t\t" + queryN +/*"\t"+err_merge+*/"\t" + err_full+"\t99%bound:\t"+1.0*errBound99/queryN+"\tPairs(TT,FF,TF,T1T2):\t"+PairT0T0+"\t"+PairFF+"\t"+PairTF+"\t"+PairT1T2);
//
//        err_result.set(RESULT_LINE, err_result.get(RESULT_LINE).concat("\t\t\t" + err_full+"\t"+1.0*errBound99/queryN));
//        time_result.set(RESULT_LINE, time_result.get(RESULT_LINE).concat("\t\t\t" + full_time));

    }


    public static void main(String[] args) throws IOException {
        long START = new Date().getTime();
        System.out.println("KLLDupliTF\t1 Domains\nTEST_CASE=" + TEST_CASE);

        MainForDomainChange main;
        int[] domainA=new int[]{0};
        int[] domainB=new int[]{1};
        for(int pairID=0;pairID<domainA.length;pairID++){
            main=new MainForDomainChange();
            main.testDomainAB(domainA[pairID],domainB[pairID],1024*256);
        }

//        for (int dataType = 0; dataType <= 0; dataType++) { // CHECK IT
//            main = new MainForDomainChange();
//            for (int queryN : new int[]{4000000})
//                for (int queryByte : new int[]{1024*256})
//                    for (double alpha : new double[]{0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0,1.1,1.2,1.3,1.4,1.5/*,1.6,1.7,1.8,1.9,2.0,2.1,2.2,2.3,2.4*/}) {
//                        if (dataType == 0) {
//                            err_result.add("N:" + queryN + ", " + "queryByte:" + queryByte + ", " + "alpha:" + alpha + "\t");
//                            time_result.add("N:" + queryN + ", " + "queryByte:" + queryByte + ", " + "alpha:" + alpha + "\t");
//                        }
//                        main.prepareZipf(dataType,alpha);
//                        main.testError(queryN, queryByte);
//                        main.RESULT_LINE++;
//                    }
//        }
//        System.out.println("\nError rate:");
//        for (String s : err_result)
//            System.out.println(s);


        System.out.println("\t\tALL_TIME:" + (new Date().getTime() - START));
    }
}
