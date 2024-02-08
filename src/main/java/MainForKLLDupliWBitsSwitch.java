import it.unimi.dsi.fastutil.doubles.Double2DoubleOpenHashMap;
import it.unimi.dsi.fastutil.doubles.Double2LongOpenHashMap;
import it.unimi.dsi.fastutil.doubles.DoubleArrayList;
import it.unimi.dsi.fastutil.longs.LongArrayList;
import it.unimi.dsi.util.XorShift1024StarPhiRandom;
import it.unimi.dsi.util.XorShift1024StarPhiRandomGenerator;
import org.apache.commons.math3.distribution.ZipfDistribution;

import java.io.IOException;
import java.text.DecimalFormat;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Date;
import java.util.Random;

public class MainForKLLDupliWBitsSwitch {
    int dataType;
    static int startType = 0, endType = 0;
    static int pageN = 8192;
    static int N = /*55000000/pageN*pageN*/10000000*3, pageNum = N / pageN; // CHECK IT
    public static int TEST_CASE = 1; // CHECK IT
    static double[] a;
    static ArrayList<String> err_result = new ArrayList<>();
    static ArrayList<String> time_result = new ArrayList<>();
    int RESULT_LINE = 0;

    Random random = new Random(233);

    public void prepareZipf(int dataType,double alpha) throws IOException {
        if (a == null) a = new double[N];
        this.dataType = dataType;

        ZipfDistribution dis = new ZipfDistribution(new XorShift1024StarPhiRandomGenerator(233), N/233, alpha);
        for (int i = 0; i < N; i++)
            a[i] = dis.sample();
    }

    public void prepareUniform(int dataType,int repeat) throws IOException {
        if (a == null) a = new double[N];
        this.dataType = dataType;

//        UniformIntegerDistribution dis = new UniformIntegerDistribution(new XorShift1024StarPhiRandomGenerator(233),0,N/repeat);
//        for (int i = 0; i < N; i++)
//            a[i] = dis.sample();
        for(int i=0;i<N;i++)a[i]=(double)(i/repeat);
        for(int i=1;i<N;i++){
            int p=random.nextInt(i);double aa=a[p];a[p]=a[i];a[i]=aa;
        }
    }

    public void prepareAOrder(int dataType){
        Double2LongOpenHashMap map0 = new Double2LongOpenHashMap(N);
        for(int i=0;i<N;i++){
            map0.putIfAbsent(a[i],0);
            map0.addTo(a[i],1);
        }
        DoubleArrayList aSet=new DoubleArrayList(map0.keySet()),aSet2=new DoubleArrayList();
        XorShift1024StarPhiRandom rand = new XorShift1024StarPhiRandom(233);
        for(int i=0;i<aSet.size();i++){
            aSet2.add(aSet.getDouble(i));
            if(i>=1) {
                int tmpP = rand.nextInt(i + 1);
                double tmpV=aSet2.getDouble(tmpP);
                aSet2.set(tmpP,aSet2.getDouble(i));
                aSet2.set(i,tmpV);
            }
        }
        Double2DoubleOpenHashMap v2v=new Double2DoubleOpenHashMap();
        for(int i=0;i<aSet.size();i++)
            v2v.put(aSet.getDouble(i),aSet2.getDouble(i));
        for(int i=0;i<N;i++)a[i]=v2v.get(a[i]);
        if(dataType==0) return;

        Double2LongOpenHashMap map = new Double2LongOpenHashMap();
        for(int i=0;i<N;i++){
            map.putIfAbsent(a[i],0);
            map.addTo(a[i],1);
        }
        aSet=new DoubleArrayList(map.keySet());
        if(dataType==1){
            aSet.sort((x,y)->(Long.compare(-map.get(x),-map.get(y))));
            int tmpN=0;
            for(double val:aSet){
                for(long i=map.get(val);i>0;i--)
                    a[tmpN++]=val;
            }
        }else
        if(dataType==2){
            aSet.sort((x,y)->(Long.compare(map.get(x),map.get(y))));
            int tmpN=0;
            for(double val:aSet){
                for(long i=map.get(val);i>0;i--)
                    a[tmpN++]=val;
            }
        }else
        if(dataType==3){
            int tmpN=0,tmpC=1;
            while(tmpN<N){
                for(double val:new DoubleArrayList(map.keySet())){
                    if(map.get(val)>=tmpC)a[tmpN++]=val;
                    else map.remove(val);
                }
                tmpC++;
            }
        }

//        for(int i=0;i<5;i++)System.out.println("\t"+a[i]);
//        for(int i=0;i<5;i++)System.out.println("\t\t"+a[N-5+i]);
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


    public void testError(int wBits,int queryN, int queryByte) throws IOException {
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
        for (int T = 0; T < TEST_CASE; T++) {
            int L = LL[T], R = RR[T];

            full_time -= new Date().getTime();
            KLLDupliWBitsSwitch worker = new KLLDupliWBitsSwitch(queryByte,wBits);
            for (int i = L; i < R; i++)
                worker.update(dataToLong(a[i]));
            errBound99=worker.queryRankErrBound(0.99);
//            if(T==0) {
//                worker.show();
//                worker.showCompact();
//            }
//            System.out.println("\t\t|worker|:\t"+worker.getMaximumMapCapacity()+"\t\tsketchM:"+sketchM);
            full_time += new Date().getTime();

            if (R - L >= 0) System.arraycopy(a, L, query_a, 0, R - L);
            Arrays.parallelSort(query_a);

            double q_add = 0.0001, q_start = q_add, q_end = 1-q_add, q_count = Math.floor((q_end - q_start - 1e-10) / q_add) + 1;
            for (double q = q_start; q < q_end + 1e-10; q += q_add) {
                int query_rank = (int) (q * queryN);

                double full_v = longToResult(worker.findMinValueWithRank(query_rank));

                int full_delta_rank = getDeltaRank(query_a, queryN, full_v, query_rank);
                double full_relative_err = 1.0 * full_delta_rank / (queryN);
                err_full += Math.abs(full_relative_err) / (q_count * TEST_CASE);

            }
        }
//        err_merge/=(queryPageNum*pageN);err_full/=(queryPageNum*pageN);
        System.out.println("\t\t" + queryN +/*"\t"+err_merge+*/"\t" + err_full+"\t"+1.0*errBound99/queryN);

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
        MainForKLLDupliWBitsSwitch main;
        System.out.println("KLLDupliWBitsSwitch\nTEST_CASE=" + TEST_CASE+"\tOnlyDealFullyDupli:"+ KLLDupliWBitsSwitch.onlyDealFullyDupli);

//        for (int dataType = startType; dataType <= endType; dataType++) { // CHECK IT
//            main = new MainForKLLDupliWBitsSwitch();
////            for (int queryN : new int[]{200})
////                for (int queryByte : new int[]{4*65})
////                    for (double alpha : new double[]{2.0}) {
//            for (int wBits : new int[]{/*4*/1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16/**/})
//            for (int queryN : new int[]{10000000})
//                for (int queryByte : new int[]{1024*128})
//                    for (double alpha : new double[]{0.7/*0.1,0.2,0.3,0.4,0.5,0.6/*,0.7,0.8,0.9,1.0,1.1,1.2,1.3,1.4,1.5,1.6,1.7,1.8,1.9,2.0,2.1,2.2,2.3,2.4*/}) {
//                        if (dataType == startType) {
//                            err_result.add("N:" + queryN + ", " + "queryByte:" + queryByte + ", " + "alpha:" + alpha + ",wBits:"+wBits+"\t");
//                            time_result.add("N:" + queryN + ", " + "queryByte:" + queryByte + ", " + "alpha:" + alpha + ",wBits:"+wBits+"\t");
//                        }
//                        main.prepareZipf(dataType,alpha);
//                        main.prepareAOrder(dataType);
//                        main.testError(wBits,queryN, queryByte);
//                        main.RESULT_LINE++;
//                    }
//        }
//        System.out.println("\nError rate:");
//        for (String s : err_result)
//            System.out.println(s);


        err_result=new ArrayList<>();
        time_result=new ArrayList<>();
        for (int dataType = startType; dataType <= endType; dataType++) { // CHECK IT
            main = new MainForKLLDupliWBitsSwitch();
            for (int wBits : new int[]{1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16})
            for (int queryN : new int[]{10000000})
                for (int queryByte : new int[]{1024*128})
                    for (int repeat : new int[]{1/*,10,100,1000,2000,5000,10000,20000,50000,100000,200000,500000,1000000*/}) {
                        if (dataType == startType) {
                            err_result.add("N:" + queryN + ", " + "queryByte:" + queryByte + ", " + "repeat:" + repeat + "\t");
                            time_result.add("N:" + queryN + ", " + "queryByte:" + queryByte + ", " + "repeat:" + repeat + "\t");
                        }
                        main.prepareUniform(dataType,repeat);
//                        main.prepareAOrder(dataType);
                        main.testError(wBits, queryN, queryByte);
                        main.RESULT_LINE++;
                    }
        }
        System.out.println("KLLDupliWBitsSwitch\nTEST_CASE=" + TEST_CASE);
        System.out.println("\nError rate:");
        for (String s : err_result)
            System.out.println(s);

        System.out.println("\t\tALL_TIME:" + (new Date().getTime() - START));
    }
}
