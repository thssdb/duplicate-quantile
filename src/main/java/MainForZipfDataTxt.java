import it.unimi.dsi.fastutil.ints.Int2DoubleOpenHashMap;
import it.unimi.dsi.fastutil.ints.IntIntPair;
import it.unimi.dsi.fastutil.objects.Object2DoubleOpenHashMap;
import it.unimi.dsi.util.XoRoShiRo128PlusRandom;
import it.unimi.dsi.util.XorShift1024StarPhiRandomGenerator;
import org.apache.commons.math3.distribution.ZipfDistribution;

import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Date;

public class MainForZipfDataTxt {
    int dataType;
    static int startType = 1, endType = 2;
    static int pageN = 8192;
    static int N = (int)3.1e7, pageNum = N / pageN; // CHECK IT
    public static int TEST_CASE = 16,TEST_CASE_M = 1; // CHECK IT
    static double[] a;
    static ArrayList<String> err_result = new ArrayList<>();
    static ArrayList<String> time_result = new ArrayList<>();
    int RESULT_LINE = 0;
    static Object2DoubleOpenHashMap<IntIntPair> typeMem2Alpha;
    static Int2DoubleOpenHashMap type2MinV;


    static void prepareZipf(double alpha) throws IOException {
        if (a == null) a = new double[N];
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
//            a[i] = v2v[dis.sample()];
            a[i] = dis.sample();
    }


    public static void main(String[] args) throws IOException {
        long START = new Date().getTime();
        for (double alpha : new double[]{0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.1, 1.2, 1.3, 1.4}) {
            prepareZipf(alpha);
            String filename = "Zipf3E7Alpha"+(int)(alpha*10+0.1)+".txt";
            System.out.println("\t\t"+filename);
            FileWriter fw = new FileWriter(filename);
            for(int i=0;i<N;i++)fw.write(a[i]+"\n");
            fw.close();
        }
        System.out.println("\t\tALL_TIME:" + (new Date().getTime() - START));
    }
}
