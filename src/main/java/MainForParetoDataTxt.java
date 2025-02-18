import it.unimi.dsi.fastutil.ints.Int2DoubleOpenHashMap;
import it.unimi.dsi.fastutil.ints.IntIntPair;
import it.unimi.dsi.fastutil.objects.Object2DoubleOpenHashMap;
import it.unimi.dsi.util.XorShift1024StarPhiRandomGenerator;
import org.apache.commons.math3.distribution.ParetoDistribution;

import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Date;

public class MainForParetoDataTxt {
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


    static void preparePareto(double alpha) throws IOException {
        if (a == null) a = new double[N];
        ParetoDistribution dis = new ParetoDistribution(new XorShift1024StarPhiRandomGenerator(233), 1, alpha);
        for (int i = 0; i < N; i++)
            a[i] = dis.sample();
        double Epsilon=1+5e-4;
        for (int i = 0; i < N; i++)a[i]=Math.pow(Epsilon,Math.ceil(Math.log(a[i])/Math.log(Epsilon)));
        double[] b= Arrays.copyOf(a,N);
        Arrays.sort(b);
        int count=0;
        for(int i=1;i<N;i++)if(b[i]!=b[i-1])count++;
        System.out.println("\t[Pareto]\talpha:\t"+alpha+"\tEpsilon:\t"+Epsilon+"\t\tcount:\t"+count+"\t\tN:\t"+N);
    }


    public static void main(String[] args) throws IOException {
        long START = new Date().getTime();
        for (double alpha : new double[]{0.5,0.75,1,1.25,1.5,1.75,2,2.25,2.5,2.75,3,3.25,3.5}) {
            preparePareto(alpha);
            String filename = "Pareto3E7Alpha"+(int)(alpha*100+0.1)+".txt";
            System.out.println("\t\t"+filename);
            FileWriter fw = new FileWriter("E:\\KLL-Dupli-Synthetic-Data\\"+filename);
            for(int i=0;i<N;i++)fw.write(a[i]+"\n");
            fw.close();
        }
        System.out.println("\t\tALL_TIME:" + (new Date().getTime() - START));
    }
}
