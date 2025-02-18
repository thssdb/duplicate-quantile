import it.unimi.dsi.fastutil.doubles.DoubleOpenHashSet;
import it.unimi.dsi.fastutil.ints.Int2DoubleOpenHashMap;
import it.unimi.dsi.fastutil.ints.IntIntPair;
import it.unimi.dsi.fastutil.objects.Object2DoubleOpenHashMap;
import it.unimi.dsi.util.XoRoShiRo128PlusPlusRandomGenerator;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Date;

public class MainForNonDupliDataTxt {
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


    static void prepareData(String filename) throws IOException {
        if (a == null) a = new double[N];
        BufferedReader reader = null;
        reader = new BufferedReader(new FileReader(new File(filename)));//4_SpacecraftThruster
        reader.readLine(); // ignore first line.
        String line;
        int cntN = 0;
        DoubleOpenHashSet hash=new DoubleOpenHashSet(N);
        while ((line = reader.readLine()) != null) {
            double tmp= Double.parseDouble(line);
            if(!hash.contains(tmp)) {
                a[cntN++] = tmp;
                hash.add(tmp);
            }
            else continue;
            if (cntN == N) break;
        }
    }
    static void prepareDHard() throws IOException {
        if (a == null) a = new double[N];
        XoRoShiRo128PlusPlusRandomGenerator random=new XoRoShiRo128PlusPlusRandomGenerator(233);
        for (int i = 0; i < N; i++) a[i] =
            Math.pow(-1, random.nextInt(2)) * Math.pow(10.0, (2 * Math.pow(random.nextDouble(), 8) - 1) * 300); // D_hard
    }


    public static void main(String[] args) throws IOException {
        long START = new Date().getTime();
        {
//            String file="Binance.txt";
//            prepareData(file);
//            String filename = "NonDupli"+file;
            prepareDHard();
            String filename = "DHard.txt";
            System.out.println("\t\t"+filename);
            FileWriter fw = new FileWriter(filename);
            for(int i=0;i<N;i++)fw.write(a[i]+"\n");
            fw.close();
        }
        System.out.println("\t\tALL_TIME:" + (new Date().getTime() - START));
    }
}
