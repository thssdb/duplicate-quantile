import it.unimi.dsi.fastutil.ints.IntOpenHashSet;
import it.unimi.dsi.util.XoRoShiRo128PlusPlusRandomGenerator;
import it.unimi.dsi.util.XorShift1024StarPhiRandomGenerator;
import org.apache.commons.math3.distribution.ParetoDistribution;
import org.apache.commons.math3.distribution.ZipfDistribution;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.util.Arrays;
import java.util.Date;

/** Compare sketchByPage & sketchByChunkDivide.
 * Different ChunkSize.
*/
public class MainForDDParameter {
    int dataType=-233;
    static int startType=3,endType=3;
    static int[] dataSetSize = new int[]{(int)3.1e7,(int)(3.1e7),(int)3.1e7,(int)3.1e7,(int)3.1e7,(int)3.1e7};
    static double[] a;


    public void prepareA(int dataType) throws IOException {
        if (a == null) a = new double[dataSetSize[dataType]];
        this.dataType = dataType;
        if (dataType == 4) {
//            LogNormalDistribution log21=new LogNormalDistribution(new XoRoShiRo128PlusPlusRandomGenerator(233),1,2);
//            for (int i = 0; i < Ns[dataType]; i++) a[i] = log21.sample();
            XoRoShiRo128PlusPlusRandomGenerator random=new XoRoShiRo128PlusPlusRandomGenerator(233);
            for (int i = 0; i < dataSetSize[dataType]; i++) a[i] =
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
            if (cntN == dataSetSize[dataType]) break;
        }
    }
    public void prepareZipf(int dataType,double alpha) throws IOException {
        if (a == null) a = new double[dataSetSize[dataType]];
        this.dataType = dataType;

        ZipfDistribution dis = new ZipfDistribution(new XorShift1024StarPhiRandomGenerator(233), dataSetSize[0]/233, alpha);
        for (int i = 0; i < dataSetSize[0]; i++)
            a[i] = dis.sample();
    }
    public void preparePareto(int dataType,double alpha) throws IOException {
        int N=dataSetSize[dataType];
        if (a == null) a = new double[N];
        this.dataType = dataType;
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

    public void prepareLognormal(int dataType,String muS,double sigma) throws IOException {
        int N=dataSetSize[dataType];
        if (a == null||a.length<N) a = new double[N];
        BufferedReader reader = new BufferedReader(new FileReader("Lognormal3E7Mu"+muS+"Sigma"+(int)(sigma*10+0.1)+".txt"));
        String line;
        int cntN = 0;
        while ((line = reader.readLine()) != null) {
            a[cntN++] = Double.parseDouble(line);
            if (cntN == N) break;
        }
    }
    public double[] getMinVAlpha(int n,int maxSeriByte){
        double minV=Double.MAX_VALUE;
        for(int i=0;i<n;i++)minV=Math.min(minV,a[i]);
//        minV=1.0;//TODO danger minV
        int maxBucketNum=maxSeriByte/DDSketchPositiveForDupli.bucketNumPerByteInMemoryForApprox;
        double alpha=Math.pow(0.95,100);// TODO danger last used alpha: 144
        double resultAlpha=alpha;
        while(true && alpha>1e-7) {
            System.out.println("\t\ttrying alpha:\t"+alpha);
            double gamma = 2 * alpha / (1 - alpha) + 1;
            double multiplier = Math.log(Math.E) / (Math.log1p(gamma - 1));
            boolean valid=true;
            int maxKey=0;
            IntOpenHashSet set=new IntOpenHashSet();
            for (int i = 0; i < n; i++) {
                double v = a[i] - minV + 1;
                int key = (int) Math.ceil(Math.log(v) * multiplier);
                maxKey=Math.max(maxKey,key);
                set.add(key);
            }
            if(set.size()>maxBucketNum) {
                valid = false;
                System.out.println("\tINVALID\tset.size():"+set.size());
            }
            if(!valid)break;
            else {
                resultAlpha=alpha;
                alpha*=0.8;
            }
        }
        return new double[]{minV,resultAlpha};
    }



    public static void main(String[] args) throws IOException {
        long START_T = new Date().getTime();

//        DDSketchPositiveForExact a=new DDSketchPositiveForExact(0.01,50);
//        DoubleArrayList va=new DoubleArrayList();
//        Random random=new Random(233);
//        for(int i=0;i<50;i++){
//            double v=1+random.nextDouble();
//            a.insert(v);
//            va.add(v);
//        }
//        System.out.println("\t\t|a|:"+a.sketch_size());
//
//        DDSketchPositiveForExact b=new DDSketchPositiveForExact(0.01,20);
//        DoubleArrayList vb=new DoubleArrayList();
//        for(int i=0;i<50;i++){
//            double v=1+random.nextDouble();
//            b.insert(v);
//            vb.add(v);
//        }
//        System.out.println("\t\t|b|:"+b.sketch_size());
//
//        b.mergeWithDeserialized(a);
//        System.out.println("\t\t|b+a|:"+b.sketch_size());
//        System.out.println("\t\t|b+a|.N:"+b.total_count());

        MainForDDParameter main=new MainForDDParameter();

//        IntArrayList Ms=new IntArrayList();//Ms.add(64);
//        for(int i = 64; i<=256; i+=16)Ms.add(i);
//        for(int dataType=startType;dataType<=endType;dataType++)
//            for(int seriKB:Ms){
//            main.prepareA(dataType);
//            double[] tmp=main.getMinVAlpha(dataSetSize[dataType], seriKB*1024);
//            System.out.println("\t\tdataset:\t"+dataType+"\tseriKB:\t"+seriKB+"\talpha:\t"+tmp[1]+"\t"+Math.log(tmp[1])/Math.log(1/0.95)+"\t\tdataSetMinV:\t"+tmp[0]);
//        }


//        for(int seriKB:new int[]{128})
//            for (double alpha : new double[]{0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0,1.1,1.2,1.3,1.4/*,1.6,1.8,2.0/*,1.5/*,1.6,1.7,1.8/*,1.9,2.0,2.1,2.2,2.3,2.4*/}) {
//                main.prepareZipf(0,alpha);
//                double[] tmp=main.getMinVAlpha(dataSetSize[0], seriKB*1024);
//                System.out.println("\t\tdataset:\t"+"Zipf"+"\tseriKB:\t"+seriKB+"\talpha:\t"+tmp[1]+"\t"+Math.log(tmp[1])/Math.log(1/0.95)+"\t\tdataSetMinV:\t"+tmp[0]);
//            }

//        for(int seriKB:new int[]{64})
//            for (double alpha : new double[]{0.5,0.75,/*1,1.25,1.5,1.75,2,2.25,2.5,2.75,3,*/3.25,3.5}) {
//                main.preparePareto(0,alpha);
//                double[] tmp=main.getMinVAlpha(dataSetSize[0], seriKB*1024);
//                System.out.println("\t\tdataset:\t"+"Pareto"+"\tseriKB:\t"+seriKB+"\talpha:\t"+tmp[1]+"\t"+Math.log(tmp[1])/Math.log(1/0.95)+"\t\tdataSetMinV:\t"+tmp[0]);
//            }

        String muS="5";
        for (int seriKB : new int[]{128})
            for (double sigma : new double[]{0.2,0.4,0.6,0.8,1,1.2,1.4,1.6,1.8,2}){
                main.prepareLognormal(5,muS,sigma);
                double[] tmp=main.getMinVAlpha(dataSetSize[5], seriKB*1024);
                System.out.println("\t\tdataset:\t"+"Lognormal\t"+muS+muS+"\tseriKB:\t"+seriKB+"\talpha:\t"+tmp[1]+"\t"+Math.log(tmp[1])/Math.log(1/0.95)+"\t\tdataSetMinV:\t"+tmp[0]);
            }

                System.out.println("\t\t\tALL_TIME:"+(new Date().getTime()-START_T));
    }
}
