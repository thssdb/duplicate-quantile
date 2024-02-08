import it.unimi.dsi.fastutil.ints.IntArrayList;
import it.unimi.dsi.fastutil.ints.IntOpenHashSet;
import it.unimi.dsi.util.XorShift1024StarPhiRandomGenerator;
import org.apache.commons.math3.distribution.ZipfDistribution;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.util.Date;

/** Compare sketchByPage & sketchByChunkDivide.
 * Different ChunkSize.
*/
public class MainForDDParameter {
    int dataType=-233;
    static int startType=3,endType=3;
    static int[] dataSetSize = new int[]{(int)3.1e7,(int)(3.1e7),(int)3.1e7,(int)1.1e8};
    static double[] a;


    public void prepareA(int dataType) throws IOException {
        if (a == null) a = new double[dataSetSize[dataType]];
        this.dataType = dataType;
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

    public double[] getMinVAlpha(int n,int maxSeriByte){
        double minV=Double.MAX_VALUE;
        for(int i=0;i<n;i++)minV=Math.min(minV,a[i]);// minV=0;//TODO danger minV
        int maxBucketNum=maxSeriByte/DDSketchPositiveForDupli.bucketNumPerByteInMemoryForApprox;
        double alpha=Math.pow(0.95,144);
        double resultAlpha=alpha;
        while(true) {
//            System.out.println("\t\ttrying alpha:\t"+alpha);
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
                alpha*=0.95;
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

        IntArrayList Ms=new IntArrayList();Ms.add(128);//for(int i = 256; i<=1024; i+=64)Ms.add(i);
        for(int dataType=startType;dataType<=endType;dataType++)
            for(int seriKB:Ms){
            main.prepareA(dataType);
            double[] tmp=main.getMinVAlpha(dataSetSize[dataType], seriKB*1024);
            System.out.println("\t\tdataset:\t"+dataType+"\tseriKB:\t"+seriKB+"\talpha:\t"+tmp[1]+"\t"+Math.log(tmp[1])/Math.log(1/0.95)+"\t\tdataSetMinV:\t"+tmp[0]);
        }


//        for(int seriKB:new int[]{128})
//            for (double alpha : new double[]{0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0,1.1,1.2,1.3,1.4/*,1.6,1.8,2.0/*,1.5/*,1.6,1.7,1.8/*,1.9,2.0,2.1,2.2,2.3,2.4*/}) {
//                main.prepareZipf(0,alpha);
//                double[] tmp=main.getMinVAlpha(dataSetSize[0], seriKB*1024);
//                System.out.println("\t\tdataset:\t"+"Zipf"+"\tseriKB:\t"+seriKB+"\talpha:\t"+tmp[1]+"\t"+Math.log(tmp[1])/Math.log(1/0.95)+"\t\tdataSetMinV:\t"+tmp[0]);
//
//            }

        System.out.println("\t\t\tALL_TIME:"+(new Date().getTime()-START_T));
    }
}
