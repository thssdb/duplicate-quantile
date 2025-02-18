import it.unimi.dsi.util.XoRoShiRo128PlusPlusRandomGenerator;
import org.apache.commons.math3.distribution.NormalDistribution;

import java.io.FileWriter;
import java.io.IOException;
import java.util.Arrays;
import java.util.Date;

public class MainForLognormalDataTxt {
    static int N = (int)3.1e7,V=100000;
    static double mu=5.0;
    static String muS="5";
    static double[] a;


    static void prepareLognormal(double sigma) throws IOException {
        if (a == null) a = new double[N+2333+V];
        NormalDistribution nd=new NormalDistribution();
        XoRoShiRo128PlusPlusRandomGenerator rnd=new XoRoShiRo128PlusPlusRandomGenerator(233);
        double[] px=new double[V+1];double totPx=0,maxPx=0;
        int cntN=0;
        for (int i = 1; i <= V; i++){
            px[i]=nd.cumulativeProbability((Math.log(i)-mu)/sigma)-nd.cumulativeProbability((Math.log(i-1)-mu)/sigma);
            maxPx=Math.max(maxPx,px[i]);
            totPx+=px[i];
            int countI=(int)Math.ceil(/*1e-9+*/px[i]*N);
            for(int j=0;j<countI;j++)a[cntN++]=i;
        }
        while(cntN<a.length){
            a[cntN++]=a[rnd.nextInt(cntN-1)+1];
        }
        for(int i=1;i<a.length;i++){// random order
            int p =rnd.nextInt(i);
            double tmp=a[p];
            a[p]=a[i];a[i]=tmp;
        }
        System.out.println("\ttotPs:\t"+totPx+"\t\ttotN:\t"+cntN);
        Arrays.sort(px,1,V+1);
        System.out.println("largest px:\t"+ Arrays.toString(Arrays.copyOfRange(px, V + 1 - 10, V + 1)));
        System.out.println("some a:\t"+ Arrays.toString(Arrays.copyOfRange(a, 0, 100)));
        for(int i=1;i<V;i++)if((int)Math.ceil(/*1e-9+*/px[i]*N)>0){
            System.out.println("smallest V rela_pos:\t"+ 1.0*i/V+"\t\tpx:\t"+px[i]+"\t\ti:\t"+i);
            break;
        }
//        System.out.println(Arrays.toString(px));
    }


    public static void main(String[] args) throws IOException {
        long START = new Date().getTime();
        for (double sigma : new double[]{0.2,0.4,0.6,0.8,1,1.2,1.4,1.6,1.8,2}) {
            prepareLognormal(sigma);
            String filename = "Lognormal3E7Mu"+muS+"Sigma"+(int)(sigma*10+0.1)+".txt";
            System.out.println("\t\t"+filename);
            FileWriter fw = new FileWriter(filename);
            for(int i=0;i<N;i++)fw.write(a[i]+"\n");
            fw.close();
        }
        System.out.println("\t\tALL_TIME:" + (new Date().getTime() - START));
    }
}
