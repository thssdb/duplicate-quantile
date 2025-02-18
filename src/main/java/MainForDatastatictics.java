import it.unimi.dsi.fastutil.doubles.Double2LongMap;
import it.unimi.dsi.fastutil.doubles.Double2LongOpenHashMap;
import it.unimi.dsi.fastutil.doubles.DoubleLongPair;
import it.unimi.dsi.fastutil.doubles.DoubleOpenHashSet;
import it.unimi.dsi.fastutil.objects.ObjectArrayList;

import java.io.BufferedInputStream;
import java.io.BufferedReader;
import java.io.File;
import java.io.FileInputStream;
import java.io.IOException;
import java.io.InputStreamReader;
import java.nio.charset.StandardCharsets;

public class MainForDatastatictics {

    public static void main(String[] args) throws IOException {
        showStatistics();
    }

    private static long dataToLong(double data) {
        long result = Double.doubleToLongBits((double) data);
        return data >= 0d ? result : result ^ Long.MAX_VALUE;
    }
    private static double longToResult(long result)  {
        result = (result >>> 63) == 0 ? result : result ^ Long.MAX_VALUE;
        return Double.longBitsToDouble(result);
    }
    public static void showStatistics() throws IOException {
        String[] filenames = new String[]{
            "Zipf3E7Alpha10.txt","DupliTorqueVoltage.txt",
            "DupliECommercePrice.txt",
            "DupliCustom.txt","DHard.txt"
        };
        int[] Ns = new int[]{
            (int)3e7,(int)3.7e7,(int)6.7e7,(int)1.1e8,(int)3e7
        };
        for (int fileid = 1; fileid < 2; fileid++) {

            double AVG=0,VAR=0,MID=0;
            double midl=0,midr;long rankl=0;
            DoubleOpenHashSet hashSet=new DoubleOpenHashSet();
            Double2LongOpenHashMap hashMap=new Double2LongOpenHashMap();
            int cntN=0;
            for(int t=0;t<2;t++) {
                KLLSketchLazyExactPriori worker=new KLLSketchLazyExactPriori(16*1024*1024);

                File file = new File(filenames[fileid]);
                BufferedInputStream fis = new BufferedInputStream(new FileInputStream(file));
                BufferedReader reader = new BufferedReader(new InputStreamReader(fis, StandardCharsets.UTF_8), 50 * 1024 * 1024);// 50M buffer
                String line = "";
                line = reader.readLine();//ignore first line
                double tmpV, sumV = 0;
                cntN = 0;
                while ((line = reader.readLine()) != null) {
                    tmpV = Double.parseDouble(line);
                    sumV += tmpV;
                    cntN++;
                    if (t==1) VAR += Math.pow(tmpV - AVG, 2);
                    if(t==0) {
                        worker.update(dataToLong(tmpV));
                        hashSet.add(tmpV);
                        hashMap.addTo(tmpV,1);
                    }
                    else{
                        if(tmpV<midl)rankl++;
                        else worker.update(dataToLong(tmpV));
                    }
                    if (cntN > Ns[fileid]) break;
                }
                reader.close();
                fis.close();

                if(t==0){
                    AVG = sumV / Ns[fileid];System.out.println("\t\tavg:"+AVG);
                    double[] res = worker.getFilter(0,0,0,0,Ns[fileid]/2,Ns[fileid]/2,0.9999);
                    midl=res[0];
                    midr=res[1];
                }else{System.out.println("\t\tVAR:"+VAR);
                    double[] res = worker.getFilter(0,0,0,0,Ns[fileid]/2-rankl,Ns[fileid]/2-rankl,0.9999);
                    MID=res[0];
                }
            }
            ObjectArrayList<DoubleLongPair> v2cList=new ObjectArrayList<>();
            for(Double2LongMap.Entry e:hashMap.double2LongEntrySet())
                v2cList.add(DoubleLongPair.of(e.getDoubleKey(),e.getLongValue()));
            v2cList.sort((x,y)->(-Long.compare(x.secondLong(),y.secondLong())));
            System.out.println("\t\tfile:\t"+filenames[fileid]+"\tN:\t"+Ns[fileid]+"\tcntN:"+cntN+"\tDistinctVal:\t"+hashSet.size()+"\tskew:\t"+3*(AVG-MID)/Math.sqrt(VAR/Ns[fileid]));
            System.out.println("\t\tv2cList:\t"+v2cList);
            double percentage=0.01;long perCount=0;
            for(int i=0;i<(int)(v2cList.size()*percentage);i++)perCount+=v2cList.get(i).secondLong();
            System.out.println("\t\t|v2cList|:\t"+v2cList.size()+"\ttotCount of top "+percentage+" values:\t"+perCount);
        }
    }

}