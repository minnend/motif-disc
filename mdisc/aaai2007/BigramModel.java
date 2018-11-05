package mdisc.aaai2007;

import java.util.*;
import java.io.*;
import kdm.data.*;
import kdm.io.*;
import kdm.util.*;

/** learn an n-gram model from a set of markup sets */
public class BigramModel
{
   protected LabelEditor led;
   protected HashMap<String, Integer> tag2i;
   protected String[] tags;
   protected int[][] tran;

   public BigramModel(LabelEditor led)
   {
      this.led = led;
   }

   /** @return set of tags used in the markup sets */
   public static Set<String> findTags(ArrayList<MarkupSet> labels, LabelEditor led)
   {
      TreeSet<String> uniq = new TreeSet<String>();
      for(MarkupSet marks : labels){
         int n = marks.size();
         for(int i = 0; i < n; i++){
            String s = marks.get(i).getTag();
            if (led != null){
               if (!led.keepLabel(s)) continue;
               s = led.adjustLabel(s);
            }
            uniq.add(s);
         }
      }
      return uniq;
   }

   public int getIndex(String tag)
   {
      return tag2i.get(tag);
   }

   public String getTag(int ix)
   {
      return tags[ix];
   }

   /** tally the bigram counts for the given markup sets */
   public boolean learn(ArrayList<MarkupSet> labels)
   {
      Set<String> tagSet = findTags(labels, led);
      int nTags = tagSet.size();
      tags = new String[nTags];
      tag2i = new HashMap<String, Integer>();
      Iterator<String> it = tagSet.iterator();
      for(int i = 0; i < nTags; i++){
         tags[i] = it.next();
         tag2i.put(tags[i], i);
         System.err.printf("%d) %s\n", i+1, tags[i]);
      }
      tran = new int[nTags][nTags];

      for(MarkupSet marks : labels){
         marks.sort();
         int nMarks = marks.size();
         for(int i = 1; i < nMarks; i++){
            TimeMarker tm1 = marks.get(i - 1);
            TimeMarker tm2 = marks.get(i);
            String s1 = tm1.getTag();
            String s2 = tm2.getTag();
            if (led != null){
               if (!led.keepLabel(s1)) continue;
               if (!led.keepLabel(s2)){
                  i++; // no need to search this again in the tm1 spot
                  continue;
               }
               s1 = led.adjustLabel(s1);
               s2 = led.adjustLabel(s2);
            }

            int ix1 = tag2i.get(s1);
            int ix2 = tag2i.get(s2);
            tran[ix1][ix2]++;
         }
      }
      return true;
   }

   /** @return normalized bigram map (rows sum to one) */
   public double[][] getNormMap()
   {
      int n = tran.length;
      double[][] p = new double[n][n];
      for(int i = 0; i < n; i++){
         int sum = Library.sum(tran[i]);
         for(int j = 0; j < n; j++)
            p[i][j] = (double)tran[i][j] / sum;
      }
      return p;
   }

   public boolean genGraphViz(File file, double thresh, boolean bSelf)
   {
      if (thresh < 0) thresh = 0;
      try{
         PrintWriter out = new PrintWriter(new BufferedWriter(new FileWriter(file)));
         
         double[][] p = getNormMap();
         int n = p.length;
         
         out.println("digraph G {");
         for(int i=0; i<n; i++)
            out.printf(" %d [label=\"%s\"];\n", i, tags[i]);         
         
         for(int i=0; i<n; i++){
            for(int j=0; j<n; j++){
               if (!bSelf && i==j) continue;
               if (p[i][j]<=thresh) continue;
               out.printf(" %d -> %d [label=\"%d\"];\n", i, j, Math.round(p[i][j]*100));
            }
         }
         out.println("}");
         
         out.close();
         return true;
      } catch (IOException e){
         return false;
      }
   }
   
   public static void main(String[] args)
   {
      //String sFile = "/home/dminn/research/kdm/data/asl-joseph/aslxyd.def";
      //String sFile = "/home/dminn/research/kdm/data/ah-data/exercise/sub8/accgyr.def";
      String sFile = "/home/dminn/asl-ijcai.def";
      ArrayList<MarkupSet> labels = LabeledDataLoader.loadLabels(new File(sFile));
      if (labels==null || labels.isEmpty()){
         System.err.printf("Error: unable to find any label data\n (%s)\n", sFile);
         System.exit(1);
      }
      System.err.printf("Found %d label files.\n", labels.size());
      
      LabelEditor led = new LabelEditor()
      {

         public String adjustLabel(String sLabel)
         {
            return sLabel.replace(" ", "\\n");
         }

         public boolean keepLabel(String sLabel)
         {
            return true;
         }
         
      };
      
      BigramModel bgm = new BigramModel(led);
      bgm.learn(labels);
      File fGraph = new File("/tmp/graph.dot");
      bgm.genGraphViz(fGraph, 0.1, true);
   }
}
