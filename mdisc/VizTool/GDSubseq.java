package mdisc.VizTool;

import java.util.*;

/////////////////////////////

/**
 * subsequence container with the actual string, description length, and number of occurrences
 */
class GDSubseq implements Comparable
{
   public double dl;
   public int count;
   public String string;

   public GDSubseq(String s, int count)
   {
      string = s;
      this.count = count;
   }

   public double calcDL(int N, int nSymbols)
   {
      HashSet<Character> uniq = new HashSet<Character>();
      int len = string.length();
      for(int i = 0; i < len; i++)
         uniq.add(string.charAt(i));
      // dl = Math.log(uniq.size()) + count*Math.log(nSymbols+1);
      dl = len + count + (N - len * count);
      //dl = len * Math.log(uniq.size()) + (count + (N - len * count)) * Math.log(nSymbols + 1);
      return dl;
   }

   public int compareTo(Object o)
   {
      GDSubseq ss = (GDSubseq)o;

      if (dl > ss.dl) return 1;
      if (dl < ss.dl) return -1;

      if (count < ss.count) return 1;
      if (count > ss.count) return -1;

      if (string.length() > ss.string.length()) return 1;
      if (string.length() < ss.string.length()) return -1;

      return string.compareTo(ss.string);
   }
   
   public String toString()
   {
      return String.format("[%s %d  %.2f]", string, count, dl);
   }   
}
