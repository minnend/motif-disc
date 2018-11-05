package mdisc.VizTool;

import java.util.*;

import kdm.data.*;
import kdm.mlpr.suffix_tree.*;
import kdm.util.*;

/** utility functions for dealing with globally discretized sequences */
public class GlobalDisc
{
   public static final char iTERM = (char)50000;
   
   /** @return string containing every unique character in the given strings */
   public static String buildAlphabet(StringRLE[] strings)
   {
      int nSeqs = strings.length;
      HashSet<Character> hashAlpha = new HashSet<Character>();
      for(int i = 0; i < nSeqs; i++){
         int n = strings[i].length();
         for(int j = 0; j < n; j++)
            hashAlpha.add(strings[i].getChar(j));
      }
      Iterator<Character> it = hashAlpha.iterator();
      StringBuffer sb = new StringBuffer();
      while(it.hasNext())
         sb.append(it.next());
      return sb.toString();
   }

   /**
    * @return suffix tree including the RLE strings; terminating characters are >= iTerm; include CountInfo
    *         annotation
    */
   public static SuffixTree buildSufTree(StringRLE[] strings)
   {
      // create suffix tree
      String sAlpha = buildAlphabet(strings);
      SuffixTree sufTree = new SuffixTree();
      TreeBuilder builder = new TreeBuilder(sufTree);
      for(int i = 0; i < strings.length; i++)
         builder.addToken(strings[i].getString() + iTERM);
      CountInfo.annotate(sufTree);
      return sufTree;
   }

   /** @return number of occurrences of each character in RLE or original string */
   public static HashMap<Character, MutableInteger> calcCounts(StringRLE[] strings, boolean bRLE)
   {
      HashMap<Character, MutableInteger> hist = new HashMap<Character, MutableInteger>();
      for(int i = 0; i < strings.length; i++){
         String s = (bRLE ? strings[i].getString() : strings[i].getOrigString());
         int n = s.length();
         for(int j = 0; j < n; j++){
            char c = s.charAt(j);
            MutableInteger mi = hist.get(c);
            if (mi == null){
               mi = new MutableInteger(0);
               hist.put(c, mi);
            }
            mi.inc();
         }
      }
      return hist;
   }
   
   /** @return number of occurrences of each ngram in RLE or original strings */
   public static HashMap<String,MutableInteger> calcNGramCounts(int ngram, StringRLE[] strings, boolean bRLE)
   {
      HashMap<String, MutableInteger> hist = new HashMap<String, MutableInteger>();
      for(int i = 0; i < strings.length; i++){
         String s = (bRLE ? strings[i].getString() : strings[i].getOrigString());         
         int n = s.length()-ngram;
         for(int j = 0; j <= n; j++){
            String ss = s.substring(j, j+ngram);
            MutableInteger mi = hist.get(ss);
            if (mi == null){
               mi = new MutableInteger(0);
               hist.put(ss, mi);
            }
            mi.inc();
         }
      }
      return hist;
   }

   /** @return probability of each character given strings */
   public static HashMap<Character, Double> calcHist(StringRLE[] strings, boolean bRLE)
   {
      HashMap<Character, Double> hist = new HashMap<Character, Double>();
      HashMap<Character, MutableInteger> counts = calcCounts(strings, bRLE);

      // calc sum
      int sum = 0;
      Iterator<Character> it = counts.keySet().iterator();
      while(it.hasNext())
         sum += counts.get(it.next()).getValue();

      // now calc frequencies
      it = counts.keySet().iterator();
      while(it.hasNext()){
         char c = it.next();
         double v = (double)counts.get(c).getValue() / sum;
         hist.put(c, v);
      }

      return hist;
   }
   
   /** @return probability of each ngram given strings */
   public static HashMap<String, Double> calcNGramHist(int ngram, StringRLE[] strings, boolean bRLE)
   {
      HashMap<String, Double> hist = new HashMap<String, Double>();
      HashMap<String, MutableInteger> counts = calcNGramCounts(ngram, strings, bRLE);

      // calc sum
      int sum = 0;
      Iterator<String> it = counts.keySet().iterator();
      while(it.hasNext())
         sum += counts.get(it.next()).getValue();

      // now calc frequencies
      it = counts.keySet().iterator();
      while(it.hasNext()){
         String ss = it.next();
         double v = (double)counts.get(ss).getValue() / sum;
         hist.put(ss, v);
      }

      return hist;
   }

   /** @return log-likelihood of string given histogram */
   public static double calcLogLik(String s, HashMap<Character, Double> hist)
   {
      double loglik = Library.LOG_ONE;
      int n = s.length();
      for(int i = 0; i < n; i++)
         loglik += Math.log(hist.get(s.charAt(i))); // TODO precalc log(hist)
      return loglik;
   }

   /** @return log-likelihood of string given histogram */
   public static double calcNGramLogLik(String s, HashMap<String, Double> hist)
   {      
      int ngram = hist.keySet().iterator().next().length();
      double loglik = Library.LOG_ONE;
      int n = s.length() - ngram;
      for(int i = 0; i <= n; i++){
         String ss = s.substring(i, i+ngram);
         loglik += Math.log(hist.get(ss)); // TODO precalc log(hist)
      }
      return loglik;
   }

   
   /** @return list of motifs in data */
   public static ArrayList<GDSubseq> findMotifs(StringRLE[] strings)
   {
      ArrayList<GDSubseq> motifs = new ArrayList<GDSubseq>();
    
      String sAlpha = buildAlphabet(strings);
      SuffixTree sufTree = buildSufTree(strings);
      collectSubseqs(sufTree, sufTree.getRoot(), "", motifs);
      
      return motifs;
   }
   
   /**
    * Collect all of the subsequences at the nodes of the suffix tree. note: subsequences that include a
    * terminating character will not be included.
    * 
    * @param tree annotated (with CountInfo) suffix tree to search
    * @param node root node at wich to start search
    * @param prefix prefix so far
    * @param list list in which to collect the results
    */
   protected static void collectSubseqs(SuffixTree tree, NodeInterface node, String prefix,
         ArrayList<GDSubseq> list)
   {
      String label = tree.getSubstring(node.getLeftIndex(), node.getLength());
      if (label.length() > 0){
         char end = label.charAt(label.length() - 1);
         if (end >= iTERM) return;
      }

      if (node instanceof InternalNode){
         InternalNode inode = (InternalNode)node;

         // if this isn't the root, store the subseq
         if (node != tree.getRoot()){
            String s = prefix + label;
            CountInfo ci = (CountInfo)inode.getInfo(CountInfo.class);
            list.add(new GDSubseq(s, ci.getCount()));
         }

         // recurse through all child nodes
         NodeInterface kid = inode.getFirstChild();
         while(kid != null){
            collectSubseqs(tree, kid, prefix + label, list);
            kid = kid.getRightSibling();
         }
      }
      else{ // we have a leaf node
         String s = prefix + label;
         CountInfo ci = (CountInfo)node.getInfo(CountInfo.class);
         list.add(new GDSubseq(s, ci.getCount()));
      }
   }
   
   public static void main(String[] args)
   {
      //String s = "aabcdxyzabcdzyxbcdxaabccdxyxyzcdcdabcd";
      String s = "abcdexabadeyabedezabcdcqa";
      StringRLE rle = new StringRLE(s, 'a');
      StringRLE[] strings = new StringRLE[]{ rle };
      
      for(int i=0; i<strings.length; i++)
         System.err.printf("%s\n", strings[i].getString());
      System.err.println();
                                   
      ArrayList<GDSubseq> motifs = findMotifs(strings);
      
      System.err.printf("p(char):\n");
      HashMap<Character, Double> hist = calcHist(strings, true);
      for(GDSubseq motif : motifs)
         motif.dl = calcLogLik(motif.string, hist) * motif.count;
      Collections.sort(motifs);
      for(GDSubseq motif : motifs)
         System.err.printf("%s\n", motif);
      System.err.println();
      
      System.err.printf("p(bigram):\n");
      HashMap<String,Double> bihist = calcNGramHist(2, strings, true);
      for(GDSubseq motif : motifs)
         motif.dl = calcNGramLogLik(motif.string, bihist) * motif.count;
      Collections.sort(motifs);
      for(GDSubseq motif : motifs)
         System.err.printf("%s\n", motif);
      System.err.println();
   }
}

