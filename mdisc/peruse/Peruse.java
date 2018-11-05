package mdisc.peruse;

import gnu.getopt.*;
import kdm.tools.SupTest;
import kdm.util.*;
import kdm.io.*;
import kdm.io.DataLoader.DLRaw;
import kdm.io.DataLoader.DataLoader;
import kdm.data.*;
import kdm.data.transform.*;
import kdm.models.*;
import kdm.models.misc.*;
import kdm.mlpr.*;

import org.apache.commons.math.stat.*;
import java.util.*;
import java.util.zip.*;
import java.io.*;
import java.util.regex.*;

import static kdm.util.TTResult.*;
import static mdisc.peruse.Params.*;

/** Control program for PERUSE style discovery */
public class Peruse
{
   public static PrintStream output = System.out;
   public static boolean bVerbose = false;
   public static DOGui dogui = null;
   public static Params params;

   // //////////////////////////////////////////////////////////

   public static void usage()
   {
      System.err.println("USAGE: java kdm.tools.DiscoverOates [options] <file patterns>");
      System.err.println();
      System.err.println(" File Patterns:");
      System.err.println("  At least one file pattern is required.  Each pattern can be a specific");
      System.err.println("  file or a wildcard expression (e.g., *dat).");
      System.err.println();
      System.err.println(" Options:");
      System.err.println("  -v                      Output verbose information");
      System.err.println("  -help                   Print this message");
      System.err.println("  -path <base path>       Specify a base path (.)");
      System.err.println("  -nlen <n>               Specify a minimum pattern length (30)");
      System.err.println("  -nskip <n>              Specify a skip size for the sliding window in ");
      System.err.println("                           exhaustive search mode (10)");
      System.err.println("  -ngrow <n>              Specify a grow length for extending a pattern (5)");
      System.err.println("  -nminmatch <n>          Specify minimum number of pattern occurrences (3)");
      System.err.println("  -nmaxinst <n>           Specify maximum number of pattern occurrences");
      System.err.println("                           to find per pattern (no max)");
      System.err.println("  -nmarkers <n>           Specify number of markers for gibbs sampler (3)");
      System.err.println("  -atime <alpha>          Specify alpha-temporal value (0.05)");
      System.err.println("  -aseries <alpha>        Specify alpha-series value (0.05)");
      System.err.println("  -nsamples <n>           Specify number of samples for stop search test (1000)");
      System.err.println("  -nmaxpats <n>           Specify maximum number of pats to find (no max)");
      System.err.println("  -out <file>             Write results to given file (stdout)");
      System.err.println("  -search <type>          Specify type of search to find exemplars (exhaustive)");
      System.err.println("                           given_i.j - use given seq (i) and location (j)");
      System.err.println("                           exhaustive - try all possible sites");
      System.err.println("                           exone - try all possible sites in first sequence");
      System.err.println("                           uniform - uniformly sample locations");
      System.err.println("                           suftree - init using a suffix tree");
      System.err.println("  -presplit <w:t:g:o>     split sequences around gaps of zero variance");
      System.err.println("                           w = length of window for computing variance");
      System.err.println("                           t = threshold for equiv to zero");
      System.err.println("                           g = min length of gap for split");
      System.err.println("                           o = max length of positive gap to overlook");
      System.err.println("  -ngarbage <n>           # random subseqs to use to form garbage model (0)");
      System.err.println("  -model <model>          Specify model to use (hmm or oates)");
      System.err.println("  -hmmeval <type>         Use \"viterbi\" (def) or \"forward\" algorithm?");
      System.err.println("  -hmmtrain <type>        Use \"bw\" (def) or \"viterbi\" algorithm?");
      System.err.println("  -hmmstates <n>          Number of states in each HMM (def: 5)");
      System.err.println("  -hmmskip <n>            Number of skip states in HMMs (def: 1)");
      System.err.println("  -gui                    Show a GUI to track progress (def: no gui)");
      System.err.println("  -fixedvar               Don't re-estimate model variance");
      System.err.println();
   }

   /**
    * Dumps data from the given nBestOMaps structure to stderr.
    */
   public static void dumpPattern(OatesMapping omap, Sequence seq)
   {
      if (seq.hasParent()) System.err.printf(" OMap: %d.%d %d (%d)  %.2f\n", seq.getParentIndex() + 1, seq
            .getParentOffset()
            + omap.getFirstDataIndex(), omap.getDataLength(), omap.getPatternLength(), omap.score);
      else System.err.printf(" OMap: %d.%d %d [%d]  %.2f\n", omap.iSeries + 1, omap.getFirstDataIndex(),
            omap.getDataLength(), omap.getPatternLength(), omap.score);
   }

   /**
    * Dumps data from the given nBestOMaps structure to stderr.
    */
   public static void dumpPattern(NBestOMaps nBestMaps, ArrayList<Sequence> vdata)
   {
      for(int i = 0; i < nBestMaps.size(); i++){
         OatesMapping omap = nBestMaps.get(i);
         Sequence seq = vdata.get(omap.iSeries);
         dumpPattern(omap, seq);
      }
   }

   public static OatesModelUSamp buildGarbageModel(ArrayList<Sequence> vdata, int nSites)
   {
      ArrayList<Sequence> vtrain = new ArrayList<Sequence>(nSites);

      // initialize the list of valid start points
      int n = vdata.size();
      SpanList span[] = new SpanList[n];
      for(int i = 0; i < n; i++)
         span[i] = new SpanList(0, vdata.get(i).length() - params.nwLen, true);

      // find the required number of subseqs
      for(int iSite = 0; iSite < nSites; iSite++){
         int i = pickRandomSeq(vdata, span);
         int x = (int)(Library.random() * span[i].size());
         assert (x >= 0 && x < span[i].size()) : "x=" + x + "   span=" + span[i].size() + "  seqlen="
               + vdata.get(i).length() + "  nwlen=" + params.nwLen;
         int j = span[i].get(x);
         assert (j >= 0 && j < vdata.get(i).length());
         // System.err.printf("iSite=%d/%d i=%d/%d x=%d/%d j=%d/%d\n",
         // iSite, nSites, i, vdata.size(), x, span[i].size(), j, vdata.get(i).length());
         span[i].sub(j, j + params.nwLen - 1);
         vtrain.add(vdata.get(i).subseq(j, j + params.nwLen));
      }

      // now train the model
      OatesModelUSamp om = new OatesModelUSamp(vtrain, params.nwLen, params.initv, params.minv);
      return om;
   }

   /**
    * Search for areas of no change (near zero variance) and remove them
    */
   public static void presplit(ArrayList<Sequence> vdata)
   {
      System.err.println("Before presplit: " + vdata.size() + " seqs");
      ArrayList<WindowLocation> vLocs = new ArrayList<WindowLocation>();
      for(int iSeq = 0; iSeq < vdata.size(); iSeq++){
         Sequence seq = vdata.get(iSeq);

         // slide window and calc variance
         int n = seq.length();
         int ndims = seq.getNumDims();
         // int step = (int)Math.max(1, params.splitW / 4);
         SpanList span = new SpanList(0, n - 1, false);

         System.err.println("seq: " + iSeq); // TODO
         double y[] = seq.extractDim(0, 0, 5);
         for(int i = 0; i < 5; i++)
            System.err.printf("%.4f ", y[i]);
         System.err.println();

         int a;
         boolean bZero = false;
         for(a = 0; a + params.splitW <= n; a++){
            // all dimensions must have zero variance
            bZero = true;
            for(int j = 0; j < ndims; j++){
               double[] x = seq.extractDim(j, a, params.splitW);
               if (j == 0 && a == 275){
                  System.err.printf("%d.%d:\n", iSeq + 1, a);
                  for(int q = 0; q < x.length; q++)
                     System.err.printf("%.4f ", x[q]);
                  System.err.println();
               }
               double v = StatUtils.variance(x);
               if (v > params.splitT) // variance too large?
               {
                  bZero = false;
                  break;
               }
            }
            if (bZero) span.add(a + params.splitW / 2);
         }

         // add beginning and end if appropriate
         if (span.contains(params.splitW / 2)) span.add(0, params.splitW / 2);
         if (bZero) span.add(a - 1 + params.splitW / 2, n);

         // "close" holes
         if (params.splitO > 0){
            System.err.println("Closing holes (" + params.splitO + ")");
            span.suffix(params.splitO);
            span.suffix(-params.splitO);
         }

         // now we can search for long zero spans
         int iStart = 0;
         n = span.getNumSpans();
         for(int i = 0; i < n; i++){
            Range r = span.getRange(i);
            if (r.length() >= params.splitG){
               vLocs.add(new WindowLocation(iSeq, r.a, r.length()));
               iStart = r.b + 1;
            }
         }
      }

      // dump matlab file for cut segments
      try{
         PrintWriter out = new PrintWriter(new FileWriter("presplit.txt"));
         for(int i = 0; i < vLocs.size(); i++){
            WindowLocation loc = vLocs.get(i);
            out.printf("\"%d\" %d %d\n", loc.iSeries + 1, loc.iStart, loc.iStart + loc.length() - 1);
         }
         out.close();
      } catch (IOException ioe){
         ioe.printStackTrace();
      }

      // we need to make a global change, so we can't just assign vsplit to vdata
      ArrayList<Sequence> vSplit = Sequence.chop(vdata, vLocs, params.nwLen);
      vdata.clear();
      vdata.addAll(vSplit);
   }

   /** @return index of random sequence considering number of free spots in each seq */
   public static int pickRandomSeq(ArrayList<Sequence> vdata, SpanList span[])
   {
      int n = vdata.size();
      int w[] = new int[n];
      int wsum = 0;
      for(int i = 0; i < n; i++){
         w[i] = span[i].size();
         wsum += w[i];
      }

      int x = (int)(Library.random() * wsum);

      int i = 0;
      int cdf = w[0];
      while(x >= cdf){
         i++;
         cdf += w[i];
      }
      assert (i >= 0 && i < n);
      return i;
   }

   /**
    * Determines if a proposed temporal extension is supported by the data.
    */
   public static boolean acceptGrowth(OatesModelUSamp om, ArrayList<Sequence> vdata, NBestOMaps nBestMaps,
         int nRealGrow)
   {
      int newStart, newStop, newLen;
      int oldStart, oldStop, oldLen;
      int wlength = om.size();
      if (nRealGrow < 0){
         newStart = 0;
         newStop = -nRealGrow;
         oldStart = newStop;
         oldStop = wlength;
      }
      else{
         oldStart = 0;
         oldStop = wlength - nRealGrow;
         newStart = oldStop;
         newStop = wlength;
      }
      newLen = newStop - newStart;
      oldLen = oldStop - oldStart;

      int nMaps = nBestMaps.size();
      GaussianDyn1D gOld = new GaussianDyn1D();
      GaussianDyn1D gNew = new GaussianDyn1D();
      for(int i = 0; i < nMaps; i++){
         OatesMapping omap = nBestMaps.get(i);
         assert (omap.getPatternLength() == om.size());
         Sequence seq = vdata.get(omap.iSeries);
         for(int j = 0; j < omap.getPatternLength(); j++){
            int iPat = j;
            int iDat = omap.imap[j];
            double v = om.eval(iPat, seq.get(iDat));
            if (iPat >= oldStart && iPat < oldStop) gOld.add(v, false);
            else gNew.add(v, false);
         }
      }
      gOld.update();
      gNew.update();

      double p = Library.tutest(gOld.getMean(), gNew.getMean(), gOld.getVar(), gNew.getVar(), gOld.getN(),
            gNew.getN(), Tails.One).p;
      boolean bAccept = (gNew.getMean() >= gOld.getMean() || p >= params.aTemporal);
      System.err.println(" accept: old=" + Library.df.format(gOld.getMean()) + "  new="
            + Library.df.format(gNew.getMean()) + "  p=" + Library.df.format(p) + "  accept=" + bAccept);

      return bAccept;
   }

   protected static ModelWindow getModelFromExemplar(ArrayList<Sequence> vdata, WindowLocation win)
   {
      if (params.model == MODEL_OATES){
         OatesModelUSamp om = new OatesModelUSamp(vdata.get(win.iSeries), win, params.initv, params.minv);

         // find the min_matches series with the best matches
         // TODO more than one mapping per seq
         NBestList<Pair<Double, Integer>> trainList = new NBestList<Pair<Double, Integer>>(params.nMinMatches);
         trainList.clear();
         for(int iToSeries = 0; iToSeries < vdata.size(); iToSeries++){
            Sequence seqTo = vdata.get(iToSeries);
            OatesMapping omap = om.findBestMapping(seqTo);
            omap.iSeries = iToSeries;
            // System.err.println("omap: "+omap);
            trainList.add(new Pair(omap.score, iToSeries));
         }
         if (trainList.size() < params.nMinMatches) return null; // must have min num of matches

         // build the vtrain list and train the model
         ArrayList<Sequence> vtrain = new ArrayList<Sequence>(params.nMinMatches);
         for(int i = 0; i < params.nMinMatches; i++)
            vtrain.add(vdata.get(trainList.get(i).second));

         om.trainEM(vtrain, params.bUpdateVar);

         // find the best mapping for each training sequence
         // TODO could find multiple mappings in one sequence
         double score = 0.0;
         NBestOMaps nBestMaps = new NBestOMaps(params.nMinMatches);
         for(int i = 0; i < params.nMinMatches; i++){
            int iSeries = trainList.get(i).second;
            OatesMapping omap = om.findBestMapping(vdata.get(iSeries));
            omap.iSeries = iSeries;
            score += omap.score;
            nBestMaps.add(omap);
         }

         OatesWindow bestWindow = new OatesWindow(win.iSeries, win.start(), params.nwLen, score, nBestMaps);
         return new ModelWindow(om, bestWindow);
      }
      else if (params.model == MODEL_HMM){
         Sequence seq = vdata.get(win.iSeries);
         Sequence sub = seq.subseq(win.start(), win.end());
         HmmLR hmm = new HmmLRa(params.nHmmStates, params.nHmmSkip, seq.getNumDims());
         hmm.setUpdateVar(false); // can't update var from single example
         hmm.setVar(params.initv);
         hmm.init_segk(sub);
         if (params.hmmTrain == HMM_VITERBI) hmm.train_viterbi(sub);
         else hmm.train_bw(sub);

         // find the best min_matches matches
         // TODO should find multiple matches per sequence
         double score = 0.0;
         NBestOMaps nBestMaps = new NBestOMaps(params.nMinMatches);
         Range rSubLen = new Range(params.nwLen, (int)Math.round(params.nwLen * 1.25));
         int nStep = Library.max(params.nwSkip / 2, 1);
         for(int iToSeries = 0; iToSeries < vdata.size(); iToSeries++){
            Sequence seqTo = vdata.get(iToSeries);
            // ScoredWindow sw = hmm.findBestSubseqSlow(seqTo, null, true, rSubLen, nStep);
            ScoredWindow sw = hmm.findBestSubseq(seqTo, true);
            if (sw == null) continue;
            OatesMapping omap = new OatesMapping(sw.score, 2, iToSeries);
            omap.imap[0] = sw.start();
            omap.imap[1] = sw.end();
            nBestMaps.add(omap);
            score += omap.score;
         }
         OatesWindow bestWindow = new OatesWindow(win.iSeries, win.start(), params.nwLen, score, nBestMaps);
         return new ModelWindow(hmm, bestWindow);
      }

      assert false : "unknown model type";
      return null;
   }

   /** find the best window within the given range */
   protected static OatesWindow findBestWindow(ArrayList<Sequence> vdata, WindowLocation wloc)
   {
      System.err.printf("find best window: %d.%d -> %d\n", wloc.iSeries, wloc.iStart, wloc.getLastIndex());
      Sequence seqFrom = vdata.get(wloc.iSeries);
      int len = seqFrom.length() - params.nwLen + 1;
      OatesWindow bestWindow = new OatesWindow();
      ArrayList<Sequence> vtrain = new ArrayList<Sequence>(params.nMinMatches);

      // loop through all window locations
      for(int iw = wloc.iStart; iw <= wloc.getLastIndex() && iw < len; iw += params.nwSkip){
         TimerMS timer = new TimerMS();
         ModelWindow mw = getModelFromExemplar(vdata, new WindowLocation(wloc.iSeries, iw, params.nwLen));

         // update the best window if necessary
         if (mw.window.score > bestWindow.score){
            bestWindow.set(wloc.iSeries, iw, params.nwLen, mw.window.score, mw.window.maps);
            System.err.println("new best: " + bestWindow);
         }

         System.err.println("searching -- iw=" + (wloc.iSeries + 1) + "." + iw + " / " + wloc.getLastIndex()
               + "  (" + timer.time() + "ms)");
      }
      return (bestWindow.loc.iSeries < 0 ? null : bestWindow);
   }

   /**
    * Finds the best exemplar window in the given list of sequences by sliding a window of the specified
    * length and skip size.
    */
   protected static OatesWindow findBestWindow(ArrayList<Sequence> vdata)
   {
      // we need to get a base exemplar that can be grown
      TimerMS timer = new TimerMS();
      OatesWindow bestWindow = null;

      if (params.kSearch == Params.SEARCH_GIVEN)
      {
         int iSeq = params.given.iSeries;
         int iStart = params.given.iStart;
         assert(vdata.get(iSeq).getParentIndex() == iSeq);
         bestWindow = findBestWindow(vdata, new WindowLocation(iSeq, iStart, 1));         
      }
      else if (params.kSearch == Params.SEARCH_EXHAUST_ONE){
         int iOrigSeq = 0; // the orig index of the sequence to search

         for(int iSeq = 0; iSeq < vdata.size(); iSeq++){
            if (vdata.get(iSeq).getParentIndex() != iOrigSeq) continue;
            OatesWindow owin = findBestWindow(vdata, new WindowLocation(iSeq, 0, vdata.get(iSeq).length()));
            if (bVerbose) System.err.println("Best local window (exone): " + owin);
            if (bestWindow == null || owin.score > bestWindow.score){
               bestWindow = owin;
               if (dogui != null) dogui.setBestWindow(bestWindow, vdata);
            }
         }
      }
      else if (params.kSearch == Params.SEARCH_EXHAUSTIVE){
         for(int iSeq = 0; iSeq < vdata.size(); iSeq++){
            OatesWindow owin = findBestWindow(vdata, new WindowLocation(iSeq, 0, vdata.get(iSeq).length()));
            if (bVerbose) System.err.println("Best local window (exone): " + owin);
            if (bestWindow == null || owin.score > bestWindow.score){
               bestWindow = owin;
               if (dogui != null) dogui.setBestWindow(bestWindow, vdata);
            }
         }
      }
      else{
         System.err.printf("Error: unrecognized exemplar search strategy (%s)\n", params.kSearch);
         System.exit(1);
      }
      return bestWindow;
   }

   /**
    * Given a set of patterns and a set of series that don't contain those patterns (ie, the patterns have
    * been chopped out of the data), determines the next best pattern location if one is supported by the data
    * for the given alpha value.
    * 
    * @param vPatsN N best patterns found so far
    * @param vdata time series data to be searched
    * @param aSeries alpha_{series} value as defined in PERUSE paper
    * @param nSamples number of samples to draw from N+1 model for comparison
    * @return OatesMapping of next best pattern location or null if none is supported
    */
   protected static OatesMapping findNextPatternSpot(ArrayList<Sequence> vPatsN, ArrayList<Sequence> vdata,
         OatesModelUSamp omGarbage, double aSeries, int nSamples)
   {
      TimerMS timer = new TimerMS();
      OatesMapping nextPatMap = null;

      if (params.model == MODEL_OATES){
         // train a model on the existing patterns
         System.err.print("fnps) Training base model (Oates)... ");
         OatesModelUSamp omBase = new OatesModelUSamp(vPatsN, params.nwLen, params.initv, params.minv);
         System.err.println("done (" + timer.time() + "ms) [len=" + omBase.length() + "].");
         int np = vPatsN.size();
         int np1 = np + 1;
         int nMaxSmaller = (int)(nSamples * (1.0 - aSeries));

         /*
          * try{ // TODO PrintWriter out = new PrintWriter(new FileWriter(String.format("model%d.txt", np)));
          * out.print(omBase); out.close(); }catch(IOException e){e.printStackTrace();}
          */

         // find the next best pattern location
         System.err.print("Searching for best match... ");
         timer.reset();
         for(int i = 0; i < vdata.size(); i++){
            OatesMapping omap = omBase.findBestMapping(vdata.get(i));
            if (omap != null){
               omap.iSeries = i;
               if (nextPatMap == null || omap.score > nextPatMap.score) nextPatMap = omap;
            }
         }
         System.err.printf("done (%dms).\n", timer.time());
         if (nextPatMap == null) return null; // no more valid spots

         System.err.println("testing pat loc: " + nextPatMap);
         Sequence nextPat = vdata.get(nextPatMap.iSeries).subseq(nextPatMap.getLoc());

         // if we have a garbage model, compare against that
         if (omGarbage != null){
            double xsum = 0.0, ysum = 0.0;
            for(int i = 0; i < np; i++){
               double x = omBase.eval(vPatsN.get(i));
               double y = omGarbage.eval(vPatsN.get(i));

               xsum += x;
               ysum += y;
            }
            double xmean = xsum / (double)np;
            double ymean = ysum / (double)np;

            double x = omBase.eval(nextPat);
            double y = omGarbage.eval(nextPat);

            double xr = x / xmean;
            double yr = y / ymean;

            System.err.printf("model: %.2f / %.2f = %.2f    garbage: %.2f / %.2f = %.2f\n", x, xmean, xr, y, ymean, yr);
            if (yr < xr) return null; // TODO
         }
         else{
            // build a list containing N+1 best patterns
            ArrayList<Sequence> vPatsNp1 = new ArrayList<Sequence>(vPatsN);
            vPatsNp1.add(nextPat);

            // train a model on the N+1 patterns
            OatesModelUSamp omNew = new OatesModelUSamp(vPatsNp1, params.nwLen);

            /*double sumA = 0, sumB = 0, lastA = 0, lastB = 0;
            for(int i = 0; i < np1; i++){
               double a = omBase.eval(vPatsNp1.get(i));
               double b = omNew.eval(vPatsNp1.get(i));
               
               sumA += a;
               sumB += b;
               if (i == np1 - 1){
                  lastA = a;
                  lastB = b;
               }
            }
            int na = omBase.length();
            int nb = omNew.length();
            System.err.printf("sum (avg): %.2f (%.2f, %.2f) | %.2f (%.2f, %.2f)\n", sumA / na, sumA
                  / (np1 * na), (sumA - lastA) / ((np1 - 1) * na), sumB / nb, sumB / (np1 * nb),
                  (sumB - lastB) / ((np1 - 1) * nb));
             */
            // TODO better to sample from N model and test new pattern location?

            // find base pattern with smallest prob given N+1 model
            double pSmall = omNew.eval(vPatsN.get(0));
            for(int i = 1; i < np; i++){
               double p = omNew.eval(vPatsN.get(1));
               if (p < pSmall) pSmall = p;
            }

            // test this prob against samples from N+1 model to see if it's too small
            int nSmaller = 0;
            for(int i = 0; i < nSamples; i++){
               Sequence sample = omNew.sample();
               double p = omNew.eval(sample);
               if (pSmall < p) nSmaller++;
            }
            System.err.printf("nSmaller: %d (max: %d)  pSmall=%f\n", nSmaller, nMaxSmaller, pSmall);

            // if too small, then there is no next pattern
            if (nSmaller >= nMaxSmaller) return null;
         }
      }
      else if (params.model == MODEL_HMM){
         // train a model on the existing patterns
         /*
          * System.err.print("fnps) Training base model (HMM)... "); timer.reset(); OatesModelUSamp omBase =
          * new OatesModelUSamp(vPatsN, params.nwLen, initv, minv); System.err.println("done
          * ("+timer.time()+"ms) [len="+omBase.length()+"]."); int np = vPatsN.size(); int np1 = np + 1; int
          * nMaxSmaller = (int)(nSamples*(1.0 - aSeries)); // find the next best pattern location
          * System.err.print("Searching for best match... "); PVMJobSpool spool = new PVMJobSpool(jpvm, tids);
          * for(int i=0; i<vdata.size(); i++) { fbm_req req = new fbm_req(i, omBase); spool.add(req,
          * MSG_PVM_FBM_REQ); } while(!spool.empty()) { msg = jpvm.pvm_recv(MSG_PVM_FBM_RET);
          * spool.pop(getWorkerIndex(msg.sourceTid)); fbm_ret fbm = new fbm_ret(msg.buffer); OatesMapping omap =
          * fbm.getMapping(); //if (bVerbose) System.err.println("Best local mapping: "+omap); if (omap!=null &&
          * (nextPatMap==null || omap.score > nextPatMap.score)) nextPatMap = omap; } System.err.println("done
          * ("+timer.time()+"ms)."); if (nextPatMap == null) return null; // no more valid spots
          * 
          * System.err.println("testing pat loc: "+nextPatMap); Sequence nextPat =
          * vdata.get(nextPatMap.iSeries).subseq(nextPatMap.getLoc());
          * 
          * double xsum = 0.0, ysum = 0.0; for(int i=0; i<np; i++) { double x = omBase.eval(vPatsN.get(i));
          * double y = omGarbage.eval(vPatsN.get(i));
          * 
          * xsum += x; ysum += y;
          * 
          * System.err.printf("%d: %.2f %.2f\n", i+1, x, y); } double xmean = xsum / (double)np; double ymean =
          * ysum / (double)np;
          * 
          * double x = omBase.eval(nextPat); double y = omGarbage.eval(nextPat);
          * 
          * double xr = x/xmean; double yr = y/ymean;
          * 
          * System.err.printf("model: %.2f / %.2f %.2f / %.2f\n", x, xr, y, yr); if (yr < xr) return null;
          */
      }
      else{
         assert false : "unknown model type: " + params.model;
         return null;
      }

      // past all the tests, so return the next pattern occurrence
      return nextPatMap;
   }

   /**
    * Finds significant, recurring patterns in a set of data sequences.
    * 
    * @param iPattern zero-based index of the pattern to find
    * @param vdata list of data sequences
    * @return true if a pattern was found, false otherwise
    */
   public static boolean discover(int iPattern, ArrayList<Sequence> vdata, OatesModelUSamp omGarbage)
   {
      System.err.printf("Discover %d: %d seqs\n", iPattern + 1, vdata.size());
      if (vdata.size() < params.nMinMatches) return false;

      // find the best window to seed the pattern
      OatesWindow bestWindow = findBestWindow(vdata);
      System.err.println("best window (raw) = " + bestWindow);
      if (bestWindow == null) return false;
      dumpPattern(bestWindow.maps, vdata);

      // TODO !! grow the sequence
      /*
       * if (params.nwGrow > 0){ // now that we have a base exemplar (the best window), we can grow it
       * temporally // first we'll try to extend the end of the example while(bestWindow.loc.end() <
       * vdata.get(bestWindow.loc.iSeries).length()){ WindowLocation newWin = new
       * WindowLocation(bestWindow.loc.iStart, Math.min(bestWindow.loc .length() + params.nwGrow,
       * vdata.get(bestWindow.loc.iSeries).length() - bestWindow.loc.iStart)); int realGrow = newWin.length() -
       * bestWindow.loc.nLength; OatesModelUSamp om = new OatesModelUSamp(vdata.get(bestWindow.loc.iSeries),
       * newWin, params.initv, params.minv);
       * 
       * TimerMS tms = new TimerMS(); om.trainEM(vdata); // TODO: need to build vtrain
       * System.err.println("train time: " + tms.time() + " ms");
       * 
       * PVMJobSpool spool = new PVMJobSpool(jpvm, tids); NBestOMaps nBestMaps = new
       * NBestOMaps(params.nMinMatches); for(int i = 0; i < vdata.size(); i++){ fbm_req req = new fbm_req(i,
       * om); spool.add(req, MSG_PVM_FBM_REQ); } while(!spool.empty()){ msg = jpvm.pvm_recv(MSG_PVM_FBM_RET);
       * spool.pop(getWorkerIndex(msg.sourceTid)); fbm_ret fbm = new fbm_ret(msg.buffer); OatesMapping omap =
       * fbm.getMapping(); if (bVerbose) System.err.println("GF: Best local mapping: " + omap); if (omap !=
       * null) nBestMaps.add(omap); }
       * 
       * if (!acceptGrowth(om, vdata, nBestMaps, realGrow)) break; System.err.println(" growing forward");
       * dumpPattern(nBestMaps, vdata); bestWindow.loc.nLength = newWin.length(); bestWindow.maps = nBestMaps; } //
       * and now try to extend the start of the window while(bestWindow.loc.iStart > 0){ int iNewStart =
       * Math.max(bestWindow.loc.iStart - params.nwGrow, 0); int realGrow = bestWindow.loc.iStart - iNewStart;
       * int iNewStop = bestWindow.loc.end(); WindowLocation newWin = new WindowLocation(iNewStart, iNewStop -
       * iNewStart); OatesModelUSamp om = new OatesModelUSamp(vdata.get(bestWindow.loc.iSeries), newWin,
       * initv, minv); om.trainEM(vdata);
       * 
       * PVMJobSpool spool = new PVMJobSpool(jpvm, tids); NBestOMaps nBestMaps = new
       * NBestOMaps(params.nMinMatches); for(int i = 0; i < vdata.size(); i++){ fbm_req req = new fbm_req(i,
       * om); spool.add(req, MSG_PVM_FBM_REQ); } while(!spool.empty()){ msg = jpvm.pvm_recv(MSG_PVM_FBM_RET);
       * spool.pop(getWorkerIndex(msg.sourceTid)); fbm_ret fbm = new fbm_ret(msg.buffer); OatesMapping omap =
       * fbm.getMapping(); if (bVerbose) System.err.println("GB: Best local mapping: " + omap); if (omap !=
       * null) nBestMaps.add(omap); }
       * 
       * if (!acceptGrowth(om, vdata, nBestMaps, -realGrow)) break; System.err.println(" growing backward:");
       * dumpPattern(nBestMaps, vdata); bestWindow.loc.iStart = newWin.start(); bestWindow.loc.nLength =
       * newWin.length(); bestWindow.maps = nBestMaps; } System.err.println("\nfinal best window: " +
       * bestWindow); }
       */

      // break the data sequence up based on the min_matches best matches
      System.err.print("Extracting base patterns and chopping sequences... ");
      ArrayList<WindowLocation> vLocs = new ArrayList<WindowLocation>();
      for(int i = 0; i < bestWindow.maps.size(); i++)
         vLocs.add(bestWindow.maps.get(i).getLoc());

      // extract patterns and chop sequence
      ArrayList<Sequence> vPats = Sequence.extractPats(vdata, vLocs);
      vdata = Sequence.chop(vdata, vLocs, params.nwLen);
      System.err.println("done.");

      // now we can search for more occurrences of the pattern
      OatesMapping omap;
      while((params.nMaxInst <= 0 || vPats.size() < params.nMaxInst)
            && (omap = findNextPatternSpot(vPats, vdata, omGarbage, params.aSeries, params.nSamples)) != null){
         System.err.print("next pat: ");
         dumpPattern(omap, vdata.get(omap.iSeries));

         Sequence par = vdata.get(omap.iSeries);
         WindowLocation loc = omap.getLoc();

         // extract the pattern
         Sequence pat = par.subseq(loc);
         vPats.add(pat);

         // chop up the sequences
         Sequence.chopInPlace(vdata, loc, params.nwLen);
         System.err.printf("num pats: %d  max: %d\n", vPats.size(), params.nMaxInst);
      }

      // at this point, we have extracted every instance of the pattern and can output the info
      // output.print("****************************************");
      // output.println("***************************************");
      // output.println("* Pattern "+(iPattern+1));
      // output.print("****************************************");
      // output.println("***************************************");
      for(int i = 0; i < vPats.size(); i++){
         Sequence pat = vPats.get(i);
         // output.println((pat.getParentIndex()+1)+" "+(pat.getParentOffset()+1)+" "+pat.length());
         String alpha = "ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789";
         output.println((pat.getParentIndex() + 1) + "  " + "\"" + alpha.charAt(iPattern) + "\"  "
               + (pat.getParentOffset() + 1) + "  " + (pat.getParentOffset() + pat.length() + 1));
      }
      output.println();

      // save the model for this pattern
      /*
       * try{ // build a model from the discovered patterns OatesModelUSamp om = new OatesModelUSamp(vPats,
       * params.nwLen, initv, minv); String sFile = String.format("model_%02d.dat", iPattern+1);
       * //FileOutputStream fos = new FileOutputStream(sFile); //GZIPOutputStream gzos = new
       * GZIPOutputStream(fos); //ObjectOutputStream out = new ObjectOutputStream(gzos);
       * //out.writeObject(om); //out.flush();
       * 
       * PrintWriter out = new PrintWriter(new FileWriter(sFile)); out.println(om); out.close(); } catch
       * (IOException e) { System.err.println(e); }
       */

      // finally, we recurse to find the next pattern
      if (params.nMaxPats < 0 || iPattern + 1 < params.nMaxPats) discover(iPattern + 1, vdata, omGarbage);

      return true;
   }

   /** */
   public static boolean validateParams(ArrayList<Sequence> vdata)
   {
      if (params.nMinMatches > vdata.size()){
         System.err.printf("Error: requested more matches (nMinMatches = %d) " + "than time series (%d)\n",
               params.nMinMatches, vdata.size());
         return false;
      }
      return true;
   }

   public static void main(String args[])
   {
      params = new Params();

      int c, iw;
      LongOpt[] longopts = new LongOpt[] { new LongOpt("help", LongOpt.NO_ARGUMENT, null, 'h'),
            new LongOpt("nlen", LongOpt.REQUIRED_ARGUMENT, null, 'l'),
            new LongOpt("nskip", LongOpt.REQUIRED_ARGUMENT, null, 'k'),
            new LongOpt("ngrow", LongOpt.REQUIRED_ARGUMENT, null, 'g'),
            new LongOpt("path", LongOpt.REQUIRED_ARGUMENT, null, 'p'),
            new LongOpt("verbose", LongOpt.NO_ARGUMENT, null, 'v'),
            new LongOpt("nminmatch", LongOpt.REQUIRED_ARGUMENT, null, 'm'),
            new LongOpt("nmaxinst", LongOpt.REQUIRED_ARGUMENT, null, 'i'),
            new LongOpt("atime", LongOpt.REQUIRED_ARGUMENT, null, 't'),
            new LongOpt("aseries", LongOpt.REQUIRED_ARGUMENT, null, 's'),
            new LongOpt("nsamples", LongOpt.REQUIRED_ARGUMENT, null, 'n'),
            new LongOpt("out", LongOpt.REQUIRED_ARGUMENT, null, 'o'),
            new LongOpt("search", LongOpt.REQUIRED_ARGUMENT, null, 'e'),
            new LongOpt("nmarkers", LongOpt.REQUIRED_ARGUMENT, null, 'b'),
            new LongOpt("nmaxpats", LongOpt.REQUIRED_ARGUMENT, null, 'a'),
            new LongOpt("unipath", LongOpt.REQUIRED_ARGUMENT, null, 1003),
            new LongOpt("presplit", LongOpt.REQUIRED_ARGUMENT, null, 1004),
            new LongOpt("ngarbage", LongOpt.REQUIRED_ARGUMENT, null, 1005),
            new LongOpt("model", LongOpt.REQUIRED_ARGUMENT, null, 1006),
            new LongOpt("hmmeval", LongOpt.REQUIRED_ARGUMENT, null, 1007),
            new LongOpt("hmmstates", LongOpt.REQUIRED_ARGUMENT, null, 1008),
            new LongOpt("hmmskip", LongOpt.REQUIRED_ARGUMENT, null, 1009),
            new LongOpt("hmmtrain", LongOpt.REQUIRED_ARGUMENT, null, 1010),
            new LongOpt("gui", LongOpt.NO_ARGUMENT, null, 1011),
            new LongOpt("fixedvar", LongOpt.NO_ARGUMENT, null, 1012) };

      Getopt g = new Getopt("Discover Oates", args, "?", longopts, true);
      while((c = g.getopt()) != -1){
         switch(c){
         case '?':
         case 'h': // help
            usage();
            System.exit(0);
            break;
         case 'l': // window width
            params.nwLen = Integer.parseInt(g.getOptarg());
            if (params.given != null) params.given.nLength = params.nwLen;
            break;
         case 'k': // window skip
            params.nwSkip = Integer.parseInt(g.getOptarg());
            break;
         case 'g': // nGrow
            params.nwGrow = Integer.parseInt(g.getOptarg());
            break;
         case 'p': // path
            params.sPath = g.getOptarg();
            File fdir = new File(params.sPath);
            if (!fdir.exists() || !fdir.isDirectory()){
               System.err.printf("Error: invalid base directory\n (%s)\n", params.sPath);
               System.exit(1);
            }
            break;
         case 'v': // verbose
            bVerbose = true;
            break;
         case 'm': // nMinMatches
            params.nMinMatches = Integer.parseInt(g.getOptarg());
            break;
         case 'i': // nMaxInst
            params.nMaxInst = Integer.parseInt(g.getOptarg());
            break;
         case 'a': // nMaxPats
            params.nMaxPats = Integer.parseInt(g.getOptarg());
            break;
         case 't': // aTime
            params.aTemporal = Double.parseDouble(g.getOptarg());
            break;
         case 's': // aSeries
            params.aSeries = Double.parseDouble(g.getOptarg());
            break;
         case 'n': // nSamples
            params.nSamples = Integer.parseInt(g.getOptarg());
            break;
         case 'e': // search
         {
            String sType = g.getOptarg();
            if (sType.startsWith("given")){
               params.kSearch = Params.SEARCH_GIVEN;
               Pattern p = Pattern.compile("^given.??(\\d+)\\.(\\d+)$");
               Matcher m = p.matcher(sType);
               if (m.matches()){
                  int iSeq = Integer.parseInt(m.group(1)) - 1;
                  int iLoc = Integer.parseInt(m.group(2));
                  params.given = new WindowLocation(iSeq, iLoc, params.nwLen);
               }
            }
            else if (Library.stricmp(sType, "exhaustive")) params.kSearch = Params.SEARCH_EXHAUSTIVE;
            else if (Library.stricmp(sType, "exone")) params.kSearch = Params.SEARCH_EXHAUST_ONE;
            else if (Library.stricmp(sType, "uniform")) params.kSearch = Params.SEARCH_UNIFORM;
            else if (Library.stricmp(sType, "suftree")) params.kSearch = Params.SEARCH_SUFFIX_TREE;
            else{
               System.err.println("Error: unrecognized search style: " + sType);
               System.exit(1);
            }
         }
            break;
         case 'b': // nMarkers
            params.nMarkers = Integer.parseInt(g.getOptarg());
            break;
         case 'o': // out
            try{
               output = new PrintStream(new File(g.getOptarg()));
            } catch (FileNotFoundException fnfe){
               fnfe.printStackTrace();
               System.exit(1);
            }
            break;
         case 1004: // presplit
         {
            StringTokenizer st = new StringTokenizer(g.getOptarg(), ":");
            assert (st.countTokens() == 4) : "-presplit parameter should have form: w:t:g:o";
            params.splitW = Integer.parseInt(st.nextToken());
            params.splitT = Double.parseDouble(st.nextToken());
            params.splitG = Integer.parseInt(st.nextToken());
            params.splitO = Integer.parseInt(st.nextToken());
         }
            break;
         case 1005: // nGarbage
            params.nGarbage = Integer.parseInt(g.getOptarg());
            break;
         case 1006: // model
         {
            String sType = g.getOptarg();
            if (Library.stricmp(sType, "hmm")) params.model = Params.MODEL_HMM;
            else if (Library.stricmp(sType, "oates")) params.model = Params.MODEL_OATES;
            else{
               System.err.println("Error: unrecognized model type: " + sType);
               System.exit(1);
            }
         }
            break;
         case 1007: // hmmeval
            if (Library.stricmp(g.getOptarg(), "viterbi")) params.hmmEval = HMM_VITERBI;
            else if (Library.stricmp(g.getOptarg(), "forward")) params.hmmEval = HMM_FORWARD;
            else{
               System.err.println("Error: unrecognized \"hmmeval\" option: " + g.getOptarg());
               System.exit(1);
            }
            break;
         case 1008: // hmmstates
            params.nHmmStates = Integer.parseInt(g.getOptarg());
            break;
         case 1009: // hmmskip
            params.nHmmSkip = Integer.parseInt(g.getOptarg());
            break;
         case 1010: // hmmtrain
            if (Library.stricmp(g.getOptarg(), "viterbi")) params.hmmTrain = HMM_VITERBI;
            else if (Library.stricmp(g.getOptarg(), "bw")) params.hmmTrain = HMM_BW;
            else{
               System.err.println("Error: unrecognized \"hmmtrain\" option: " + g.getOptarg());
               System.exit(1);
            }
            break;
         case 1011: // gui
            dogui = new DOGui();
            break;
         case 1012: // fixedvar
            params.bUpdateVar = false;
         default:
            System.exit(1);
            break;
         }
      }

      // add file patterns to param object
      if (g.getOptind() >= args.length){
         System.err.println("Error: no file patterns found!");
         usage();
         System.exit(1);
      }
      for(int i = g.getOptind(); i < args.length; i++)
         params.addFilePat(args[i]);

      if (bVerbose){
         System.err.println("Parameters (v" + Params.VER + "):");
         System.err.println(" base path:                          " + params.sPath);
         System.err.println(" sliding window size:                " + params.nwLen);
         System.err.println(" sliding skip size:                  " + params.nwSkip);
         System.err.println(" temporal growth increment:          " + params.nwGrow);
         System.err.println(" given exemplar locatoin:            " + params.given);
         System.err.println(" minimum # occurrences:              " + params.nMinMatches);
         System.err.println(" a_temporal:                         " + params.aTemporal);
         System.err.println(" a_series:                           " + params.aSeries);
         System.err.println(" num samples for new occurence test: " + params.nSamples);
         System.err.println(" max number of patterns:             " + params.nMaxPats);
         System.err.println(" # occs for garbage model:           " + params.nGarbage);
         System.err.println(" knn:                                " + params.knn);
         System.err.println(" model:                              " + params.model);
         if (params.model == MODEL_HMM){
            System.err.println(" HMM states:                         " + params.nHmmStates);
            System.err.println(" HMM skip:                           " + params.nHmmSkip);
         }
         if (params.splitW > 0){
            System.err.println(" splitW:                             " + params.splitW);
            System.err.printf(" splitT:                             %.4f\n", params.splitT);
            System.err.println(" splitG:                             " + params.splitG);
            System.err.println(" splitO:                             " + params.splitO);
         }
         System.err.println(" num markers:                        " + params.nMarkers);
         System.err.println(" File Patterns:                      " + params.filePats.size());
         for(int i = 0; i < params.filePats.size(); i++)
            System.err.println("  " + params.filePats.get(i));
      }

      TimerMS timerTotal = new TimerMS();
      TimerMS timer = new TimerMS();

      // load the data
      DataLoader loader = new DLRaw();
      ArrayList<Sequence> vdata = new ArrayList<Sequence>();
      for(String filePat : params.filePats){
         if (filePat.endsWith(".def")){
            TreeMap<String, ArrayList<Sequence>> data = LabeledDataLoader.load(new File(params.sPath,
                  filePat));
            if (data == null){
               System.err.printf("Warning: failed to load labeled data\n (%s)\n", filePat);
               continue;
            }

            int N = LabeledDataLoader.tseries.size();
            for(int i = 0; i < N; i++){
               Sequence seq = LabeledDataLoader.tseries.get(i);
               seq.setParentIndex(vdata.size());
               vdata.add(seq);
            }
         }
         else{
            File[] files = Library.getFilesWild(params.sPath, filePat);
            if (files == null){
               System.err.printf("Warning: no files found from file pattern (%s)\n", filePat);
               continue;
            }
            for(File file : files){
               System.err.println("Loading: " + file.getPath());
               Sequence data = loader.load(file.getPath());
               data.setParentIndex(vdata.size());
               vdata.add(data);
            }
         }
      }
      if (vdata.size() < 1){
         System.err.printf("Error: no data sequences!\n");
         System.exit(1);
      }
      System.err.printf("Time to load data: %dms\n", timer.time());

      // compute global stats
      params.nDims = vdata.get(0).getNumDims();
      params.gmSeries = SupTest.calcGauss(vdata.toArray(new Sequence[0]));
      params.initv = SupTest.initv;
      params.minv = SupTest.minv;

      if (bVerbose){
         System.err.printf("Found %d data files:\n", vdata.size());
         for(int i = 0; i < params.nDims; i++)
            System.err.printf(
                  " Dim %d: mean=%.2f   var=%.2f (%.2f)   initv=%.2f (%.2f)   minv=%.2f (.2f)\n", i + 1,
                  params.gmSeries[i].getMean(), params.gmSeries[i].getVar(), params.gmSeries[i].getSDev(),
                  params.initv.get(i), Math.sqrt(params.initv.get(i)), params.minv.get(i), Math
                        .sqrt(params.minv.get(i)));
      }

      // start the discovery process
      if (validateParams(vdata)){
         // presplit the data if requested
         if (params.splitW > 0){
            presplit(vdata);
            if (bVerbose) System.err.printf("After presplit: %d sequences\n", vdata.size());
         }

         // build a garbage model
         OatesModelUSamp omGarbage = (params.nGarbage > 0 ? buildGarbageModel(vdata, params.nGarbage) : null);

         // create gui if requested
         if (dogui != null) dogui.setVisible(true);

         // now we can discover recurring patterns
         discover(0, vdata, omGarbage);
      }

      System.err.printf("Total run time: %dms\n", timerTotal.time());
   }
}
