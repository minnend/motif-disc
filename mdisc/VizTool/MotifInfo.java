package mdisc.VizTool;

import java.awt.*;
import java.io.*;
import java.util.*;

import kdm.data.*;
import mdisc.io.*;
import mdisc.data.*;

/**
 * Stores information about a motif.
 */
public class MotifInfo extends Motif implements MotifMetaIO
{
   protected String sHash;
   protected int nCountFixed, nCountDTW;
   protected double mdl;   
   
   public MotifInfo()
   {
      sHash = "ERROR";
      nCountFixed = nCountDTW = -1;
      mdl = 0;
   }
   
   public MotifInfo(String _sHash, int _nCountFixed, int _nCountDTW, double _mdl, ArrayList<WindowLocation> _occs)
   {
      super(_occs);
      sHash = _sHash;
      nCountFixed = _nCountFixed;      
      nCountDTW = _nCountDTW;
      mdl = _mdl;      
   }   
   
   public int getCountFixed(){ return nCountFixed; }
   public int getCountDTW(){ return (occs==null ? 0 : occs.size()); }
   public double getMDL(){ return mdl; }
   
   public String getHash(){ return sHash; }
   public int getLength(){ return sHash.length(); }
      
   public String toString(){ return String.format("%s (%d,%d) %.2f", sHash, nCountFixed, getCountDTW(), mdl); }

   public void saveMeta(PrintWriter out, Object meta)
   {
      MotifInfo mi = (MotifInfo)meta;
      out.printf("%s %d %d %f\n", mi.sHash, mi.nCountFixed, mi.nCountDTW, mi.mdl);
   }

   public Object loadMeta(BufferedReader in) throws IOException
   {
      StringTokenizer st = new StringTokenizer(in.readLine());
      sHash = st.nextToken();
      nCountFixed = Integer.parseInt(st.nextToken());
      nCountDTW = Integer.parseInt(st.nextToken());
      mdl = Double.parseDouble(st.nextToken());
      return this;
   }
   
   public void setMeta(Object _meta)
   {
      if (_meta != this)
      {
         MotifInfo mi = (MotifInfo)_meta;
         sHash = mi.sHash;
         nCountFixed = mi.nCountFixed;
         nCountDTW = mi.nCountDTW;
         mdl = mi.mdl;
      }
   }
   
   public Object getMeta()
   {
      return this;
   }

   public Motif construct()
   {
      return new MotifInfo();
   }

   public Motif construct(ArrayList<WindowLocation> _occs)
   {
      MotifInfo mi = new MotifInfo();
      mi.setOccs(_occs);
      return mi;
   }

   public boolean hasOccCount()
   {
      return true;
   }

   public int getOccCount()
   {
      return nCountDTW;
   }
}
