package mdisc.data;

import java.util.*;
import kdm.data.*;

/** Encapsulates a motif (a set of similar subsequences) */
public class Motif
{
   protected ArrayList<WindowLocation> occs;
   protected Object meta;
   
   public Motif()
   {
      this(null);
   }
   public Motif(ArrayList<WindowLocation> occs)
   {
      this.occs = occs;
   }     
   
   public void setOccs(ArrayList<WindowLocation> _occs){ occs = _occs; }
   public ArrayList<WindowLocation> getOccs(){ return occs; }
   public int getNumOccs(){ return occs.size(); }
   public WindowLocation getOcc(int i){ return occs.get(i); }
   public void setMeta(Object _meta){ meta = _meta; }
   public Object getMeta(){ return meta; }   
}
