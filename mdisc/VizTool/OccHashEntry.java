package mdisc.VizTool;

import kdm.data.*;
import java.util.*;

/**
 * Holds info that serves as the entry for (fixed) motif occurrences during hashing
 */
public class OccHashEntry
{
   public int nOccsFixed, nTrans, nUniq;
   public double mdl;
   public String sMotif;
   public ArrayList<WindowLocation> occs;
   
   public OccHashEntry()
   {
      this(null,0,0,0);
   }
   
   public OccHashEntry(String _sMotif, int _nOccsFixed, int _nUniq, int _nTrans)
   {
      sMotif = _sMotif;
      nOccsFixed = _nOccsFixed;
      nUniq = _nUniq;
      nTrans = _nTrans;      
      mdl = Double.NaN;
   }
   
   public int getNumOccsDTW(){ return (occs==null ? nTrans : occs.size()); }
   
   public int getNumFramesUsed()
   {
      if (occs == null) return 0;
      int n = 0;
      for(WindowLocation wloc : occs) n += wloc.nLength;
      return n;
   }
   
   public static Comparator getComparator(int iSortBy)
   {
      if (iSortBy == HashView.SORTBY_NUM_OCCS)
      {
         return new Comparator()
         {
            public int compare(Object o1, Object o2)
            {
               OccHashEntry occA = (OccHashEntry)o1;
               OccHashEntry occB = (OccHashEntry)o2;
               if (occA.getNumOccsDTW() > occB.getNumOccsDTW()) return -1;
               if (occA.getNumOccsDTW() < occB.getNumOccsDTW()) return 1;
               if (occA.nOccsFixed > occB.nOccsFixed) return -1;
               if (occA.nOccsFixed < occB.nOccsFixed) return 1;
               if (occA.nTrans > occB.nTrans) return -1;
               if (occA.nTrans < occB.nTrans) return 1;
               if (occA.nUniq > occB.nUniq) return -1;
               if (occA.nUniq < occB.nUniq) return 1;
               if (!Double.isNaN(occA.mdl) && !Double.isNaN(occB.mdl))
               {
                  if (occA.mdl < occB.mdl) return -1;
                  if (occA.mdl > occB.mdl) return 1;
               }
               return occA.sMotif.compareTo(occB.sMotif);
            }
            
            public boolean equals(Object o1, Object o2)
            {
               return compare(o1, o2) == 0;
            }
         };
      }
      else if (iSortBy == HashView.SORTBY_NUM_TRANS)
      {
         return new Comparator()
         {
            public int compare(Object o1, Object o2)
            {
               OccHashEntry occA = (OccHashEntry)o1;
               OccHashEntry occB = (OccHashEntry)o2;
               if (occA.nTrans > occB.nTrans) return -1;
               if (occA.nTrans < occB.nTrans) return 1;
               if (occA.nUniq > occB.nUniq) return -1;
               if (occA.nUniq < occB.nUniq) return 1;
               if (occA.getNumOccsDTW() > occB.getNumOccsDTW()) return -1;
               if (occA.getNumOccsDTW() < occB.getNumOccsDTW()) return 1;
               if (occA.nOccsFixed > occB.nOccsFixed) return -1;
               if (occA.nOccsFixed < occB.nOccsFixed) return 1;
               if (!Double.isNaN(occA.mdl) && !Double.isNaN(occB.mdl))
               {
                  if (occA.mdl < occB.mdl) return -1;
                  if (occA.mdl > occB.mdl) return 1;
               }
               return occA.sMotif.compareTo(occB.sMotif);
            }
            
            public boolean equals(Object o1, Object o2)
            {
               return compare(o1, o2) == 0;
            }
         };
      }
      else if (iSortBy == HashView.SORTBY_MDL || iSortBy == HashView.SORTBY_MDL_RLE)
      {
         return new Comparator()
         {
            public int compare(Object o1, Object o2)
            {
               OccHashEntry occA = (OccHashEntry)o1;
               OccHashEntry occB = (OccHashEntry)o2;
               if (!Double.isNaN(occA.mdl) && !Double.isNaN(occB.mdl))
               {
                  if (occA.mdl < occB.mdl) return -1;
                  if (occA.mdl > occB.mdl) return 1;
               }
               if (occA.nTrans > occB.nTrans) return -1;
               if (occA.nTrans < occB.nTrans) return 1;
               if (occA.nUniq > occB.nUniq) return -1;
               if (occA.nUniq < occB.nUniq) return 1;
               if (occA.getNumOccsDTW() > occB.getNumOccsDTW()) return -1;
               if (occA.getNumOccsDTW() < occB.getNumOccsDTW()) return 1;
               if (occA.nOccsFixed > occB.nOccsFixed) return -1;
               if (occA.nOccsFixed < occB.nOccsFixed) return 1;
               return occA.sMotif.compareTo(occB.sMotif);
            }
            
            public boolean equals(Object o1, Object o2)
            {
               return compare(o1, o2) == 0;
            }
         };
      }
      else if (iSortBy == HashView.SORTBY_STRING)
      {
         return new Comparator()
         {
            public int compare(Object o1, Object o2)
            {               
               OccHashEntry occA = (OccHashEntry)o1;
               OccHashEntry occB = (OccHashEntry)o2;
               int ret = occA.sMotif.compareTo(occB.sMotif);
               if (ret!=0) return ret;
               if (occA.nOccsFixed > occB.nOccsFixed) return -1;
               if (occA.nOccsFixed < occB.nOccsFixed) return 1;
               if (occA.getNumOccsDTW() > occB.getNumOccsDTW()) return -1;
               if (occA.getNumOccsDTW() < occB.getNumOccsDTW()) return 1;
               if (!Double.isNaN(occA.mdl) && !Double.isNaN(occB.mdl))
               {
                  if (occA.mdl < occB.mdl) return -1;
                  if (occA.mdl > occB.mdl) return 1;
               }
               if (occA.nTrans > occB.nTrans) return -1;
               if (occA.nTrans < occB.nTrans) return 1;
               if (occA.nUniq > occB.nUniq) return -1;
               if (occA.nUniq < occB.nUniq) return 1;
               return 0;
            }
            
            public boolean equals(Object o1, Object o2)
            {
               return compare(o1, o2) == 0;
            }
         };
      }
      assert false : "Invalid SORTBY value: "+iSortBy;
      return null;
   }
}
