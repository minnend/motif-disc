package mdisc.aaai2007;

import kdm.data.*; 
import kdm.util.*;

/** stores information about a nearest neighbor */
public class NNInfo implements Comparable, Conflicter
{
   /** index of the nearest neighbor */
   public int index;
   
   /** distance to the nearest neighbor */
   public double dist;
   
   /** window location of this subsequence */
   public WindowLocation wloc;
   
   public NNInfo(int index, double dist, WindowLocation wloc)
   {
      this.index = index;
      this.dist = dist;
      this.wloc = wloc;
   }
   
   public int compareTo(Object o)
   {
      NNInfo nni = (NNInfo)o;
      if (dist < nni.dist) return -1;
      if (dist > nni.dist) return 1;
      return 0;
   }
   
   @Override
   public boolean equals(Object o)
   {
      return (dist == ((NNInfo)o).dist);
   }

   public boolean hasConflict(Object o)
   {
      return wloc.overlaps(((NNInfo)o).wloc);
   }

   @Override
   public String toString()
   {
      return String.format("[%d|%.2f %s]", index, dist, wloc);
   }
}
