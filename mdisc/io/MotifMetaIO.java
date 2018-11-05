package mdisc.io;

import java.io.*;
import kdm.data.*;
import java.util.*;
import mdisc.data.*;

/** interface for class that loads or saves motif meta data */
public interface MotifMetaIO
{
   public void saveMeta(PrintWriter out, Object meta);
   public Object loadMeta(BufferedReader in) throws IOException;   
   public Motif construct();
   public Motif construct(ArrayList<WindowLocation> _occs);
   public boolean hasOccCount();
   public int getOccCount();
}
