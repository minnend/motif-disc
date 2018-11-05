package mdisc.io;

import java.io.*;
import java.util.*;

import javax.swing.*;

import kdm.data.*;
import kdm.util.*;
import mdisc.data.*;

/** Provides file I/O for motifs */
public abstract class MotifIO
{
   /**
    * Save motifs in the given file
    * 
    * @param list of motifs motifs to save
    * @param f file in which to save motifs
    * @param mmio meta info io object
    * @return true if successful
    */
   public static boolean saveMotifs(ArrayList<Motif> motifs, File f, MotifMetaIO mmio)
   {
      try
      {
         PrintWriter out = new PrintWriter(new FileWriter(f));
         out.printf("%d\n", motifs.size());
         for(int i = 0; i < motifs.size(); i++)
         {
            Motif motif = motifs.get(i);
            if (mmio != null) mmio.saveMeta(out, motif.getMeta());
            if (mmio==null || !mmio.hasOccCount()) out.printf("%d\n", motif.getNumOccs());
            for(int j = 0; j < motif.getNumOccs(); j++)
               out.println(motif.getOcc(j).toText());
         }
         out.close();
         return true;
      } catch (IOException e)
      {
         return false;
      }
   }

   /**
    * Load motifs from the given file
    * 
    * @param f the file from which to load motifs
    * @return list of motifs (null on error)
    */
   public static ArrayList<Motif> loadSeeds(File f, MotifMetaIO mmio)
   {
      try
      {
         BufferedReader in = new BufferedReader(new FileReader(f));
         int nMotifs = Integer.parseInt(in.readLine());
         ArrayList<Motif> motifs = new ArrayList<Motif>();
         for(int i = 0; i < nMotifs; i++)
         {
            Motif motif;
            ArrayList<WindowLocation> occs = new ArrayList<WindowLocation>();
            
            if (mmio != null)
            {
               motif = mmio.construct();
               Object meta = mmio.loadMeta(in);
               if (meta != null) motif.setMeta(meta);
            }
            else{
               motif = new Motif();
            }
            int nOccs = (mmio != null && mmio.hasOccCount() ? mmio.getOccCount() : Integer.parseInt(in
                  .readLine()));
            
            for(int j = 0; j < nOccs; j++)
            {
               String line = in.readLine();
               WindowLocation wloc = WindowLocation.createFromText(line);
               if (wloc == null)
               {                  
                  System.err.println("Failed to parse wloc (" + line + ")");
                  throw new NullPointerException("Invalid wloc: " + line);
               }
               occs.add(wloc);
            }

            motif.setOccs(occs);
            motifs.add(motif);
         }
         in.close();
         return motifs;
      } catch (Exception e)
      {
         System.err.println("Error: failed to load seed motif locations in file:");
         System.err.println(" " + Library.getCanonical(f));
         return null;
      }
   }

}
