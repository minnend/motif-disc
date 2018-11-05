package mdisc.kdd2007;

import java.awt.*;
import java.awt.image.*;
import java.io.*;
import java.util.*;
import javax.imageio.*;
import gnu.getopt.*;
import kdm.data.*;
import kdm.io.DataSaver.*;
import kdm.util.*;
import kdm.models.*;

/** generate synthetic character data */
public class GenCharData
{
   public enum Fill { white, uniform, randwalk }
   
   protected static char[] alpha = "abcdefghijklmnopqrstuvwxyz".toCharArray();
   
   protected static void usage()
   {
      System.err.println();
      System.err.println("USAGE: GenCharData [Options]");
      System.err.println("Options:");
      System.err.println(" -fontsize <size>    specify font size (32)");
      System.err.println(" -nchars <N>         number of characters per line (20)");
      System.err.println(" -nlines <N>         number of lines (10)");
      System.err.println(" -nfeats <N>         number of features per scan line (-1 = one per pixel)");
      System.err.println(" -shot <percent>     percent shot noise (none)");
      System.err.println(" -gauss <sdev>       std dev for Gaussian noise [0..1 scale] (none)");
      System.err.println(" -title <title>      title of output files = title_X.[png,txt,etc]");
      System.err.println(" -hscale <f>         divider for horizontal resolution (none)");
      System.err.println(" -hgap <w>           # additional pixels between each character (0)");
      System.err.println(" -hfill <fill>       method for filling gaps (white,uniform,randwalk)");
      System.err.println(" -alpha <chars>      specify which characters to use (abc..xyz)");
      System.err.println(" -font <type>        font style (serif,sans,mono)");
      System.err.println();
   }
   
   protected static int alphaIndex(char c)
   {
      for(int i=0; i<alpha.length; i++)
         if (alpha[i]==c) return i;
      return -1;
   }

   public static void main(String[] args) throws Exception
   {
      boolean bVerbose = false;
      String sTitle = null;
      int fontSize = 32;
      int nChars = 20;
      int nLines = 10;
      int nFeats = -1;
      double shot = 0;
      double gauss = 0;
      int hscale = 1;
      int hgap = 0;
      Fill fill = Fill.white;      
      String fontName = Font.SERIF;

      int c;
      LongOpt[] longopts = new LongOpt[] { new LongOpt("help", LongOpt.NO_ARGUMENT, null, 'h'),
            new LongOpt("v", LongOpt.NO_ARGUMENT, null, 'v'),
            new LongOpt("fontsize", LongOpt.REQUIRED_ARGUMENT, null, 1001),
            new LongOpt("nchars", LongOpt.REQUIRED_ARGUMENT, null, 1002),
            new LongOpt("nlines", LongOpt.REQUIRED_ARGUMENT, null, 1003),
            new LongOpt("nfeats", LongOpt.REQUIRED_ARGUMENT, null, 1004),
            new LongOpt("title", LongOpt.REQUIRED_ARGUMENT, null, 1005),
            new LongOpt("shot", LongOpt.REQUIRED_ARGUMENT, null, 1006),
            new LongOpt("gauss", LongOpt.REQUIRED_ARGUMENT, null, 1007),
            new LongOpt("hscale", LongOpt.REQUIRED_ARGUMENT, null, 1008),
            new LongOpt("hgap", LongOpt.REQUIRED_ARGUMENT, null, 1009),
            new LongOpt("fill", LongOpt.REQUIRED_ARGUMENT, null, 1010),
            new LongOpt("alpha", LongOpt.REQUIRED_ARGUMENT, null, 1011),
            new LongOpt("font", LongOpt.REQUIRED_ARGUMENT, null, 1012)
      };

      Getopt getopt = new Getopt("SupTest", args, "?", longopts, true);
      while((c = getopt.getopt()) != -1){
         String sArg = getopt.getOptarg();
         switch(c){
         case '?':
         case 'h': // help
            usage();
            System.exit(0);
            break;
         case 'v': // verbose
            bVerbose = true;
            break;
         case 1001: // fontsize
            fontSize = Integer.parseInt(sArg);
            break;
         case 1002: // nchars
            nChars = Integer.parseInt(sArg);
            break;
         case 1003: // nlines
            nLines = Integer.parseInt(sArg);
            break;
         case 1004: // nfeats
            nFeats = Integer.parseInt(sArg);
            break;
         case 1005: // title
            sTitle = sArg;
            break;
         case 1006: // shot
            shot = Double.parseDouble(sArg);
            break;
         case 1007: // gauss
            gauss = Double.parseDouble(sArg);
            break;
         case 1008: // hscale
            hscale = Integer.parseInt(sArg);
            break;
         case 1009: // hgap
            hgap = Integer.parseInt(sArg);
            break;
         case 1010: // fill
            if (sArg.equals("white")) fill = Fill.white;
            else if (sArg.equals("uniform")) fill = Fill.uniform;
            else if (sArg.equals("randwalk")) fill = Fill.randwalk;
            else{
               System.err.printf("\nError: unrecognized fill method (\"%s\")\n\n", sArg);
               System.exit(1);
            }
            break;
         case 1011: // alpha
            alpha = sArg.toCharArray();
            break;
         case 1012: // font
            if (sArg.equals("serif")) fontName = Font.SERIF;
            else if (sArg.equals("sans")) fontName = Font.SANS_SERIF;
            else if (sArg.equals("mono")) fontName = Font.MONOSPACED;
            else{
               System.err.printf("\nError: unrecognized font style (\"%s\")\n\n", sArg);
               System.exit(1);
            }
            break;
         }
      }

      if (sTitle == null){
         System.err.println();
         System.err.println("Error: must specify a title for data and images (use -title).");
         System.err.println();
         System.exit(1);
      }

      Font font = new Font(fontName, Font.PLAIN, fontSize);
      Graphics2D g = (Graphics2D)new BufferedImage(1, 1, BufferedImage.TYPE_3BYTE_BGR).getGraphics();
      FontMetrics fm = g.getFontMetrics(font);
      int h = fm.getHeight();
      if (nFeats <= 0) nFeats = h;

      if (bVerbose){
         System.err.println();
         System.err.println("Settings:");
         System.err.printf(" Font size: %d\n", fontSize);
         System.err.printf("  -> height: %d\n", h);
         System.err.printf(" Font style: %s\n", fontName);
         System.err.printf(" # characters per line: %d\n", nChars);
         System.err.printf(" # lines: %d\n", nLines);
         System.err.printf(" # features: %d\n", nFeats);
         System.err.printf(" Gaussian noise: %f\n", gauss);
         System.err.printf(" Shot noise: %.1f%%\n", shot);
         System.err.printf(" File title: %s\n", sTitle);
         System.err.printf(" Horizontal scale: %d\n", hscale);
         System.err.printf(" Horizontal gap: %d\n", hgap);
         System.err.printf(" Fill method: %s\n", fill);         
         System.err.printf(" Alphabet: %s\n", new String(alpha));
         System.err.println();
      }

      PrintWriter outInfo = new PrintWriter(new BufferedWriter(new FileWriter(new File(sTitle + ".info"))));
      PrintWriter outDef = new PrintWriter(new BufferedWriter(new FileWriter(new File(sTitle + ".def"))));
      outDef.printf("# synthetic character data: #chars=%d  fontSize=%d  gauss=%.2f  shot=%.2f\n\n",
            nChars, fontSize, gauss, shot);
      
      // generate the data
      DSRaw saver = new DSRaw();
      saver.setFormat("%.2f ");
      int nAlpha = alpha.length;
      int[] count = new int[nAlpha];
      int[] charWidth = new int[nAlpha];
      for(int i = 0; i < nAlpha; i++)
         charWidth[i] = fm.stringWidth("" + alpha[i]);

      for(int iLine = 0; iLine < nLines; iLine++){
         // generate info for this line
         char[] chars = new char[nChars];
         int[] xchar = new int[nChars];
         int w = (nChars-1)*hgap;
         for(int i = 0; i < nChars; i++){
            int iSym = Library.random(nAlpha);
            chars[i] = alpha[iSym];
            count[iSym]++;
            w += charWidth[iSym];
         }

         // generate the image of this line
         BufferedImage img = new BufferedImage(w, h, BufferedImage.TYPE_BYTE_GRAY);
         g = (Graphics2D)img.getGraphics();
         g.setFont(font);
         Library.setAntiAlias(g, false);
         g.setColor(Color.white);
         g.fillRect(0, 0, w, h);

         g.setColor(Color.black);
         Library.setAntiAlias(g, true);
         g.setRenderingHint(RenderingHints.KEY_TEXT_ANTIALIASING, RenderingHints.VALUE_TEXT_ANTIALIAS_ON);
         // g.drawString(new String(chars), 0, h - fm.getDescent() - 2);
         int xc = 0;
         for(int i = 0; i < nChars; i++){
            String s = "" + chars[i];
            g.drawString(s, xc, h - fm.getDescent() - 2);
            xchar[i] = xc;
            xc += charWidth[alphaIndex(chars[i])];
            if (hgap > 0){
               xc += hgap;
               if (fill == Fill.uniform){
                  // TODO
               }
               else if (fill == Fill.randwalk){
                  // TODO
               }
            }            
         }
         g.dispose();

         // corrupt with shot noise?
         if (shot > 0){
            float[] a = new float[1];
            WritableRaster wr = img.getRaster();
            int nPixels = w * h;
            int nShot = (int)Math.round(shot * nPixels / 100.0);
            int[] shotPixels = Library.selectRandomIndices(nShot, nPixels);
            for(int i = 0; i < nShot; i++){
               int x = shotPixels[i] % w;
               int y = shotPixels[i] / w;
               a[0] = (float)Library.random(256);
               wr.setPixel(x, y, a);
            }
         }

         // corrupt with Gaussian noise?
         if (gauss > 0){
            float[] a = new float[1];
            Gaussian1D g1 = new Gaussian1D(0, gauss * gauss);
            WritableRaster wr = img.getRaster();
            for(int y = 0; y < h; y++)
               for(int x = 0; x < w; x++){
                  double dc = g1.sample1();
                  wr.getPixel(x, y, a);
                  a[0] = (float)Math.max(0, Math.min(255, a[0] + dc));
                  wr.setPixel(x, y, a);
               }
         }
         
         // scale image?
         if (hscale > 1){
            BufferedImage img2 = new BufferedImage(w/hscale, h, img.getType());
            g = (Graphics2D)img2.getGraphics();
            g.setRenderingHint(RenderingHints.KEY_INTERPOLATION, RenderingHints.VALUE_INTERPOLATION_BICUBIC);
            g.drawImage(img, 0, 0, w/hscale, h, null);
            g.dispose();
            img = img2;
            w = img.getWidth();
            h = img.getHeight();
         }

         // save the image
         String sFile = String.format("%s_%04d.png", sTitle, iLine + 1);
         File file = new File(sFile);
         if (bVerbose) System.err.printf("Saving image: %s\n", sFile);
         try{
            ImageIO.write(img, "png", file);
         } catch (IOException e){
            System.err.printf("Error: failed to save image: %s\n", file.getAbsolutePath());
            System.exit(1);
         }

         // save the data
         sFile = String.format("%s_%04d.txt", sTitle, iLine + 1);
         if (bVerbose) System.err.printf("Saving data: %s\n", sFile);
         Sequence seq = new Sequence(sFile);
         Raster raster = img.getData();
         float[] col = new float[h];
         float[] feats = (nFeats < h ? new float[nFeats] : col);
         for(int i = 0; i < w; i++){
            FeatureVec fv = new FeatureVec(nFeats);
            raster.getPixels(i, 0, 1, h, col);
            if (nFeats < h){
               for(int d=0; d<nFeats; d++){
                  int a = (int)Math.round((double)d * h / nFeats);                  
                  int b = (int)Math.max(Math.round((double)(d+1) * h / nFeats), a+1);                  
                  int n = b-a;
                  float sum = 0;
                  while(a<b) sum += col[a++];
                  feats[d] = sum / n;
               }               
            }            
            for(int d=0; d<nFeats; d++) fv.set(d, feats[d]);
            seq.add(fv);            
         }
         if (!saver.save(seq, sFile)){
            System.err.printf("Error: failed to save data: %s\n", file.getAbsolutePath());
            System.exit(1);
         }

         // save meta data
         outInfo.printf("%s  %d\n", sFile, nChars);
         for(int i = 0; i < nChars; i++)
            outInfo.printf("%c ", chars[i]);
         outInfo.println();
         for(int i = 0; i < nChars; i++)
            outInfo.printf("%d ", xchar[i]/hscale);
         outInfo.println();
         for(int i = 0; i < nChars; i++)
            outInfo.printf("%d ", charWidth[alphaIndex(chars[i])]/hscale);
         outInfo.println();
         
         // save seq in data def file
         outDef.println("<data>");
         outDef.printf("name = %s_%04d\n", sTitle, iLine+1);
         outDef.printf("data = %s\n", sFile);
         outDef.println("dataLoader = DLRaw");
         outDef.printf("labels = %s_%04d.labels\n", sTitle, iLine+1);
         outDef.println("labelLoader = MLGeneral");
         outDef.println("</data>");
         outDef.println();
         
         // save label file
         sFile = String.format("%s_%04d.labels", sTitle, iLine+1);
         if (bVerbose) System.err.printf("Saving labels: %s\n", sFile);
         PrintWriter out = new PrintWriter(new BufferedWriter(new FileWriter(new File(sFile))));
         for(int i=0; i<nChars; i++)
            out.printf("\"%s\" %d %d\n", chars[i], xchar[i]/hscale, (xchar[i]+charWidth[alphaIndex(chars[i])])/hscale);
         out.close();
      }
      outInfo.close();
      outDef.close();
      if (bVerbose){
         System.err.print("Counts: ");
         for(int i=0; i<nAlpha; i++) System.err.printf("%d ", count[i]);
         System.err.println();
         System.err.print("Widths: ");
         for(int i=0; i<nAlpha; i++) System.err.printf("%d ", charWidth[i]/hscale);
         System.err.println();
      }
   }
}
