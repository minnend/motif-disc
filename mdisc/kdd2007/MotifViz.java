package mdisc.kdd2007;

import java.awt.*;
import kdm.data.*;
import kdm.gui.*;
import kdm.models.*;

/** visualize a motif HMM */
public class MotifViz extends JMyComponent
{
   protected AbstractHMM hmm;
   protected FeatureVec vMin, vMax;

   public MotifViz(AbstractHMM hmm, FeatureVec vMin, FeatureVec vMax)
   {
      setOpaque(true);
      this.hmm = hmm;
      this.vMin = vMin;
      this.vMax = vMax;
   }

   public void setHMM(AbstractHMM hmm)
   {
      this.hmm = hmm;
      repaint();
   }

   @Override
   public Dimension getPreferredSize()
   {
      return new Dimension(120, 120);
   }

   @Override
   public void paintComponent(Graphics2D g, int cw, int ch)
   {
      if (hmm == null){
         g.setColor(Color.white);
         g.fillRect(0, 0, cw, ch);
         return;
      }

      int nStates = hmm.getNumStates();
      int nDims = hmm.getNumDims();

      FeatureVec vRange = vMax.sub(vMin);
      for(int iState = 0; iState < nStates; iState++){
         int x = (int)Math.round((double)iState * cw / nStates);
         int x2 = (int)Math.round((double)(iState + 1) * cw / nStates);
         int w = x2 - x;
         ProbFVModel pm = hmm.getState(iState);
         for(int d = 0; d < nDims; d++){
            int y = (int)Math.round((double)d * ch / nDims);
            int y2 = (int)Math.round((double)(d + 1) * ch / nDims);
            int h = y2 - y;

            if (pm instanceof GaussianDiagonal){
               FeatureVec v = ((GaussianDiagonal)pm).getMean();

               float f = (float)((v.get(d) - vMin.get(d)) / vRange.get(d));
               if (f > 1) f = 1;
               else if (f < 0) f = 0;
               g.setColor(new Color(f, f, f));
               g.fillRect(x, y, w, h);
            }
            else if (pm instanceof GMM){
               GMM gmm = (GMM)pm;
               int nMix = gmm.getNumMix();
               for(int iMix=0; iMix<nMix; iMix++){
                  FeatureVec v = ((GaussianDiagonal)gmm.getComp(iMix)).getMean();
                  int xofs = (int)Math.round((double)iMix*w/nMix);
                  float f = (float)((v.get(d) - vMin.get(d)) / vRange.get(d));
                  g.setColor(new Color(f, f, f));                  
                  g.fillRect(x+xofs, y, Math.round(w-xofs), h);
               }
            }
            else{
               System.err.printf("Error: unsupported HMM observation distribution (%s)\n", pm.getClass()
                     .getName());
               g.setColor(Color.white);
               g.fillRect(x, y, w, h);
            }
         }
      }

   }

}
