package mdisc.VizTool;

/** error caused when a collision matrix has too many or too few hits */ 
public class ColMatrixException extends Exception
{
   public final static int tooMany = 1;
   public final static int tooFew = 2;
   
   public final static String[] msg = new String[]{
     "Too many hits",
     "Too few hits"
   };
   
   public final int type;
   
   public ColMatrixException(int _type)
   {
      super(msg[_type]);
      type = _type;
   }
}
