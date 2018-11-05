package mdisc.aaai2007;

import java.io.*;
import java.nio.*;
import java.nio.channels.*;

/** tests an i/o issue in java 1.6 beta */
public class Test
{
   public static float getVal(int i)
   {
      return (float)(i + Math.PI*Math.pow(2.41,i+.92));
   }
   
   public static void gen(ByteBuffer bb, int N)
   {
      bb.putInt(1234);
      bb.putInt(4321);
      bb.putShort((short)1234);
      bb.putShort((short)4321);
      for(int i=0; i<N; i++){
         bb.putFloat(getVal(i));
      }
   }
   
   public static void save(ByteBuffer bb, File file) throws Exception
   {
      FileChannel fc = new FileOutputStream(file).getChannel();
      bb.flip();
      fc.write(bb);
      fc.close();
   }

   public static void test(String name, File file, int N) throws Exception
   {
      System.err.printf("File size (%s): %d\n", name, file.length());
      DataInputStream in = new DataInputStream(new FileInputStream(file));
      if (in.readInt() != 1234){
         System.err.printf("Test \"%s\" failed: first int failed\n", name);
         System.exit(1);
      }
      if (in.readInt() != 4321){
         System.err.printf("Test \"%s\" failed: secont int failed\n", name);
         System.exit(1);
      }
      if (in.readShort() != 1234){
         System.err.printf("Test \"%s\" failed: first short failed\n", name);
         System.exit(1);
      }
      if (in.readShort() != 4321){
         System.err.printf("Test \"%s\" failed: second short failed\n", name);
         System.exit(1);
      }
      for(int i=0; i<N; i++){
         float f = in.readFloat();
         float x = getVal(i);
         if (f != x){
            System.err.printf("Test \"%s\" failed: expected %f, found %f!\n", name, x, f);
            System.exit(1);
         }
      }
      in.close();
   }

   public static void main(String[] args) throws Exception
   {
      File file = File.createTempFile("test", "bin");

      int N = 10000;

      ByteBuffer bb = ByteBuffer.allocate(N * 4 + 1024);
      bb.order(ByteOrder.BIG_ENDIAN);
      gen(bb, N);

      save(bb, file);
      test("Normal", file, N);

      bb = ByteBuffer.allocateDirect(N * 4 + 1024);
      bb.order(ByteOrder.BIG_ENDIAN);
      gen(bb, N);
      
      save(bb, file);
      test("Direct", file, N);
   }
}
