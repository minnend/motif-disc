package mdisc.io;

import kdm.util.*;
import kdm.io.*;
import java.io.*;
import java.util.*;

public class PatternFileLoader
{
    public ArrayList<ArrayList<PatData>> pats;

    public PatternFileLoader()
    {
        pats = new ArrayList<ArrayList<PatData>>();
    }

    public boolean load(String sFile)
    {
        // each line of the pattern file should look like this:
        //  Occurrence 1: 8.1019  60
        // <occurrence index>, <series index>, <start position>, <length>
        
        pats.clear();
        try{
            CommentedFileReader reader = new CommentedFileReader(new FileReader(sFile), '*');
            String line;
            ArrayList<PatData> occs = new ArrayList<PatData>();

            // now we can read the pattern data
            while((line = reader.readLine()) != null)
            {
                StringTokenizer st = new StringTokenizer(line, " \t");
                String header = st.nextToken();
                if (Library.stricmp(header, "pattern"))
                {
                    if (!occs.isEmpty()) pats.add(occs);
                    occs = new ArrayList<PatData>();
                    continue;
                }
                int iSeries = Integer.parseInt(header);
                int pos = Integer.parseInt(st.nextToken());
                int len = Integer.parseInt(st.nextToken());
                PatData pd = new PatData(iSeries, pos, len);
                occs.add(pd);
            }
            if (!occs.isEmpty()) pats.add(occs);
            reader.close();
        }
        catch(Exception e){ e.printStackTrace(); return false; }
        return true;
    }

    public ArrayList<PatData> get(int i){ return pats.get(i); }
    public PatData get(int i, int j){ return pats.get(i).get(j); }

////////////////////////////////////////////////////////////

public class PatData
{
    public int iSeries, pos, len;
    
    public PatData(int _iSeries, int _pos, int _len)
    {
        iSeries = _iSeries;
        pos = _pos;
        len = _len;
    }

    public String toString()
    {
        return "[PatData: "+iSeries+"."+pos+"  |"+len+"]";
    }
}

}
