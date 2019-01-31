import java.io.*;
import java.util.*;
import java.util.ArrayList;

public class parseInSilicoClonesToUnique
{
    public static void main(String args[]) throws IOException
    {
        Scanner fin = new Scanner(new File("~/DALEC/Documents and Settings/firelab/Desktop/test_in/test.txt"));
        
        HashMap<Oligo, String[]> element = new HashMap<Oligo, String[]>();
        DpnIFrag frag = new DpnIFrag();
        
        int count = 0;
        ArrayList<Integer> processed = new ArrayList<Integer>(20000);
        
        while(fin.hasNext())
        {
            frag = DpnIFrag.toDpnIObj(fin.next());
            ArrayList<DpnIFrag> temp = new ArrayList<DpnIFrag>();
            String[] coordinates = new String[temp.size()];
            Oligo oligo = frag.getSeq();
            int key = frag.getKey();
            
            if(!processed.contains((Integer)key))
            {
                processed.add(key);
                
            }
        }
    }//end main()
}
