import java.io.*;
import java.util.*;

public class objectRAFv1
{
    public static void main(String[] args) throws IOException
    {
        //open inputstream and define DpnIFrag object
        Scanner fin = new Scanner(new File("~/DALEC/test/test_input.txt"));
        DpnIFrag dpn = new DpnIFrag();
        
        //obtain byte array for temp storage of objects
        ByteArrayOutputStream bos = new ByteArrayOutputStream();
        ObjectOutput oos = new ObjectOutputStream(bos);
        
        ArrayList objectList = new ArrayList<int[]>();
        
        RandomAccessFile raf = new RandomAccessFile("~/DALEC/test/testRAF.bin", "rw");
        
        int bufLength;
        int position = 0;
        
        //serialize objects to file, one at a time
        dpn = DpnIFrag.toDpnIObj(fin.nextLine());
        oos.writeObject(dpn);
            System.out.println("size of bos: " + bos.size());
        
        //convert byte array to bytes
        byte[] buf = bos.toByteArray();
        int[] objectInfo = new int[2];
        
        bufLength = buf.length;
        objectInfo[0] = position;    //get the identifier of present object
        objectInfo[1] = bufLength;      //get the length of present object
        objectList.add(objectInfo);
            System.out.println("position: " + objectInfo[0] + " bufLength: " + objectInfo[1]);
            System.out.println(dpn.toString());
        
        //write to file
        raf.seek(position);
        raf.write(buf);
        bos.reset();
        fin.close();
        raf.close();
        System.out.println("finished serializing");
    }//end main()
}

