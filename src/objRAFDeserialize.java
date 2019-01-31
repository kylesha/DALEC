import java.io.*;

public class objRAFDeserialize
{
    public static void main(String[] args) throws IOException, ClassNotFoundException
    {
        RandomAccessFile raf = new RandomAccessFile("~/DALEC/test/testRAF2.bin", "r");
        byte[] buf = new byte[(int)raf.length()];
        raf.seek(0);
        raf.readFully(buf);
        
        ByteArrayInputStream bis = new ByteArrayInputStream(buf);
        ObjectInputStream ois = new ObjectInputStream(bis);
        DpnIFrag dpn = new DpnIFrag();
        dpn = (DpnIFrag)ois.readObject();
        System.out.println(dpn.toString());
        ois.close();
        raf.close();
    }//end main()
}
