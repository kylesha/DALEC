public class SmithWaterman
{
    public static void main(String[] args)
    {
        String source = "xGGGTTACGAATTGCCCTTNCCTCCTCGCGCCATCAGAGTTCAAGGTCGTACCTCATTCTCTCCGACTCTACGTTTTCAGAATCTGGATCACCAGCTGCTCATAGACTGTACCAGTCTGAGCGGGCTGGCAAGGCAAGGGCGAATTCGCGGCCGCTAAATTCAATTCGCCCTATAGTGAGTCGTATTACAATTCACTGGCCGGTCGTTTTACAA";
        String query  = "xGTTCAAGTCGTTCCTCATTCTC";
    }//end main()
    
    public static Oligo SmithWaterman(Oligo source, Oligo query)
    {
        char[] s = source.toCharArray();
        char[] q = query.toCharArray();
        float[][] matrix = new float[q.length + 1][s.length + 1];
        
        //initialize row 0 to 0.0
        for(int i = 0; i <= q.length; i++)
            matrix[i][0] = (float)0.0;
        
        //initialize column 0 to 0.0
        for(int j = 0; j <= s.length; j++)
            matrix[0][j] = (float)0.0;
        
        float max = 0;
        float best = 0; //highest score in each comparison
        float iScore = 0;
        float hScore = 0;
        float diagScore = 0; //nucleotide-nucleotide alignment score: match = 1.0; mismatch = -0.3
        int imax = 0;
        int jmax = 0;
        
        //construct scores matrix
        for(int i = 1; i <= q.length - 1; i++)
        {
            for(int j = 1; j <= s.length - 1; j++)
            {
                //calculate diagonal score
                diagScore = (q[i] == s[j]) ? (float)1.0 : (float)-0.3;
                best = matrix[i-1][j-1] + diagScore; //initially, asume diagonal is best score
                
                //calc max vertical score
                for(int vGap = i; vGap >= 0; vGap--)
                {
                    iScore = matrix[i-vGap][j] - (float)(1 + 0.3*vGap);
                    best = (iScore > best) ? iScore : best;
                }
                
                //calc max horizontal score
                for(int hGap = j; hGap >= 0; hGap--)
                {
                    hScore = matrix[i][j-hGap] - (float)(1 + 0.3*hGap);
                    best = (hScore > best) ? hScore : best;
                }
                
                matrix[i][j] = (best > 0) ? best : 0;
                if(matrix[i][j] > max)
                {
                    max = matrix[i][j];
                    imax = i;
                    jmax = j;
                }
            }//for(int j...) loop
        } //for(int i...) loop
        
        for(int i = 0; i <= q.length - 1; i++)
        {
            for(int j = 0; j <= s.length - 1; j++)
                System.out.printf("%2.1f ", matrix[i][j]);
            System.out.println();
        }
        System.out.println("max: " + max);
        System.out.println("imax: " + imax);
        System.out.println("jmax: " + jmax);
        
        //reconstruct alignment
        float above = 0;
        float left = 0;
        float diag = 0;
        boolean isAbove = false; //best score is cell above
        boolean isDiag = false;  //best score is diagonal cell
        boolean isLeft = false;  //best score is cell to left
        StringBuilder key = new StringBuilder();
        key = key.append(q[imax]);  //initialize key to ch associated with max score
        
        do
        {
//            System.out.println("imax = " + imax + " jmax = " + jmax + "  matrix[imax][jmax] = " + matrix[imax][jmax]);
            
            above = matrix[imax-1][jmax];
            diag = matrix[imax-1][jmax-1];
            left = matrix[imax][jmax-1];
            
            isAbove = (above >= diag) && (above >= left);
            isDiag = (diag >= above) && (diag >= left);
            isLeft = (left >= diag) && (left >= above);
            
            if(isAbove && !isDiag) //deletion in source sequence
                imax--;
            
            else if(isDiag) //nucleotide-nucleotide match
            {
                imax--;
                jmax--;
                
                if(imax > 0) //prevent q[0] from being appended to key
                    key.append(q[imax]);
            }
            
            else if(isLeft && !isDiag) //insertion in source sequence; insert 'n' into query to compensate
            {
                jmax--;
                key.append('n');
            }
            
        } while(matrix[imax][jmax] > 0);
        
        return new Oligo(key.reverse());
    }
} //end class SmithWaterman
