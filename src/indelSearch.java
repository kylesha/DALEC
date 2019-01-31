public class indelSearch
{
    public static void main(String[] args)
    {
        Oligo keyA = new Oligo                                   ("GTTCAAGGTCGTTCCTCnATaTC");
        Oligo A01 = new Oligo("GTTACGAATTgCCCTTNCCTcCTCGCGCCAtCAGAGTcTCAAGGTCGTTCNTCgATCTCTCCGACTCTACGTTTTCAGAATCTGGATCACCAGCTGCTCATAGACTGTACCAGTCTGAGCGGGCTGGCAAGGCAAGGGCGAATTCGCGGCCGCTAAATTCAATTCGCCCTATAGTGAGTCGTATTACAATTCACTGGCCGGTCGTTTTACAA");
        
        System.out.println(indelSearch(A01, keyA, 2));
    } //end main()
    
    
//////////////////////// indelSearch() accessory methods \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
//---------------- returns the source index that gives the best match ---------------//
    private static int indelSearch(Oligo s, Oligo t, int am)
    {
        String oligo = s.toString().toLowerCase();
        String key = t.toString().toLowerCase();
        char ignoredChar = 'n';
        final int OLIGO_LENGTH = oligo.length();
        final int KEY_LENGTH = key.length();
        float[] scoreKeeper = new float[KEY_LENGTH];
        float[] maxScoreArray = new float[KEY_LENGTH]; //keeps score array of best match
        int allowedMismatches = am;
        int nnMismatches = 0;
        int ntMatches = 0; //ntMatches = number of nucleotide-nucleotide (GC or AT) matches
        int ignoredCharMatches = 0; //ignoredCharMatches = number of nucleotide-ignoredChar (i.e. 'N') matches
        int bestMatchIndex = 0;
        float total = 0;
        float maxScore = 0;
        
        for(int i = 0; i + KEY_LENGTH <= OLIGO_LENGTH; i++) //find bestMatchIndex
        {
            for(int j = 0; j <= KEY_LENGTH - 1; j++)
            {
                boolean nnMatch = oligo.charAt(i+j) == key.charAt(j);
                boolean niMatch = oligo.charAt(i+j) == ignoredChar || key.charAt(j) == ignoredChar;
                
                if(nnMatch && !niMatch) //count only nucleotide-nucleotide matches
                {
                    scoreKeeper[j] = 1;
                    ntMatches++;
                }
                else if(niMatch)  //count only nucleotide-ignoredChar matches
                {
                    scoreKeeper[j] = 0;
                    ignoredCharMatches++;
                }
                else if(!nnMatch && !niMatch) //count only nucleotide-nucleotide mismatches
                {
                    scoreKeeper[j] = -1;
                    nnMismatches++;
                }
            } //inner for() loop
            
            //tally up scores for current sliding window, rewarding each consecutive match with an extra 0.5 points and penalizing non-consecutive matches
            for(int c = 0; c <= scoreKeeper.length - 2; c++)
            {
                if(scoreKeeper[c] >= 0 && (scoreKeeper[c+1] == 1 || scoreKeeper[c+1] == 0))
                    total += scoreKeeper[c] + (float)0.5;
                else total += scoreKeeper[c] - (float)0.5;
            }

            if(total > maxScore)
            {
                maxScore = total;
                bestMatchIndex = i;
                maxScoreArray = scoreKeeper.clone();
            }
            else //reset all counters
            {
                ntMatches = 0;
                ignoredCharMatches = 0;
                nnMismatches = 0;
                total = 0;
            }
        }
        
        //obtain restricted searchWindow with 2nt buffer on each side
        String searchWindow = new String();
        boolean withinRangeLeft = bestMatchIndex - 2 >= 0;
        boolean withingRangeRight = bestMatchIndex + KEY_LENGTH + 2 <= OLIGO_LENGTH - 1;
        
        if(withinRangeLeft && withingRangeRight)
        {
            int leftBoundary = bestMatchIndex - 2;
            int rightBoundary = bestMatchIndex + KEY_LENGTH + 2;
            searchWindow = new String(oligo.substring(leftBoundary, rightBoundary));
        }
        else
        {
            int leftBoundary = bestMatchIndex;
            int rightBoundary = bestMatchIndex + KEY_LENGTH;
            searchWindow = new String(oligo.substring(leftBoundary, rightBoundary));
        }
        
        System.out.println("bestMatchIndex: " + bestMatchIndex);
        System.out.println("\nsearchWindow: " + searchWindow);
        
        for(float elements: maxScoreArray)
            System.out.print(elements + " ");
        
        return demarcatorIndex(maxScoreArray);
    } //end indelSearch() method
    
    
    //------------------------ demarcatorIndex() -----------------------------//
    private static int demarcatorIndex(float[] array)
    {
        float[] maxScoreArray = array;
        int demarcatorIndex = 0;
        float tempLeft = 0;
        float tempRight = 0;
        float maxSum = 0;
        float larger = 0;
        
        for(int i = 0; i <= maxScoreArray.length - 1; i++)
        {
            for(int j = i; j >= 0; j--) //get sum left of i in maxScoreIndex
            {
                if(j > 0)
                {
                    if(maxScoreArray[j] >= 0 && (maxScoreArray[j-1] == 1 || maxScoreArray[j-1] == 0))
                        tempLeft += maxScoreArray[j] + (float)0.5;
                    else tempLeft += maxScoreArray[j] - (float)0.5;
                }
                else tempLeft += maxScoreArray[j];
            }
            
            for(int j = i; j <= maxScoreArray.length - 1; j++)  //get sum right of i in maxScoreIndex
            {
                if(j < maxScoreArray.length - 1)
                {
                    if(maxScoreArray[j] >= 0 && (maxScoreArray[j+1] == 1 || maxScoreArray[j+1] == 0))
                        tempRight += maxScoreArray[j] + (float)0.5;
                    else tempRight += maxScoreArray[j] - (float)0.5;
                }
                else tempRight += maxScoreArray[j];
            }
            
            System.out.println("\ni: " + i);
            System.out.println("tempLeft: " + tempLeft);
            System.out.println("tempRight: " + tempRight);
            
            if(tempLeft != tempRight)
            {
                larger = (tempLeft > tempRight) ? tempLeft : tempRight;
                if(larger > maxSum)
                {
                    maxSum = larger;
                    demarcatorIndex = i;
                    System.out.println("larger: " + larger);
                    System.out.println("maxSum: " + maxSum);
                    System.out.println("demarcatorIndex: " + demarcatorIndex);
                }
            }
            tempLeft = 0;
            tempRight = 0;
        } //end for() loop
        return demarcatorIndex;
    } //end demarcatorIndex()
} //end indelSearch class
