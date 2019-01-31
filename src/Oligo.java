/**
 *@author Ky Sha
 *A class dedicated to manipulating/searching Oligo objects
 *VERSIONING: xx.yy.zz [xx = major revision/addition/deletion; yy = addition of new method(s); zz = minor modifications (to class/methods/etc.)]
 *@version 2.0 (07/30/2007) - Addition of exception handling to class.
 *@version 2.1 (08/05/2007) - addition of isFuzzySmithWatermanMatch()
 *-------------------------------------------------------------------------------------------------
 *METHODS ON REQUEST
 * replace()
 * mutateINS() - insertional mutation
 * mutateDEL() - deletion mutation
 * mutateTranspose() - transpositional insertion
 *label5Prime() - adds adduct (i.e. biotin) to 5 prime end
 *label3Prime() - adds adduct to 3 primer end
 *explicitMode() - show 5'/3' labels of specified oligo
 *constructors strip away punctuation marks (including space) embedded in input oligo
 *-------------------------------------------------------------------------------------------------
 */

import java.util.*;
import java.lang.*;
import java.io.*;
import java.util.regex.*;

public class Oligo implements Serializable
{
    private String oligo;
    private String mnf = "MatchNotFound"; //exception
    private String soob = "StartIndexOutOfBounds"; //exception
    private String eoob = "EndIndexOutOfBounds"; //exception
    private char ignoredChar = 'n'; //default character to ignore in oligo
    private final int OLIGO_LENGTH;
    
    
//============================| CONSTRUCTORS |================================//
    /**
     *Constructor: Constructs an oligo object of arbitrary length
     *@param oligo (type: String, StringBuilder, StringBuffer, void)
     *@return Oligo object
     */
    public Oligo(String oligo)
    {
        this.oligo = oligo.toLowerCase();
        OLIGO_LENGTH = oligo.length();
    }
    
    
    
    public Oligo(StringBuilder oligo)
    {
        this.oligo = oligo.toString().toLowerCase();
        OLIGO_LENGTH = oligo.length();
    }
    
    
    
    public Oligo(StringBuffer oligo)
    {
        this.oligo = oligo.toString().toLowerCase();
        OLIGO_LENGTH = oligo.length();
    }
    
    
    
    public Oligo(Oligo oligo)
    {
        this.oligo = oligo.toString().toLowerCase();
        OLIGO_LENGTH = oligo.length();
    }
    
    
    
    public Oligo()
    {
        oligo = "";
        OLIGO_LENGTH = oligo.length();
    }
    
    
    
//=========================| OVER RIDDEN METHODS |============================//
    /**
     *Converts oligo to character array
     *@param void
     *@return CharArray representation of oligo
     */
    public char[] toCharArray()
    {
        return oligo.toCharArray();
    }
    
    
    
    /**
     *Converts all characters in oligo to lowercase characters
     *@return Oligo object - oligo in lowercase
     */
    public Oligo toLowerCase()
    {
        return new Oligo(oligo.toLowerCase());
    }
    
    
    
    /**
     *Converts oligo to String object
     *@param void
     *@return String object representation of oligo
     */
    public String toString()
    {
        return oligo;
    }
    
    
    
    /**
     *Converts oligo to StringBuilder object
     *@param void
     *@return StringBuilder object representation of oligo
     */
    public StringBuilder toStringBuilder()
    {
        return new StringBuilder(oligo);
    }
    
    
    
    
//=================================| METHODS |================================//
    /**
     *Returns the reverse complement of the input oligo
     */
    public Oligo antiparallel()
    {
        Oligo tempOligo = new Oligo(oligo);
        Oligo complement = new Oligo(tempOligo.complement());
        return complement.reverse();
    }
    
    
    
    //--------------------------returns the base compositions of oligo (as a String()---------------------------------//
    public String baseContent()
    {
        int guanine = 0;
        int adenine = 0;
        int thymine = 0;
        int cytosine = 0;
        char[] source = oligo.toCharArray();
        
        for(int i = 0; i <= OLIGO_LENGTH - 1; i++)
        {
            if(source[i] == 'G' || source[i] == 'g')
                guanine++;
            else if(source[i] == 'A' || source[i] == 'a')
                adenine++;
            else if(source[i] == 'T' || source[i] == 't')
                thymine++;
            else if(source[i] == 'C' || source[i] == 'c')
                cytosine++;
        }
        return new String("G=" + guanine + " A=" + adenine + " T=" + thymine + " C=" + cytosine);
    }
    
    
    
    /**
     *Performs search using binary algorithm - attemps to incorporate indels and/or mismatches
     *@param query The query sequence to search for
     *@return Oligo object representing best match; otherwise, throws MatchNotFoundException
     *@throws OligoException
     */
    public Oligo binarySearch(Oligo query) throws OligoException
    {
        int nBothNotFound = 0; //counts number of times both keys not found
        int nRightKeyFound = 0; //counts number of times right key is not found
        int nLeftKeyFound = 0; //counts number of times right key is not found
        int rightBoundary = 0;
        int leftBoundary = 0;
        int keySize = 0;
        
        Oligo searchWindow = new Oligo(oligo);
        Oligo searchKey = query;
        Oligo[] splitKey = new Oligo[2];
        Oligo leftKey = new Oligo();
        Oligo rightKey = new Oligo();
        Oligo remainingKey = new Oligo();
        
        ArrayList tempKey = new ArrayList();    //stores temporary, unassembled key
        ArrayList matchList = new ArrayList();
        
        boolean foundLeftKey;
        boolean foundRightKey;
        boolean foundRightKeyIsDownstream;
        boolean willNeverFind;
        
        do
        {
            splitKey = halve(searchKey);
            leftKey = splitKey[0];
            rightKey = splitKey[1];
            
            //conditions that will never lead to a match
            if(willNeverFind(searchKey, searchWindow) || willNeverFind(searchWindow, searchKey))
            {
                searchKey = searchKey.resetToNull();
                searchWindow = searchWindow.resetToNull();
            }
            else
            {
                matchList = searchWindow.getAllMatchCoordinates(rightKey, 0);
                
                //define search criteria
                if(matchList.size() > 1)
                    foundRightKeyIsDownstream = searchWindow.lastIndex(leftKey, 0) < searchWindow.getLastMatchCoordinate(rightKey, 0); //rightKey match must be downstream of leftKey match
                else foundRightKeyIsDownstream = searchWindow.lastIndex(leftKey, 0) < searchWindow.getFirstMatchCoordinate(rightKey, 0);
                
                foundLeftKey = searchWindow.isFuzzyMatch(leftKey, 0);
                foundRightKey = searchWindow.isFuzzyMatch(rightKey, 0) && foundRightKeyIsDownstream;
                
                if(foundLeftKey ^ foundRightKey)
                {
                    if(foundLeftKey)
                    {
                        ++nLeftKeyFound;
                        tempKey.add("\\" + leftKey);
                        
                        //if searchWindow is less than either leftKey or rightKey, then maximum possible match is <50% ==> condition for a failed search
                        if(searchWindow.length() >= leftKey.length() || searchWindow.length() >= rightKey.length())
                        {
                            leftBoundary = searchWindow.lastIndex(leftKey, 0) + 1;
                            
                            //rightKey hasn't been found ==> set rightBoundary to extend 2nt beyond lenth of rightKey
                            if(nRightKeyFound == 0)
                            {
                                int lastIndexSearchWindow = searchWindow.length() - 1;
                                if(lastIndexSearchWindow - searchWindow.lastIndex(leftKey, 0) >= rightKey.length() + 2)
                                    rightBoundary = searchWindow.lastIndex(leftKey, 0) + rightKey.length() + 2;
                                else if(lastIndexSearchWindow - searchWindow.lastIndex(leftKey, 0) == rightKey.length() + 1)
                                    rightBoundary = searchWindow.lastIndex(leftKey, 0) + rightKey.length() + 1;
                                else if(lastIndexSearchWindow - searchWindow.lastIndex(leftKey, 0) <= rightKey.length())
                                    rightBoundary = searchWindow.length() - 1;
                            }
                            //rightKey has been found at least once; rightBoundary is simply the last index
                            else if(nRightKeyFound >= 1)
                                rightBoundary = searchWindow.length() - 1;
                            
                            searchWindow = searchWindow.extractSequence(leftBoundary, rightBoundary);
                            searchKey = rightKey;   //rightKey is now the new searchKey
                            
//                            System.out.println("FOUND LEFT KEY");
//                            System.out.println("LEFT...leftKey: " + leftKey);
//                            System.out.println("LEFT...rightKey: " + rightKey);
//                            System.out.println("LEFT...leftBoundary: " + leftBoundary);
//                            System.out.println("LEFT...rightBoundary: " + rightBoundary);
//                            System.out.println("LEFT...searchWindow: " + searchWindow);
//                            System.out.println("LEFT...searchKey: " + searchKey);
//                            System.out.println("nRightKeyFound: " + nRightKeyFound);
//                            System.out.println();
                        }
                        else searchKey = searchKey.resetToNull();
                    } //end LEFT found
                    
                    else if(foundRightKey)
                    {
                        ++nRightKeyFound;
                        
                        //store found key
                        if(searchWindow.length() == 1)
                        {
                            tempKey.add("/" + rightKey);
                            searchWindow = searchWindow.resetToNull(); //searchWindow is null ==> end search
                        }
                        
                        else if(searchWindow.length() >= leftKey.length() || searchWindow.length() >= rightKey.length())
                        {
                            tempKey.add("/" + rightKey);
                            
                            if(searchWindow.getFirstMatchCoordinate(rightKey, 0) >= 1)
                                rightBoundary = searchWindow.getFirstMatchCoordinate(rightKey, 0) - 1;
                            else /*if(searchWindow.getFirstMatchCoordinate(rightKey, 0) < 1)*/
                                rightBoundary = searchWindow.getFirstMatchCoordinate(rightKey, 0);
                            
                            if(nLeftKeyFound == 0) //leftKey hasn't been found ==> set leftBoundary to be 2nt outside of itself
                            {
                                if(rightBoundary - leftKey.length() > 2)
                                    leftBoundary = rightBoundary - leftKey.length() - 3;
                                else if(rightBoundary - leftKey.length() >= 2)
                                    leftBoundary = rightBoundary - leftKey.length() - 2;
                                else if(rightBoundary - leftKey.length() == 1)
                                    leftBoundary = rightBoundary - leftKey.length() - 1;
                                else if(rightBoundary - leftKey.length() < 1)
                                    leftBoundary = 0;
                            }
                            else if(nLeftKeyFound >= 1)
                                leftBoundary = 0;
                            
                            searchWindow = searchWindow.extractSequence(leftBoundary, rightBoundary);
                            searchKey = leftKey;    //define new searchKey
                        }
//                        System.out.println("FOUND RIGHT KEY");
//                        System.out.println("RIGHT...leftKey: " + leftKey);
//                        System.out.println("RIGHT...rightKey: " + rightKey);
//                        System.out.println("RIGHT...leftBoundary: " + leftBoundary);
//                        System.out.println("RIGHT...rightBoundary: " + rightBoundary);
//                        System.out.println("RIGHT...searchWindow: " + searchWindow);
//                        System.out.println("RIGHT...searchKey: " + searchKey);
//                        System.out.println("nLeftKeyFound: " + nLeftKeyFound);
//                        System.out.println();
                    } //end RIGHT found
                } //end ONE OR THE OTHER
                
                else if(foundLeftKey && foundRightKey)
                {
                    if(searchWindow.getAllMatchCoordinates(rightKey, 0).size() > 1)
                    {
                        //reset left and right borders and set new searchWindow
                        leftBoundary = searchWindow.getFirstMatchCoordinate(leftKey, 0);
                        rightBoundary = searchWindow.getLastMatchCoordinate(rightKey, 0) + rightKey.length() - 1;
                        searchWindow = searchWindow.extractSequence(leftBoundary, rightBoundary);
                        
                        if(nBothNotFound == 0)
                            tempKey.add("|" + searchWindow);  //store as a middle key
                        else if(nBothNotFound == 1)
                            tempKey.add("\\" + searchWindow); //store as a left key
                        else if(nBothNotFound == 2)
                            tempKey.add("/" + searchWindow);  //store as a right key
                        
                        searchWindow = searchWindow.resetToNull();
                    }
                    else
                    {
                        //reset left and right borders and set new searchWindow
                        leftBoundary = searchWindow.getFirstMatchCoordinate(leftKey, 0);
                        rightBoundary = searchWindow.lastIndex(rightKey, 0);
                        searchWindow = searchWindow.extractSequence(leftBoundary, rightBoundary);
                        
                        if(nBothNotFound == 0)
                            tempKey.add("|" + searchWindow);  //store as a middle key
                        else if(nBothNotFound == 1)
                            tempKey.add("\\" + searchWindow); //store as a left key
                        else if(nBothNotFound == 2)
                            tempKey.add("/" + searchWindow);  //store as a right key
                        
                        searchWindow = searchWindow.resetToNull();
                    }
                    
                    if(!remainingKey.isEmpty())
                    {
                        searchKey = remainingKey;  //transfer remainingKey to current searchKey
                        searchWindow = new Oligo(oligo);  //new searchWindow is entire the initial entire source sequence
                        remainingKey = remainingKey.resetToNull(); //reset remainingKey to null
                    }
                    else searchKey = searchKey.resetToNull();   //DONE...end of search!!
                    
//                    System.out.println("FOUND BOTH KEYS");
//                    System.out.println("BOTH...leftBoundary: " + leftBoundary);
//                    System.out.println("BOTH...rightBoundary: " + rightBoundary);
//                    System.out.println("BOTH...leftKey: " + leftKey);
//                    System.out.println("BOTH...rightKey: " + rightKey);
//                    System.out.println("BOTH...searchWindow AFTER: " + searchWindow);
//                    System.out.println("BOTH...tempKey Array: " + tempKey);
//                    System.out.println();
                }  //end BOTH found
                
                else if(!foundLeftKey && !foundRightKey)
                {
                    ++nBothNotFound;
                    
                    //1st time in through this loop: use leftKey as searchKey
                    if(nBothNotFound == 1)
                    {
                        searchKey = leftKey;
                        remainingKey = rightKey;
                        nLeftKeyFound = 0; //reset left most boundary for new search
                    }
                    //2nd time in through loop: use rightKey as searchKey
                    else if(nBothNotFound == 2)
                    {
                        splitKey = halve(query);
                        searchKey = splitKey[1];
                        nRightKeyFound = 0; //reset right most boundary for new search
                    }
                    //3rd time through loop: search failed
                    else if(nBothNotFound == 3)
                        return new Oligo();
                    
//                    System.out.println("BOTH NOT FOUND");
//                    System.out.println("BOTH NOT FOUND...searchWindow: " + searchWindow);
//                    System.out.println("NOT BOTH...searchKey: " + searchKey);
//                    System.out.println("NOT BOTH...leftKey: " + leftKey);
//                    System.out.println("NOT BOTH...rightKey: " + rightKey);
//                    System.out.println();
                } //end BOTH NOT found
            } //if willNeverFind()
        } while(!searchKey.isEmpty() && !searchWindow.isEmpty());
        
        //search fails if less than half of key is found
        for(int i = 0; i <= tempKey.size() - 1; i++)
            keySize += tempKey.get(i).toString().length();
        
        if(keySize - tempKey.size() < query.length()/2) //tempKey.size() => # key frags and hence, equals number of slashes (\, |, /)
            throw new OligoException(mnf, "binarySearch()");
        else return assembleKey(tempKey);
        
    }  //end binarySearch() method
    
    
    
    /**
     *binarySearch() accessory method: assembles the found key for binarySearch()
     *@param inputKey key to be assembled and returned to binarySearch()
     *@@return ArrayList object
     */
    private static Oligo assembleKey(ArrayList inputKey)
    {
        ArrayList key = inputKey;
        StringBuilder leftKey = new StringBuilder();
        StringBuilder rightKey = new StringBuilder();
        StringBuilder middleKey = new StringBuilder();
        
        while(!key.isEmpty())
        {
            if(key.get(0).toString().startsWith("\\"))
            {
                leftKey.append(key.get(0).toString().substring(1));
                key.remove(0);
            }
            else if(key.get(0).toString().startsWith("/"))
            {
                rightKey.insert(0, key.get(0).toString().substring(1));
                key.remove(0);
            }
            else if(key.get(0).toString().startsWith("|"))
            {
                middleKey.append(key.get(0).toString().substring(1));
                key.remove(0);
            }
        }
        return new Oligo(leftKey.append(middleKey.append(rightKey)));
    } //end assembleKey()
    
    
    
    /**
     *binarySearch() accessory method: divides query into two halves
     *@param query sequence to be halved
     *@@return Oligo object
     */
    private static Oligo[] halve(Oligo query)
    {
        Oligo target = query;
        Oligo[] output = new Oligo[2];
        output[0] = new Oligo(target.toString().substring(0, target.toString().length()/2));    //obtain left half of string
        output[1] = new Oligo(target.toString().substring(target.toString().length()/2));   //obtain right half of string
        return output;
    }
    
    
    
    /**
     *binarySearch() accessory method
     */
    private static boolean willNeverFind(Oligo target, Oligo source)
    {
        String t = "[" + target.toString() + "]";
        String s = "[" + source.toString() + "]";
        
        //if source does not contain a base in the target
        Pattern pat = Pattern.compile(t);
        Matcher mat = pat.matcher(s);
        boolean cond1 = mat.find();
        
        //vice versa
        pat = Pattern.compile(s);
        mat = pat.matcher(t);
        boolean cond2 = mat.find();
        
        boolean cond3 = target == source.reverse() && !target.equals(source);
        
        if(cond1 || cond2 || cond3)
            return false;
        else return true;
    } //end willNeverFind()
    
    
    
    /**
     *Returns the complement of the input oligo object
     *@return Oligo object
     *
     */
    public Oligo complement()
    {
        String tempOligo = "";
        char[] source = oligo.toCharArray();
        
        for(int i = 0; i <= OLIGO_LENGTH - 1; i++)
        {
            if(source[i] == 'a')
                tempOligo += 't';
            else if(source[i] == 'c')
                tempOligo += 'g';
            else if(source[i] == 'g')
                tempOligo += 'c';
            else if(source[i] == 't')
                tempOligo += 'a';
            else
                tempOligo += source[i];
        }
        return new Oligo(tempOligo);
    }
    
    
    public boolean contains(String s)
    {
        return (oligo.contains(s)) ? true: false;
    }
    
    
    public boolean contains(Oligo query)
    {
        return (oligo.contains(query.toString())) ? true: false;
    }
    
    
    
    /**
     *Counts the specified character in oligo object
     *@param inputChar character to be counter
     *@return int represent number of character counted
     */
    public int countChar(char inputChar)
    {
        char ch = inputChar;
        int count = 0;
        for(int i = 0; i <= OLIGO_LENGTH - 1; i++)
        {
            if(oligo.charAt(i) == ch || oligo.charAt(i) == ch)
                count++;
        }
        return count;
    }
    
    
    
    /**
     *Excises everthing left of query, including the query sequence
     *@param query search key
     *@param mismatches max number of allowed mismatches
     *@return [Oligo object] Resulting oligo sequence after excision; otherwise throws OligoException
     *@throws OligoException
     */
    public Oligo exciseLeftFrom(Oligo query, int mismatches) throws OligoException
    {
        if(isFuzzyMatch(query, mismatches))
        {
            int index = getFirstMatchCoordinate(query, mismatches);
            return new Oligo(oligo.substring(index + query.length()));
        }
        else throw new OligoException(mnf, "exciseLeftFrom()"); // exhausted all possibilities, no matches found
    } // end exciseLeftFrom()
    
    
    
    /**
     *Excises everthing left of query, including the query sequence
     *@param query query sequence
     *@param mismatches maximum number of allowedMismatches
     *@param ins maximum allowed number of inserts in SOURCE sequence
     *@param del maximum allowed number of deletes in SOURCE sequence
     *@param minKeyLength minimum length of found key in order for search to be considered successful
     *@return [Oligo object] Resulting oligo sequence after excision; otherwise throws OligoException
     *@throws OligoException
     */
    public Oligo exciseLeftFrom(Oligo query, int mismatches, int ins, int del, int minKeyLength) throws OligoException
    {
        if(isFuzzySmithWatermanMatch(query, mismatches, ins, del, minKeyLength))
        {
            Oligo swKey = getFuzzySWkey(query, mismatches, ins, del, minKeyLength);
            return new Oligo(exciseLeftFrom(swKey, mismatches));
        }
        else throw new OligoException(mnf, "exciseLeftFrom()"); // exhausted all possibilities, no matches found
    } // end exciseLeftFrom()
    
    
    
    /**
     *Excises everthing left of query
     *@param query search key
     *@param mismatches max number of allowed mismatches
     *@return [Oligo object] Resulting oligo sequence after excision; otherwise throws OligoException
     *@throws OligoException
     */
    public Oligo exciseLeftOf(Oligo query, int mismatches) throws OligoException
    {
        if(isFuzzyMatch(query, mismatches))
        {
            int index = getFirstMatchCoordinate(query, mismatches);
            return new Oligo(oligo.substring(index));
        }
        else throw new OligoException(mnf, "exciseLeftOf()"); // exhausted all possibilities, no matches found
    }//end exciseLeftOf()
    
    
    
    /**
     *Excises everthing left of query
     *@param query query sequence
     *@param mismatches maximum number of allowedMismatches
     *@param ins maximum allowed number of inserts in SOURCE sequence
     *@param del maximum allowed number of deletes in SOURCE sequence
     *@param minKeyLength minimum length of found key in order for search to be considered successful
     *@return [Oligo object] Resulting oligo sequence after excision; otherwise throws OligoException
     *@throws OligoException
     */
    public Oligo exciseLeftOf(Oligo query, int mismatches, int ins, int del, int minKeyLength) throws OligoException
    {
        int index = getFirstMatchCoordinate(query, mismatches, ins, del, minKeyLength);
        if(index != -1)
            return new Oligo(oligo.substring(index));
        else throw new OligoException(mnf, "exciseLeftOf()");
    }//end exciseLeftOf()
    
    
    
    /**
     *Excises everthing right of query, including the query sequence
     *@param query search key
     *@param mismatches max number of allowed mismatches
     *@return [Oligo object] Resulting oligo sequence after excision; otherwise throws OligoException
     *@throws OligoException
     */
    public Oligo exciseRightFrom(Oligo query, int mismatches) throws OligoException
    {
        if(isFuzzyMatch(query, mismatches))
        {
            int index = getFirstMatchCoordinate(query, mismatches);
            return new Oligo(oligo.substring(0, index));
        }
        else throw new OligoException(mnf, "exciseRightFrom()"); // exhausted all possibilities, no matches found
    } //end xciseRightFrom()
    
    
    
    /**
     *Excises everthing right of query, including the query sequence
     *@param query query sequence
     *@param mismatches maximum number of allowedMismatches
     *@param ins maximum allowed number of inserts in SOURCE sequence
     *@param del maximum allowed number of deletes in SOURCE sequence
     *@param minKeyLength minimum length of found key in order for search to be considered successful
     *@return [Oligo object] Resulting oligo sequence after excision; otherwise throws OligoException
     *@throws OligoException
     */
    public Oligo exciseRightFrom(Oligo query, int mismatches, int ins, int del, int minKeyLength) throws OligoException
    {
        int index = getFirstMatchCoordinate(query, mismatches, ins, del, minKeyLength);
        if(index != -1)
            return new Oligo(oligo.substring(0, index));
        else throw new OligoException(mnf, "exciseRightFrom()");
    }//end xciseRightFrom()
    
    
    
    /**
     *Excises everthing right of query
     *@param query search key
     *@param mismatches max number of allowed mismatches
     *@return [Oligo object] Resulting oligo sequence after excision; otherwise throws OligoException
     *@throws OligoException
     */
    public Oligo exciseRightOf(Oligo query, int mismatches) throws OligoException
    {
        if(isFuzzyMatch(query, mismatches))
        {
            int index = getFirstMatchCoordinate(query, mismatches);
            return new Oligo(oligo.substring(0, index + query.length()));
        }
        else throw new OligoException(mnf, "exciseRightOf()"); // exhausted all possibilities, no matches found
    } //end exciseRightOf()
    
    
    
    /**
     *Excises everthing right of query
     *@param query query sequence
     *@param mismatches maximum number of allowedMismatches
     *@param ins maximum allowed number of inserts in SOURCE sequence
     *@param del maximum allowed number of deletes in SOURCE sequence
     *@param minKeyLength minimum length of found key in order for search to be considered successful
     *@return [Oligo object] Resulting oligo sequence after excision; otherwise throws OligoException
     *@throws OligoException
     */
    public Oligo exciseRightOf(Oligo query, int mismatches, int ins, int del, int minKeyLength) throws OligoException
    {
        if(isFuzzySmithWatermanMatch(query, mismatches, ins, del, minKeyLength))
        {
            int index = getFirstMatchCoordinate(query, mismatches, ins, del, minKeyLength);
            Oligo swKey = getFuzzySWkey(query, mismatches, ins, del, minKeyLength);
            return new Oligo(oligo.substring(0, index + swKey.length()));
        }
        else throw new OligoException(mnf, "exciseRightOf()");
    }//end exciseRightOf()
    
    
    
    /**
     *Extracts the first best match sequence from oligo, based on fuzzy-match algorithm
     *@param query the sequence to extract
     *@param mismatches the maximum number of allowed mismatches
     *@return [Oligo object] Resulting oligo sequence after excision; otherwise throws OligoException
     *@throws MatchNotFoundException
     */
    public Oligo extractSequence(Oligo query, int mismatches) throws OligoException
    {
        if(isFuzzyMatch(query, mismatches))
        {
            if(getFirstMatchCoordinate(query, mismatches) != -1)
            {
                int start = getFirstMatchCoordinate(query, mismatches);
                int stop = start + query.length();
                return new Oligo(oligo.substring(start, stop));
            }
            else throw new OligoException(mnf, "extractSequence()");
        }
        else throw new OligoException(mnf, "extractSequence()");
    } //end extractSequence()
    
    
    
    public Oligo extractSequence(Oligo query, int mismatches, int ins, int del, int minKeyLength) throws OligoException
    {
        Oligo swKey = new Oligo();
        try
        {
            swKey = getFuzzySWkey(query, mismatches, ins, del, minKeyLength);
        }
        catch (OligoException error)
        {error.getStackTrace();}
        return extractSequence(swKey, mismatches);
    } //end extractSequence()
    
    
    
    /**
     *Extracts a sequence based on given start and end indices, inclusively
     *@param start start index
     *@param end end index
     *@return [Oligo object] Extracted source sequence; otherwise throws OligoException
     *@throws OligoException
     */
    public Oligo extractSequence(int start, int end) throws OligoException
    {
        String s = new String();
        boolean startIndexOutofBounds = start < 0;
        boolean endIndexOutofBounds = end > oligo.length() - 1;
        
        if(start >= 0 && end <= oligo.length() - 1)
            s = oligo.substring(start, end + 1);
        else if(startIndexOutofBounds ^ endIndexOutofBounds)
        {
            if(startIndexOutofBounds)
                throw new OligoException(soob, "extractSequence()");
            else throw new OligoException(eoob, "extractSequence()");
        }
        else if (startIndexOutofBounds && endIndexOutofBounds)
            throw new OligoException("(Start&End)IndexOutofBounds", "extractSequence()");
        return new Oligo(s);
    } //end extractSequence()
    
    
    
    /**
     *Generates a random base
     *@return char
     */
    public static char generateRandomBase()
    {
        char ch = ' ';
        int base = (int)(Math.random() * 100);
        
        if(base >= 0 && base <= 24)
            ch = 'a';
        else if(base >= 25 && base <= 49)
            ch = 'c';
        else if(base >= 50 && base <= 74)
            ch = 'g';
        else if(base >= 75 && base <= 99)
            ch = 't';
        return ch;
    } //end generateRandomBase()
    
    
    
    /**
     *Generate a random oligo of a specified length
     *@param length length of random oligo
     *@return Oligo object
     */
    public static Oligo generateRandomOligo(int length)
    {
        String randomOligo = new String();
        
        for(int i = 0; i <= length - 1; i++)
        {
            int base = (int)(Math.random() * 100);
            
            if(base >= 0 && base <= 24)
                randomOligo += 'a';
            else if(base >= 25 && base <= 49)
                randomOligo += 'c';
            else if(base >= 50 && base <= 74)
                randomOligo += 'g';
            else if(base >= 75 && base <= 99)
                randomOligo += 't';
        }
        return new Oligo(randomOligo);
    } //end generateRandomOligo()
    
    
    
    /**
     *Returns an ArrayList containing the indices of all occurences of the query
     *@param query query sequence
     *@param mismatches maximum number of allowedMismatches
     *@return ArrayList object (ArrayList is empty if zero matches are found)
     */
    public ArrayList getAllMatchCoordinates(Oligo query, int mismatches)
    {
        Oligo searchWindow = new Oligo();
        ArrayList hits = new ArrayList();
        final int QUERY_LENGTH = query.length();
        
        for(int i = 0; i + QUERY_LENGTH <= OLIGO_LENGTH; i++)
        {
            searchWindow = new Oligo(oligo.substring(i, i + QUERY_LENGTH));
            if(searchWindow.isFuzzyMatch(query, mismatches))
                hits.add(i);
        }
        return hits;
    }//end getAllMatchCoordinates()
    
    
    
    /**
     *Returns an ArrayList containing the indices of all occurences of the query
     *@param query query sequence
     *@param mismatches maximum number of allowedMismatches
     *@param ins maximum allowed number of inserts in SOURCE sequence
     *@param del maximum allowed number of deletes in SOURCE sequence
     *@param minKeyLength minimum length of found key in order for search to be considered successful
     *@return ArrayList object (ArrayList is empty if zero matches are found)
     *@see smithWaterman()
     *@see isFuzzSmithWatermanMatch()
     */
    public ArrayList getAllMatchCoordinates(Oligo query, int mismatches, int ins, int del, int minKeyLength)
    {
        Oligo searchWindow = new Oligo();
        Oligo swKey = new Oligo();
        ArrayList hits = new ArrayList();
        final int QUERY_LENGTH = query.length(); //net query length must include allowed number of inserts in search window
        
        for(int i = 0; i + QUERY_LENGTH + ins <= OLIGO_LENGTH; i++)
        {
            try
            {
                searchWindow = new Oligo(oligo.substring(i, i + QUERY_LENGTH + ins));   //searchWindow must compensate for number of allowed inserts
                swKey = searchWindow.smithWaterman(query, ins, del);
                searchWindow = new Oligo(oligo.substring(i, i + swKey.length()));       //narrow searchWindow to length of potential key
            }
            catch (OligoException e)
            {continue;}
            if(searchWindow.isFuzzySmithWatermanMatch(query, mismatches, ins, del, minKeyLength) && searchWindow.isFuzzyMatch(swKey, mismatches))
                hits.add(i);
        }
        return hits;
    }//end getAllMatchCoordinates()
    
    
    /**
     *Returns the index of the first occurence of query
     *@param query query sequence
     *@param mismatches maximum number of allowedMismatches
     *@return int index if found; otherwise returns -1
     */
    public int getFirstMatchCoordinate(Oligo query, int mismatches)
    {
        ArrayList list = getAllMatchCoordinates(query, mismatches);
        return (list.size() > 0) ? Integer.parseInt(list.get(0).toString()) : -1;
    }//end getFirstMatchCoordinate()
    
    
    
    /**
     *Returns the index of the first occurence of query
     *@param query query sequence
     *@param mismatches maximum number of allowedMismatches
     *@param ins maximum allowed number of inserts in SOURCE sequence
     *@param del maximum allowed number of deletes in SOURCE sequence
     *@param minKeyLength minimum length of found key in order for search to be considered successful
     *@return int index if found; otherwise returns -1
     */
    public int getFirstMatchCoordinate(Oligo query, int mismatches, int ins, int del, int minKeyLength)
    {
        ArrayList list = getAllMatchCoordinates(query, mismatches, ins, del, minKeyLength);
        return (list.size() > 0) ? Integer.parseInt(list.get(0).toString()) : -1;
    }//end getFirstMatchCoordinate()
    
    
    /**
     *Returns the found key (modified query) based on the given parameters. This method is based on
     *isFuzzySmithWatermanMatch(Oligo query, int mismatches, int ins, int del, int minKeyLength) and simply returns the key equivalent.
     *@param query query sequence
     *@param mismatches maximum number of allowedMismatches
     *@param ins maximum allowed number of inserts in SOURCE sequence
     *@param del maximum allowed number of deletes in SOURCE sequence
     *@param minKeyLength minimum length of found key in order for search to be considered successful
     *@return [Oligo object] Found key
     *@throws OligoException
     *@see also isFuzzySmithWatermanMatch(Oligo query, int mismatches, int ins, int del, int minKeyLength)
     */
    public Oligo getFuzzySWkey(Oligo query, int mismatches, int ins, int del, int minKeyLength) throws OligoException
    {
        if(isFuzzySmithWatermanMatch(query, mismatches, ins, del, minKeyLength))
            return smithWaterman(query, ins, del);
        else throw new OligoException(mnf, "getFuzzySWkey()");
    }//end getFuzzySWkey()
    
    
    
    /**
     *Returns the currenly set character to ignore
     *@return char character
     */
    public char getIgnoredChar()
    {
        return ignoredChar;
    }
    
    
    
    /**
     *Returns the start index of the last occurence of a match. If only one match found then result is same is getFirstMatchCoordinate()
     *@param query The query to search for
     *@param mismatches The maximum number of allowed mismatches
     *@return int representing the found index; otherwise returns -1
     */
    public int getLastMatchCoordinate(Oligo query, int mismatches)
    {
        ArrayList list = getAllMatchCoordinates(query, mismatches);
        return (list.size() > 0) ? Integer.parseInt(list.get(list.size() - 1).toString()) : -1;
    }
    
    
    
    /**
     *Returns the start index of the last occurence of a match. If only one match found then result is same is getFirstMatchCoordinate()
     *@param query query sequence
     *@param mismatches maximum number of allowedMismatches
     *@param ins maximum allowed number of inserts in SOURCE sequence
     *@param del maximum allowed number of deletes in SOURCE sequence
     *@param minKeyLength minimum length of found key in order for search to be considered successful
     *@return int representing the found index; otherwise returns -1
     */
    public int getLastMatchCoordinate(Oligo query, int mismatches, int ins, int del, int minKeyLength)
    {
        ArrayList list = getAllMatchCoordinates(query, mismatches, ins, del, minKeyLength);
        return (list.size() > 0) ? Integer.parseInt(list.get(list.size() - 1).toString()) : -1;
    }
    
    
    
    /**
     *Inserts the specified insert immediately BEFORE the specified index
     *@param insert sequence to be inserted
     *@param index The index where insertion is to occur
     *@return Oligo object with inserted sequence
     */
    public Oligo insert(Oligo insert, int index)
    {
        String insertSeq = insert.toString();
        StringBuffer targetOligo = new StringBuffer(oligo);
        return new Oligo(targetOligo.insert(index, insertSeq));
    }
    
    
    
    /**
     *Returns true if the input sequence contains only 'A/a', 'C/c', 'G/g', 'T/t' and the ignoredChar; Otherwise, returns false
     *@return boolean
     */
    public boolean isDNA()
    {
        char[] s = oligo.toCharArray();
        for(char base: s)
        {
            if(base != 'a' && base != 'c' && base != 'g' && base != 't' && base != ignoredChar)
                return false;
        }
        return true;
    }
    
    
    
    /**
     * Determines whether an oligo has zero length
     *@param none
     *@return boolean
     */
    public boolean isEmpty()
    {
        return (OLIGO_LENGTH == 0) ? true: false;
    }
    
    
    
    /**
     *Determines whether query sequence is contained in the oligo sequence, given the allowed mismatches. Will ignore N's. Method is not case sensitive
     *@param inputQuery The query to search for
     *@param mimatches - The maximum number of allowed mismatches
     *@return boolean
     */
    public boolean isFuzzyMatch(Oligo inputQuery, int mismatches)
    {
        char[] query = inputQuery.toCharArray();
        char[] source = oligo.toCharArray();
        int ntMatches = 0; //ntMatches = number of nucleotide-nucleotide (GC or AT) matches
        int ignoredCharMatches = 0; //ignoredCharMatches = number of nucleotide-ignoredChar (i.e. 'N') matches
        final int QUERY_LENGTH = query.length;
        
        for(int i = 0; i + QUERY_LENGTH <= OLIGO_LENGTH; i++)
        {
            for(int j = 0; j <= QUERY_LENGTH - 1; j++)
            {
                boolean nnMatch = source[i+j] == query[j];
                boolean niMatch = source[i+j] == ignoredChar || query[j] == ignoredChar;
                
                if(nnMatch && !niMatch) //count only nucleotide-nucleotide matches
                    ntMatches++;
                else if(niMatch)  //count only nucleotide-ignoredChar matches
                    ignoredCharMatches++;
            } //inner for() loop
            
            if(ntMatches >= QUERY_LENGTH - mismatches - ignoredCharMatches)
                return true;
            else
            {
                ntMatches = 0; //must reset for each unsuccesful loop
                ignoredCharMatches = 0; //must reset for each unsuccesful loop
            }
        } //outer for() loop
        return false; // exhausted all possibilities, no matches found
    } //end isFuzzyMatch() method
    
    
    
    /**
     *Combined isFuzzyMatch() and smithWaterman() to allow search that incorporates mismatches and/or indels
     *@param query The query sequence to search for
     *@param mismatches The maximum allowed number of mismatches between query and source
     *@param ins The maximum allowed number of inserts in the SOURCE sequence (i.e. # deletions in query)
     *@param del The maximum allowed number of deletions in the SOURCE sequence (i.e. # insertions in query)
     *@return boolean TRUE if query is found within source, given the specified conditions; FALSE otherwise
     *@throws OligoException
     *@see isFuzzyMatch()
     *@see smithWaterman()
     */
    public boolean isFuzzySmithWatermanMatch(Oligo query, int mismatches, int ins, int del)
    {
        Oligo swKey = new Oligo();
        
        try
        { swKey = smithWaterman(query, ins, del); }
        catch (OligoException error)
        { error.getStackTrace(); }
        
        return (isFuzzyMatch(swKey, mismatches)) ? true : false;
    } //end isFuzzySmithWatermanMatch()
    
    
    
    public boolean isFuzzySmithWatermanMatch(Oligo query, int mismatches, int ins, int del, int minKeyLength)
    {
        Oligo swKey = new Oligo();
        
        try
        {
            swKey = smithWaterman(query, ins, del);
        }
        catch (OligoException error)
        {
            error.getStackTrace();
        }
        
        return (isFuzzyMatch(swKey, mismatches) && swKey.length() >= minKeyLength) ? true : false;
    } //end isFuzzySmithWatermanMatch()
    
    
    
    /**
     *Returns the last index of the first match, if found. Uses isFuzzyMatch() algorithm
     *@param query The sequence to search for
     *@param mismatches The max number of allowed mismatches
     *@return int index (if found); otherwise, returns -1
     */
    public int lastIndex(Oligo query, int mismatches)
    {
        return (isFuzzyMatch(query, mismatches)) ? getFirstMatchCoordinate(query, mismatches) + query.length() - 1 : -1;
    }
    
    
    
    /**
     *Returns the length of oligo
     *@param void
     *@return int - length of oligo
     */
    public int length()
    {
        return OLIGO_LENGTH;
    }
    
    
    
    //------------------ligate the 5' (upstream) end of oligo2 to the 3' (downstream) end of oligo--------------------//
    public Oligo ligate(Oligo input)
    {
        return new Oligo(oligo + input.toString());
    }
    
    
    
    /**
     *Induces base substitution mutations in oligo
     *@param percent The percentage of mutation to be induced
     *@return Oligo object - mutated oligo sequence
     */
    public Oligo mutate(int percent)
    {
        int probability = (int)(Math.random() * (100+1)); //probability that the oligo will be mutated
        final int LAST_INDEX = oligo.length() - 1;
        StringBuffer tempOligo = new StringBuffer(oligo);
        
        if(probability > percent)
            return new Oligo(tempOligo);
        else if(probability <= percent)
        {
            int index = (int)(Math.random() * LAST_INDEX);
            char ch = generateRandomBase();
            if(oligo.charAt(index) != ch)
                tempOligo.setCharAt(index, ch);
        }
        return new Oligo(tempOligo);
    } //end mutate()
    
    
    
    /**
     *Randomizes the given oligo sequence
     *@return [Oligo object] randomized oligo sequence
     */
    public Oligo randomizeOligo()
    {
        StringBuffer tempOligo = new StringBuffer(oligo);
        
        for(int i = 0; i <= 2 * oligo.length(); i++)
        {
            //randomly pick a character and move it to the end
            int randomIndex = (int)(Math.random()*(oligo.length() - 1));
            char tempChar = tempOligo.charAt(randomIndex);
            tempOligo.deleteCharAt(randomIndex);
            tempOligo.append(tempChar);
        }
        return new Oligo(tempOligo);
    } //end randomizeOligo()
    
    
    
    /**
     *Resets the size of an oligo to zero length
     *@return Oligo object: oligo of zero length
     */
    public Oligo resetToNull()
    {
        return new Oligo();
    }
    
    
    
    /**
     *Returns the reverse of the oligo
     *@return Oligo object: the reversed oligo
     */
    public Oligo reverse()
    {
        StringBuilder tempOligo = new StringBuilder(oligo);
        return new Oligo(tempOligo.reverse());
    }
    
    
    
    /**
     *Sets the character to ignore during manipulation operations. Default is 'N'
     */
    public void setIgnoredChar(char inputChar)
    {
        this.ignoredChar = inputChar;
    }
    
    
    
    /**
     *Returns a modified query sequence based on successful of search. Allows user to specify maximum number of insertions and deletions in the source sequence
     *@param query The query sequence to search for
     *@param ins The maximum number of inserts allowed in the source sequence
     *@param del The maximum number of deletes allowed in the source sequence
     *@return Oligo object: transformed query sequence, if within specified conditions. Otherwise, returns "smithWaterman(): NotFoundException"
     *@throws OligoException
     */
    public Oligo smithWaterman(Oligo query, int ins, int del) throws OligoException
    {
        char[] s = ("x" + oligo).toCharArray(); //pad 'x' as first char of oligo
        char[] q = ("x" + query).toCharArray(); //pad 'x' as first char of query
        float[][] matrix = new float[q.length + 1][s.length + 1];
        float max = 0;
        float best = 0; //highest score in each comparison
        float iScore = 0;
        float jScore = 0;
        float diagScore = 0; //nucleotide-nucleotide alignment score: match = 1.0; mismatch = -0.3
        int imax = 0;
        int jmax = 0;
        
        //initialize row 0 to 0.0
        for(int i = 0; i <= q.length; i++)
            matrix[i][0] = (float)0.0;
        
        //initialize column 0 to 0.0
        for(int j = 0; j <= s.length; j++)
            matrix[0][j] = (float)0.0;
        
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
                    jScore = matrix[i][j-hGap] - (float)(1 + 0.3*hGap);
                    best = (jScore > best) ? jScore : best;
                }
                
                matrix[i][j] = (best > 0) ? best : 0;
                if(matrix[i][j] > max)
                {
                    max = matrix[i][j];
                    imax = i;
                    jmax = j;
                }
            } //for(int j...) loop
        } //for(int i...) loop
        
        //reconstruct alignment
        int nDel = 0;               //number of deletions in source
        int nIns = 0;               //number of insertions in source
        float above = 0;            //score of cell matrix[imax - 1][jmax]
        float left = 0;             //score of cell matrix[imax][jmax - 1]
        float diag = 0;             //score of cell matrix[imax - 1][jmax - 1]
        boolean isAbove = false;    //best score is cell above
        boolean isDiag = false;     //best score is diagonal cell
        boolean isLeft = false;     //best score is cell to left
        StringBuilder key = new StringBuilder();
        key = key.append(q[imax]);  //initialize key to ch associated with max score
        
        do
        {
            //System.out.println("imax = " + imax + " jmax = " + jmax + "  matrix[imax][jmax] = " + matrix[imax][jmax]);
            above = matrix[imax-1][jmax];
            diag = matrix[imax-1][jmax-1];
            left = matrix[imax][jmax-1];
            isAbove = (above >= diag) && (above >= left);
            isDiag = (diag >= above) && (diag >= left);
            isLeft = (left >= diag) && (left >= above);
            
            if(isAbove && !isDiag) //deletion in source sequence
            {
                imax--;
                nDel++;
            }
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
                nIns++;
            }
        } while(matrix[imax][jmax] > 0);
        
        //evaluate success of search based upon specified conditions
        if(nDel <= del && nIns <= ins)
            return new Oligo(key.reverse());
        else throw new OligoException(mnf, "smithWaterman()");
    } //end smithWaterman() method
    
    
    
    /**
     *Deletes the first occurence of the target sequence based on a best fuzzy-match search. Overloaded to include the option
     *of including indels
     *@param target
     *@param mismatches The maximum number of allowed mismatches
     *@param indel The maximum number of allowed insertions and/or deletions
     *@return Oligo object: spliced oligo, if target is found
     */
    public Oligo spliceOut(Oligo query, int mismatches) throws OligoException
    {
        if(isFuzzyMatch(query, mismatches))
        {
            int index = getFirstMatchCoordinate(query, mismatches);
            return new Oligo(oligo.substring(0, index) + oligo.substring(index + query.length(), oligo.length() - 1));
        }
        else throw new OligoException(mnf, "spliceOut()");
    } //end spliceOut()
    
    
    
    /**
     *Splices out a given sequence based on specified start and end indices, inclusively
     *@param start start index
     *@param end end index
     *@return Oligo object representing remaining sequence after splicing
     *@throws OligoException
     */
    public Oligo spliceOut(int start, int end) throws OligoException
    {
        String s = "";
        boolean startIndexOutofBounds = start < 0;
        boolean endIndexOutofBounds = end > oligo.length() - 1;
        
        if(start >= 0 && end <= oligo.length() - 1)
            s = oligo.substring(0, start) + oligo.substring(end + 1, oligo.length() - 1);
        else if(startIndexOutofBounds ^ endIndexOutofBounds)
        {
            if(startIndexOutofBounds)
                throw new OligoException(soob, "extractSequence()");
            else throw new OligoException(eoob, "extractSequence()");
        }
        else if (startIndexOutofBounds && endIndexOutofBounds)
            throw new OligoException("(Start&End)IndexOutofBounds", "extractSequence()");
        return new Oligo(s);
    } //end spliceOut()
    
    
    
    //----------------------------parseFAST(): extracts id tag and sequence from a FASTA file-------------------------//
//    public static String[] readFromFASTA(File inputFile) throws IOException
//    {
//        String[] element = new String[2];
//        InputStream fin = new BufferedInputStream(new FileInputStream(inputFile));
//
//        //extract the id only
//        String id = new String();
//        char ch = (char)fin.read();
//
//        if(ch == '>')
//        {
//            do
//            {
//                id += ch;
//            } while((ch = (char)fin.read()) != '\r');
//        }
//        else
//        {
//            fin.close();
//            System.out.println("Tag not found. Program terminated");
//            return element;
//        }
//
//        //if current ch read is a CR, then title is found ==> proceed to extract the DNA sequence
//        String seq = new String();
//        for(int i = 0; fin.available() > 0; i++)
//        {
//            if((ch = (char)fin.read()) != '\r')
//                seq += ch;
//        }
//
//        fin.close(); //close current input file
//        element[0] = id + '\r';
//        element[1] = seq;
//        return element;
//    } //end parseFASTA()
////////////////////////////////////////////////////////////////////////////////////////
} //end Oligo class
