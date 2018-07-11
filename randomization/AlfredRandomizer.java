import java.io.BufferedReader;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import java.io.InputStreamReader;
import java.io.RandomAccessFile;
import java.nio.ByteBuffer;
import java.nio.channels.FileChannel;
import java.nio.channels.FileLock;
import java.nio.channels.OverlappingFileLockException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.BitSet;
import java.util.Formatter;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Locale;
import java.util.Map;
import java.util.Random;
import java.util.TreeMap;
import java.util.regex.Matcher;
import java.util.regex.Pattern;
import java.util.zip.GZIPInputStream;

/**
 * A randomization test to determine genes with increased burden of variants in cases, compared to controls; 
 * alternatively, performs the ALFRED test for co-occurence of damaging variants and LOH events in the same genes. <p>
 * 
 * Can account for population structure and/or batch effects due to technical biases by randomizing within 
 * clusters of patients, which have to be defined beforehand, typically by running PCA on a matrix of common variants. <p>
 * 
 * See AlfredRandomizer.main() for description of parameters and examples of usage. <p>
 * 
 * See reference: Park, Supek and Lehner (2018) Nature Communications 
 * "Systematic discovery of germline cancer predisposition genes through the identification of somatic second hits".
 * https://www.nature.com/articles/s41467-018-04900-7
 * 
 * @author Fran Supek, CRG.
 */
public class AlfredRandomizer {

  public int numGenes;
  
  /**
   * Names of the genotype-containing table columns, ordered as given in the input file. 
   * The ordering is the same as the ordering in the outer array in genotypes.
   */
  public List<String> genotypeNames = new ArrayList<>();
  
  /**
   * Indexed first by gene, then by patient. (handy for randomizing within each gene).
   * Typically contains 1.0 for mutant or 0.0 for wildtype or Float.NaN for unknown.
   * May contain also continuous values, eg. the LRT scores which weigh each SNP
   * differently, depending on its deleteriousness.
   */
  public float[][] genotypes;
  
  /**
   * Indexed by cluster # (see this.clusterNamesList); the inner array contains the indices of the individuals within that cluster.
   * This is very handy for the purposes of randomization. For lookup by patient, take "clusterMemberIndicesRepacked".
   */
  public int[][] clusterMemberIndices;
  
  /** 
   * The ordering is important - it corresponds to cluster indices in this.clusterMemberIndices.  The ordering is 
   * established by the order of occurrence in the input file -- the clusters are currently NOT sorted by ID.
   * Zero-based indices!
   */
  public List<String> clusterNamesList;
  
  /** Created by copying clusterMemberIndices; indexed by individuals, stores the cluster index for each individual. */
  public int[] clusterMemberIndicesRepacked;
  
  /** May be directly supplied by user, or read from file (if user provided an index). */
  public String clusterColumnName;
  /** May be directly supplied by user, or read from file (if user provided an index). */
  public String classColumnName;
  
  
  /**
   * Ordered as supplied in input table, and useful to retrieve patient index for this.genotypes.
   * Note: the patients that had class value other than "0" or "1" are NOT included here - they're skipped.
   */ 
  public List<String> patients = new ArrayList<>();
  
  
  /**
   * Contains TRUE for positive class (cancer patients) and false for controls.
   * May NOT contain a missing value under any circumstances - these cases are simply not loaded from input file.
   */
  public BitSet classLabels = new BitSet();
  
  /**
   * Number of patients; includes both positive class (cancer/disease) and control set (healthy) people;
   * or equivalently the LOH and no-LOH samples.
   */
  public int numPatients;
  
  
  /**
   * Parses a tab-separated text table with input data, preparing for randomization later.
   * 
   * This table describes one patient per row, and it must contain a header row. The parameters below allow the columns
   * to be recognized.
   * 
   * @param patientIdColumnName One column with patient ID (which may contain any kind of string; not used in calculations).
   * 
   * @param classColumnNameOrIdx A single column containing class labels; Typically this means "cancer patients" 
   * (encoded as 1) versus healthy controls (encoded as 0). 
   * 
   * Alternatively, the 1s and the 0s may be used to compare one cancer type to others, or young cancer patients to old,
   * or similar. In all cases, "1" are the positive cases, and "0" are the controls. If another value is present here, 
   * this row will be completely ignored.
   * 
   * @param clusterColumnNameOrIdx A single column containing the clustering of the data points; clusters may be given 
   * as any String. 
   * 
   * IMPORTANT: rows containing a "?" or "NA" for the cluster ID will be ignored!
   * 
   * The clustering is typically done on the principal components of the common variants, and is meant to reflect
   * population structure and (possibly) technical biases.
   * 
   * @param geneColumnsRegexOrIdx A regular expression that matches many columns, one per examined gene (possibly a 
   * pathway/complex), stating if a gene is mutated in that specific patient. <p>
   * 
   * Alternatively, if you supply a SINGLE integer, will not treat this as a regex but will instead just examine this 
   * one gene column (marked with this 1-based index). <p>
   * 
   * Finally, you can supply a range of integers, eg. 55-170, then it will test all columns in range. <p>
   * 
   * These genotype columns must contain numbers: typically 0 for not mutated, and 1 for mutated), or "?" for a missing
   * value; they may also contain fractional numbers in case we use LRT scores.<p>
   */
  public AlfredRandomizer(String filePath,
          String patientIdColumnName, 
          String classColumnNameOrIdx,
          String clusterColumnNameOrIdx,
          String geneColumnsRegexOrIdx
          //String princompColumnsRegex
          ) throws IOException {
    
    // store indices of columns
    int classCol = -1, idCol = -1, clusterCol = -1;
    List<Integer> geneCols = new ArrayList<Integer>(1000);
    // List<Integer> pcCols = new ArrayList<Integer>(1000);
    
    // temporary storage for genotypes, first indexed by patient, and then
    // contains genotypes for all genes in the patient
    List<float[]> tempGeno = new ArrayList<>();
    
    // temporary storage for pcs, first indexed by patient, and then contains all PCs in the patient
    // List<float[]> tempPcs = new ArrayList<>();
    
    // cluster membership for every patient (will later be converted to an array of cluster indices)
    List<String> tempClusterPerPatient = new ArrayList<>(); 
    Map<String, Integer> tempClusterSize = new LinkedHashMap<>();  // IMPORTANT: LinkedHashMap, so it keeps the order in which they appear in input file
    
    BufferedReader br;
    if ( filePath.toLowerCase().endsWith(".gz") ) {
      GZIPInputStream gzip = new GZIPInputStream(
              new FileInputStream(filePath) );
      br = new BufferedReader(new InputStreamReader(gzip));
    } else {
      br = new BufferedReader(new FileReader(filePath));
    }    
    String line;
    boolean isHeader = true;
    int curPatientIdx = 0;
    
    while ( (line = br.readLine()) != null ) {
      line = line.trim();
      if ( line.length() == 0 ) {
        continue;
      }
      String[] cols = line.split("\t"); // this may be huge
      if ( isHeader ) {
        isHeader = false;
        List<String> colsList = Arrays.asList(cols);
        idCol = colsList.indexOf(patientIdColumnName);
        
        //int classCol;
        try {
          classCol = Integer.parseInt(classColumnNameOrIdx)-1;  // user supplied a 1-based index
        } catch (NumberFormatException numberFormatException) {
          classCol = colsList.indexOf(classColumnNameOrIdx);    // .. or a column name
        }
        try {
          clusterCol = Integer.parseInt(clusterColumnNameOrIdx)-1;  // user supplied a 1-based index
        } catch (NumberFormatException numberFormatException) { 
          clusterCol = colsList.indexOf(clusterColumnNameOrIdx);  // .. or a column name
        }
        this.classColumnName = colsList.get(classCol);
        this.clusterColumnName = colsList.get(clusterCol);
        
        int tryGeneCol;
        try {
          tryGeneCol = Integer.parseInt(geneColumnsRegexOrIdx);  // did we supply a single integer? (1-based column index)
          geneCols.add(tryGeneCol-1);  // in this case, this will be the only column!  (add 0-based index to the list)
          this.genotypeNames.add( String.valueOf(cols[tryGeneCol-1]) );
        } catch (NumberFormatException numberFormatException) {
          tryGeneCol = -1;    // we supplied a standard regex; will try to match it now
        }
        
        if ( geneCols.isEmpty() ) {  // is it maybe a range of columns?
          Matcher mColRange = Pattern.compile(geneColumnsRegexOrIdx).matcher("(\\d+)-(\\d+)");
          if ( mColRange.matches() ) {
            for ( int cc = Integer.parseInt(mColRange.group(1)); cc <= Integer.parseInt(mColRange.group(2)); cc++ ) {
              geneCols.add( cc-1 );  // this list contains 0-based indices
            }
          }
        }
        
        if ( geneCols.isEmpty() ) {  // it's not a single column nor a range, try a regex
          Matcher mGene = Pattern.compile(geneColumnsRegexOrIdx).matcher("blah");
          // Matcher mPc = Pattern.compile(geneColumnsRegex).matcher("blah");
          for ( int i=0; i<cols.length; i++ ) { // examine header, column by column
            if ( i==idCol || i==classCol || i==clusterCol ) {
              continue;   
            }
            mGene.reset(cols[i]);
            if ( mGene.matches() ) {
              geneCols.add(i);
              this.genotypeNames.add( String.valueOf(cols[i]) );
            }
            //mPc.reset();
            // if ( mPc.matches() ) pcCols.add(i);
          }
        }
        
        
        if ( classCol==-1 || clusterCol==-1     // note: we don't enforce that the patientID column must be found - we don't use it anyway
                || geneCols.isEmpty() //|| pcCols.isEmpty()
                ) {
          throw new IllegalArgumentException("Could not find some of the necessary columns in the table.");
        }
        
        this.numGenes = geneCols.size();
        // this.numPcs = pcCols.size();
        continue;
      }  // if is header?
      
      // ......... conditions for excluding rows
      // (a) if the 'class' is undefined (given as anything but 1 or 0), will skip this row completely
      if ( !cols[classCol].equals("1")  && !cols[classCol].equals("0") ) continue;
      // (b) rows with the cluster assigned to "NA" or "?" are completely skipped!
      if ( cols[clusterCol].equals("NA") || cols[clusterCol].equals("?")
           // !cols[clusterCol].equals("1") // @@@ DEBUG - DON'T USE @@@
              ) continue;
      
      this.classLabels.set(curPatientIdx, cols[classCol].equals("1"));
      if ( idCol > -1 )
        this.patients.add(String.valueOf( cols[idCol] ));
      else
        this.patients.add("noPatientIdDefined");
      float[] genoForPatient = new float[this.numGenes];   // loads row by row (patient by patient)
      for ( int i=0; i<this.numGenes; i++ ) {
        try {
          genoForPatient[i] = Float.parseFloat(cols[ geneCols.get(i)]);
        } catch (NumberFormatException numberFormatException) {
          genoForPatient[i] = Float.NaN;  // we allow missing values in genotypes
        }
      }
      tempGeno.add(genoForPatient);
      String aClust = String.valueOf( cols[clusterCol] );
      tempClusterPerPatient.add(aClust);
      if ( !tempClusterSize.containsKey(aClust) )
        tempClusterSize.put( aClust, 1 );
      else 
        tempClusterSize.put( aClust, tempClusterSize.get(aClust) + 1 );
      
      //double[] pcsForPatient = new double[this.numPcs];
      // for ( int i=0; i<this.numPcs; i++ ) pcsForPatient[i] = Double.parseDouble(cols[ pcCols.get(i) ] );
      //this.pcs.add(pcsForPatient);
      curPatientIdx++;
    }  // ... line by line in input file
    this.numPatients = curPatientIdx;
    
    // transpose the matrix with the genotypes, so that they're
    // first indexed by gene -- because we'll randomize across patients, for speed
    this.genotypes = new float[ this.numGenes ][ this.numPatients ]; 
    for ( int i=0; i<this.numPatients; i++ ) {
      for ( int j=0; j<this.numGenes; j++ ) {
        this.genotypes[j][i] = tempGeno.get(i)[j];
      }
    }
    
    // ... create the 2D matrix of cluster membership indices; indexed by cluster #, the inner array contain
    // the patient indices (positive or negative) that belong to that cluster 
    // Set<String> tempClusterSet = new HashSet<>( tempClusterPerPatient );  // a set of unique cluster IDs ... 
    
    // transform the set of unique cluster IDs to a list to give them ordering
    this.clusterNamesList = new ArrayList<>( tempClusterSize.keySet() );  
    this.clusterMemberIndices = new int[this.clusterNamesList.size()][];

    for (int i = 0; i < clusterNamesList.size(); i++) {
      this.clusterMemberIndices[i] =  
              new int[ tempClusterSize.get( this.clusterNamesList.get(i) ) ];  // how many members in each cluster
    }
    
    int[] tempClusterI = new int[ this.clusterNamesList.size() ];  // where we are (current idx) in each cluster
    
    // store also the clusterMemberIndices indexed per individual, storing indices of clusters (instead of vice versa as in clusterMemberIndices).
    this.clusterMemberIndicesRepacked = new int[ tempClusterPerPatient.size() ];
    
    for ( int patientIdx=0; patientIdx<tempClusterPerPatient.size(); patientIdx++ ) {
      
      int whichClust = this.clusterNamesList.indexOf(
              tempClusterPerPatient.get(patientIdx) );  // this is a bit slow, but at least it's a short list
      this.clusterMemberIndices[ whichClust ][  tempClusterI[whichClust] ] = patientIdx;
      tempClusterI[whichClust]++;  // advance the counter for this cluster
      
      this.clusterMemberIndicesRepacked[ patientIdx ] = whichClust;  // also, fill an array which is indexed by cluster (not by patient)
      
    }
    
    // just a check .. 
    for (int i = 0; i < clusterNamesList.size(); i++) {
      if ( tempClusterI[i] != tempClusterSize.get( clusterNamesList.get(i) ) )
        throw new IllegalArgumentException("Whoops!");  // this should not be a problem if code is bug-free
    }
    
  }
  
  
  
  /**
   * Calculates the difference of 'average' genotype value for the positive class (cases) versus the 
   * negative class (controls). Genotypes are typically given as 0 and 1 - then, this test statistic amounts
   * to difference in frequencies - but continuous values for genotype are also supported here. <p>
   * 
   * Genotypes marked with Float.NaN are completely ignored. <p>
   * 
   * "clusterMemberIndices" Is ignored in this particular function. If you require per-cluster statistics, use
   * the overloaded variant of this function  calcTestStatistic( float[], BitSet, int[][] ). 
   * 
   * @param genotypes For a single gene, across all patients.
   * @param classLabels "True" for positive class, "false" for negative. (this may be the actual, observed labels, 
   * or also randomized ones).
   * 
   * @return 
   */
  public static double calcTestStatistic( float[] genotypes, BitSet classLabels ){ 
    
    int numPos=0, numNeg=0;
    double sumPos=0.0, sumNeg=0.0;
    for (int i = 0; i < genotypes.length; i++) {
      float f = genotypes[i];
      if ( Float.isNaN( f ) )
        continue;
      if ( classLabels.get(i) ) {  // positive class (cancer patients)
        numPos++;
        sumPos += f;  // this will be +0.0 for persons not having the gene mutated
      } else { // negative class (healthy controls)
        numNeg++;
        sumNeg += f;
      }
    }
    return (  (  sumPos / numPos  ) - ( sumNeg / numNeg ) );
    
 }
  
  
  public static class RandOutStats {
    final double testStat;
    final int numPos, numNeg;
    final double sumPos, sumNeg;
    public RandOutStats(double testStat, double sumPos, int numPos, double sumNeg, int numNeg) {
      this.testStat = testStat;
      this.numPos = numPos; this.numNeg = numNeg;
      this.sumPos = sumPos; this.sumNeg = sumNeg;
    }
  }
  
  
  
  /**
   * Calculates the difference of 'average' genotype value for the positive class (cases) versus the 
   * negative class (controls). Determined for each cluster separately. <p>
   * 
   * Genotypes are typically given as 0 and 1; in that case then, this test statistic denotes the 
   * difference in frequencies between groups.  However, fractional values for the genotype are also supported here
   * allowing a more general burden test. <p>
   * 
   * Genotypes marked with Float.NaN (missing values in the data) are completely ignored.
   * 
   * @param genotypes For a single gene, across all patients.
   * @param classLabels "True" for positive class, "false" for negative. (this may be the actual, observed labels, 
   * or also randomized ones).
   * @param clusterMemberIndices The outer array is indexed by the cluster #; the inner array contains indices of individuals.
   * 
   * @return Array, indexed by the cluster_# + 1; the index[0] corresponds to the overall data (all clusters combined;
   * as would be obtained from the standard function calcTestStatistic(float[], BitSet). <p>
   */
  public static RandOutStats[] calcTestStatisticsPerClust( float[] genotypes, BitSet classLabels, int[][] clusterMemberIndices ){ 
    
    // overall sums, for all clusters (they will not include instances that are not a member of any cluster - this is not supported)
    int numPos=0, numNeg=0;
    double sumPos=0.0, sumNeg=0.0;
    
    int numClust = clusterMemberIndices.length;  // the outer array is indexed by the cluster #
    
    RandOutStats[] result = new RandOutStats[ clusterMemberIndices.length + 1 ];  // note "+1" - the last member contains the summary stats for all clusters
    
    // operate within each cluster separately
    for ( int i=0; i < clusterMemberIndices.length; i++ ) { 
      
      int[] clusterMembers = clusterMemberIndices[i];
      
      int numPosCl=0, numNegCl=0;
      double sumPosCl=0.0, sumNegCl=0.0;
      
      for ( int aPatient : clusterMembers ) {
         
        float f = genotypes[ aPatient ];
        if ( Float.isNaN( f ) ) continue;
        
        if ( classLabels.get(aPatient) ) {  // positive class (cancer patients)
          numPosCl++;
          sumPosCl += f;  // this will be +0.0 for persons not having the gene mutated
        } else { // negative class (healthy controls)
          numNegCl++;
          sumNegCl += f;
        }        
        
      } // loop over patients in *this* cluster
      
      double testStatCl = ( sumPosCl / numPosCl ) - ( sumNegCl / numNegCl );
      result[i+1] = new RandOutStats(testStatCl, sumPosCl, numPosCl, sumNegCl, numNegCl);
      numPos += numPosCl; sumPos += sumPosCl;
      numNeg += numNegCl; sumNeg += sumNegCl; 
      
    } // cluster by cluster
    
    // int i=clusterMemberIndices.length; // the last index -- sums across all clusters 

    // [0][] stores the pooled data for ALL clusters
    double testStat = (  sumPos / numPos  ) - ( sumNeg / numNeg );
    result[0] = new RandOutStats(  testStat, sumPos, numPos, sumNeg, numNeg );
    
    return result;
    
 }  
  
  
  
  
  /**
   * Counts the total number of variants (or sum of weights) in either the cases (positive class) 
   * or the controls (negative class). <p>
   * 
   * Note that the 'patients' vs 'healthy' may actually be made to refer to one 
   * cancer type versus another (depending on how +ve and -ve classes were defined). <p>
   * 
   * Actually computes the sum of genotypes - if they're encoded as 0/1, it's a count
   * of 1s. (they could be fractional quantities as well in principle). <p>
   * 
   * If a genotype is NaN, skips it - it counts neither towards positive nor negative class. <p>
   * 
   * @param genotypesForGene For a single gene, across all patients.
   * @param countPatients "True" to count positive class, "false" to count negative class.
   * @return A floating point number, since genotypes can also be fractional.
   */
  public double countVariantsForGene( float[] genotypesForGene, boolean countPatients ) {
    
    double sumPos=0.0, sumNeg=0.0;
    for (int i = 0; i < genotypesForGene.length; i++) {  // loop over patients
      float f = genotypesForGene[i];
      if ( Float.isNaN( f ) )
        continue;
      if ( classLabels.get(i)  ) {  // positive class (cancer patients)
        sumPos += f;  // this will be +0.0 for persons not having the gene mutated
      } else { // negative class (healthy controls)
        sumNeg += f;
      }
    }
    return countPatients ? sumPos : sumNeg;
    
  }
  
  
  
  /**
   * Alters the supplied BitSet by randomizing the positive (cases) vs negative (controls) labels. <p>
   * 
   * Observes the cluster memberships stored in this.clusterMemberIndices, taking care 
   * to randomize only within a cluster. This is meant to control for population stratification.
   * 
   * @param toShuffle WARNING - alters this bitset permanently! Remember to clone before randomizing.
   * @param r Recommended to use the included XorShiftRandom (fast). If running multithreaded, each thread
   * should have its own random generator.
   * @return A new bitset, which is a randomization of the existing one.
   */
  public void randomizeLabels( BitSet toShuffle, Random r ) {
    
    for (int[] clusterMembers : this.clusterMemberIndices) {  // operate within clusters
      
      // exchange the class label of this patient, with another patient from the same cluster
      for (int aPatient : clusterMembers) {
        // do this for all patients, thus shuffling the cluster
        int swapWith = clusterMembers[r.nextInt(clusterMembers.length)];
        if ( toShuffle.get(aPatient) != toShuffle.get(swapWith) ) {
          toShuffle.flip(aPatient);   // exchange the two (without using a helper variable)
          toShuffle.flip(swapWith);
        }
        //boolean temp = toShuffle.get( aPatient );
        //toShuffle.set( aPatient, toShuffle.get( swapWith ) );
        //toShuffle.set( swapWith, temp );
      }
      
    }
    
  }
  
  
  
  /**
   * Contains the results of many randomizations for a single gene.
   */
  static class GeneRandomizationResult {
    
    /** Used only for printing. */
    String genotypeName; 
    /**
     * In how many randomizations the test statistic is higher or equal than the actual value.
     * Normally, pRight + pLeft is slightly larger than 1.0; the cases where the 
     * randomized test statistic is == to the observed test statistic count in both left and right tail.
     */
    double pRight;
    /** In how many randomizations the test statistic is lower or equal than the actual value. */
    double pLeft;
    /**
     * The pRight, corrected for multiple testing by the Westfall & Young maxT algorithm. Will be NaN if the correction was not performed.
     * CURRENTLY NOT SUPPORTED/USED.
     */
    double pRightMaxT = Double.NaN;
    /** The observed value of the test statistic. */
    double observed;
    /** The observed value of the test statistic and its constituent components, for all clusters (in [0]) and additionally for each cluster separately. */
    RandOutStats[] observedAllStats;
    /** The minimum, maximum and average and median of the randomized test statistic. */
    double randMin, randMax, randAvg, randMedian;
    /**
     * The lower bounds (inclusive) of 100 bins in the histogram of randomized test statistic.
     * The upper bound of the highest (100th) bin is given in the 101th index (this correspond to this.randMax).
     */
    double[] randHistoBinsLowBound = new double[101];
    /** The histogram of randomized test statistic is fixed to contain 100 bins. */
    int[] randHistoCounts = new int[100];
    
    /**
     * Always of length 101; stores the n-th percentile of the randomized values in each field.
     * Useful to obtain confidence intervals.
     */
    double[] randPercentiles = new double[101]; 
    
    int nIter;
    
    /**
     * The randomized values of the test statistic, index is the # randomizations. <p>
     * 
     * This field is normally null to save memory, but is not null if we run the randomizeOneGene() with the switch
     * keepRandomized=TRUE. Then it could be used for a multiple testing correction based on the distribution of the
     * test statistic. It could also be used to output the details of individual randomizations.
     */
    double[] randomizedStat = null;
    
    /**
     * The randomized values of the test statistic, outer index is the # randomizations, inner index is the #clusters+1.
     * Indexed under [][0], the stats are stored for all the clusters pooled, while the individual clusters are at [][1], [][2] etc.s
     */
    RandOutStats[][] randomizedAllStats = null;
    
    /**
     * The count of variant alleles in each of the two groups of people that are compared.
     * Actually, this is the sum of variant weights, in case you use the fractional variants.
     * PosClass are typically the cancer patients, NegClass are the healthy controls. (although this could
     * be a comparison between one and other cancer types).
     */
    double nVariantsPosClass, nVariantsNegClass;
    
    /** Makes a tab-separated header for a text table that will contain GeneRandomizationResult as rows. */
    public static String getTabSepHeader( ) {
      StringBuilder sb = new StringBuilder();
      sb.append("nVariantPosClass\tnVariantNegClass\t");
      sb.append("pValRight\tpValLeft\t");
      // sb.append("maxTCorr_pValRight\t" );
      sb.append("excessOverRandMedian\texcessOverRand95CI\t");
      sb.append("observed\trandMedian\trandMin\trandMax");
      // the 0th percentile = minimum; the 100th percentile = maximum
      for ( int i=0; i<=100; i++ ) sb.append( String.format("\trandPerc_%03d", i ) );
      for ( int i=0; i<100; i++ ) sb.append( String.format("\trandHist_%03d", i ) );
      return sb.toString();
    }

    /** 
     * Prints out a tab-delimited row with all the statistics, including the distribution (both as percentiles
     * and as a histogram with 100 fixed-width bins. Uses the US locale (decimal point, not comma).
     * 
     * @return Does not contain a newline character.
     */
    public String toString( ) {
      Formatter fm = new Formatter(Locale.US);
      fm.format( "%.2f\t%.2f\t",
              this.nVariantsPosClass, this.nVariantsNegClass
              );
      // fm.format( "%.2e\t%.2e\t%.2e\t", this.pRight, this.pLeft, this.pRightMaxT);
      fm.format( "%.2e\t%.2e\t", this.pRight, this.pLeft);
      
      fm.format( "%.4f\t[%.4f..%.4f]\t",    // the excess frequencyDifference (between tumor and healthy), or whatever other statistic you have, over the randomized case
              this.observed-this.randMedian, 
              this.observed-(this.randPercentiles[2]+this.randPercentiles[3])/2.0,   // and the 95% C.I. derived from percentiles
              this.observed-(this.randPercentiles[97]+this.randPercentiles[98])/2.0
              );
      fm.format( "%.4f\t%.4f\t%.4f\t%.4f",
              this.observed, this.randMedian, this.randMin, this.randMax
              );
      for ( int i=0; i<=100; i++ ) fm.format( "\t%.4f", this.randPercentiles[i] );
      for ( int i=0; i<100; i++ ) {  // now histogram, formatted as text fields
        fm.format( "\t[%.2e..%.2e" + (i==99 ? "]" : ">") + "_%d",
                this.randHistoBinsLowBound[i], this.randHistoBinsLowBound[i+1], this.randHistoCounts[i]  );
      }
      return fm.toString();
    }
    
  }
  
    
  
  /**
   * Performs many randomizations, and finds the test statistic for the original data & for randomizations,
   * as well as an empirical p-value. Outputs a histogram of the randomized test statistic, in 100 equal-width bins
   * that span the whole range of the randomized values
   * 
   * @param genotypes For a single gene, across many patients. Typically 0s (not mutated)
   * and 1s (mutated), although may also contain fractional values.
   * @param classLabels A bitset containing 'true' for positive class (cancer patients) and
   * 'false' for negative class.
   * @param nIter How many iterations. For the Bonferroni corrected p of 2.5e-6 (0.05/20000 genes)
   * we need 5*2.5e6 = 1.25e7 (13 million) randomizations. Thus recommended 20_000_000.
   * The maximum you can ask is what fits in a Java int (2e9, 2 thousand million).
   * @param randSeed Creates a new XorShiftRandom using this seed.
   * 
   * @param keepRandDistInMem Stores the distribution randomized values of the test statistic in memory - 
   * uses a bit more RAM. If this is used, enables the Westfall & Young (1993) maxT algorithm for multiple testing
   * correction to be used, or also the statistics of individual randomization interations to be output. <p>
   * 
   * This is stored in the result.randomizedStat, which keeps the distribution of the test statistic
   * (difference in frequencies) and also the individual components: [][0] is the test statistic, 
   * (  sumPos / numPos  ) - ( sumNeg / numNeg ); [][1] is the sumPos; [][2] is the numPos; 
   * [][3] is the sumNeg; [][4] is the numNeg.
   * 
   * @return A GeneRandomizationResult object, which is quite compact in memory (it does not contain
   * the entire randomized distribution).
   */
  public GeneRandomizationResult randomizeOneGene( float[] genotypes, BitSet classLabels, int nIter, long randSeed,
          boolean keepRandDistInMem ) {

    GeneRandomizationResult result = new GeneRandomizationResult();
    result.nIter = nIter;
    
    double[] randTestStat = new double[nIter];  // huge array, discarded after running this function (needed for histogram)
    
    if ( keepRandDistInMem )  {
      result.randomizedAllStats =  new RandOutStats[nIter][];   // huge array, kept after running this function only if neccesary (for output
    }
    // (or not discarded if keepRandomized!!)
    
    result.observed = calcTestStatistic(genotypes, classLabels); 
    result.observedAllStats = calcTestStatisticsPerClust(genotypes, classLabels, clusterMemberIndices);
    
    result.nVariantsPosClass = this.countVariantsForGene(genotypes, true);
    result.nVariantsNegClass = this.countVariantsForGene(genotypes, false);    
    
    // this is the case when we have 0 examples of genotype "1" 
    // if ( Double.isNaN(result.observed) ) {  
    // }
    BitSet forShuffling = (BitSet) classLabels.clone();
    int pRightInt=0, pLeftInt=0;
    //double minHere = Double.MAX_VALUE, maxHere = -Double.MAX_VALUE;
    double sumHere = 0.0;
    
    Random r = new XorShiftRandom(randSeed);
    for ( int i=0; i<nIter; i++ ) {
      randomizeLabels(forShuffling, r);
      double oneRandStat; 
      if ( keepRandDistInMem ) {  // the per-cluster stats are calculated only if they will be kept in memory
        RandOutStats[] oneRandAllStats = calcTestStatisticsPerClust(genotypes, forShuffling, clusterMemberIndices); 
        result.randomizedAllStats[i] = oneRandAllStats; 
        oneRandStat = oneRandAllStats[0].testStat;
      } else {  // otherwise, only the overall statistics are calculated
        oneRandStat = calcTestStatistic(genotypes, forShuffling);
      }
      if ( oneRandStat >= result.observed ) pRightInt ++;
      if ( oneRandStat <= result.observed ) pLeftInt ++;
      //if ( oneRandStat > maxHere ) maxHere = randTestStat[i];
      //if ( oneRandStat > minHere ) minHere = randTestStat[i];
      sumHere += oneRandStat;
      randTestStat[i] = oneRandStat;
    }

    if ( keepRandDistInMem ) 
      result.randomizedStat = Arrays.copyOf(randTestStat, randTestStat.length);  // keep this in RAM? takes up a bit of memory
    else 
      result.randomizedStat = null;    
    
    Arrays.sort(randTestStat);
    result.randMin = randTestStat[0];
    result.randMax = randTestStat[ randTestStat.length-1 ];
    result.pLeft = (double) pLeftInt / nIter;
    result.pRight = (double) pRightInt / nIter;
    result.randAvg = (double) sumHere / nIter;
    
    // find percentiles
    for ( int i=0; i<=100; i++ ) {
      int index = Math.round( (float)i/100f * nIter );
      result.randPercentiles[i] = randTestStat[Math.min( index, randTestStat.length-1 )];
    }
    result.randMedian = result.randPercentiles[50];
    // make a histogram with 100 equal-width bins
    double binWidth = (result.randMax - result.randMin)/100.0;
    for ( int i=0; i<=100; i++ ) {
      result.randHistoBinsLowBound[i] = result.randMin + i*binWidth;
    }
    result.randHistoBinsLowBound[100] = result.randMax;   // they should already be nearly the same (save for rounding errors)
    int currBin = 0;
    for ( int i=0; i<nIter; i++) {  // the randTestStat is sorted here; use this to divide into bins elegantly
      // this could support the variable-length binning as well
      while ( currBin<99   // check if the value matches or exceeds the lower bound of the next bin
              && randTestStat[i] >= result.randHistoBinsLowBound[ currBin+1 ] ) {
        currBin++;
      }
      result.randHistoCounts[currBin]++;   // check is this works okay
    }

    return result;
  }

  
  
  
  
  
  
  
  
  /**
   * A randomization test for excess frequency of a genotype in a disease versus a control cohort; it is also used in
   * the ALFRED method where co-occurence of LOH and damaging variants is checked. Controls for population structure.
   * 
   * See reference: Park, Supek and Lehner (2018) Nature Communications 
   * "Systematic discovery of germline cancer predisposition genes through the identification of somatic second hits".
   * https://www.nature.com/articles/s41467-018-04900-7
   * 
   * @param args See description in code below.
   */
  public static void main(String[] args) throws IOException {
    
    if ( args.length <= 6 ) {
      System.out.println("A randomization test for frequency of a genotype in a case vs a control cohort. Can control for population stratification.");
      System.out.println("");
      System.out.println("Usage:\n");
      System.out.println("java -jar CancerGeneticsCode.jar $INFILE $LABEL $CLUSTERING $GENOTYPES $NUMITER $RANDSEED $MINVARIANTS $MAXT [$OUTFILESUFFIX] [$OUTFILE_RAND_LOG]\n");
      System.out.println("Output will be appended to the file named \"$LABEL__$CLUSTERING__geneRandomizerOutput.txt\" in the current working directory.");
      System.out.println("\n$INFILE is the file path to the input tab-separated table.");
      System.out.println("\n$LABEL is the name of the SINGLE column containing class labels, where \"1\" is cancer patient, \"0\" is healthy control, ");
      System.out.println("  IMPORTANT: Rows with anything except 0 or 1 are completely ignored (same as being absent from the table).");
      System.out.println("  Instead of the column name, you may also supply a 1-based column index here.");
      System.out.println("\n$CLUSTERING is the name of the SINGLE column with the cluster IDs for randomization (they may be any kind of String).");
      System.out.println("  Instead of the column name, you may supply a 1-based column index here.");
      System.out.println("  IMPORTANT: rows containing a \"?\" or \"NA\" for the cluster ID will be ignored!");
      System.out.println("\n$GENOTYPES is used to find the MULTIPLE genotype columns, which typically contain 0 or 1 ");
      System.out.println("  (for wt and mutated, respectively) but may contain any kind of floating point number. If this column contains a");
      System.out.println("  non-numeric value (such as \"?\" or \"NA\"), this person/row will be marked as missing value and skipped when");
      System.out.println("  randomizing. If you want this to match all the columns except the ones used above, just put \".+\"");
      System.out.println("  \nHere, you can supply either:");
      System.out.println("  -> a regular expression, such as \"^gene_.+\" to match all column that start with \"gene_\"; or just \".+\"");
      System.out.println("  -> a single number, then this column is looked up by its 1-based index. (script will randomize a single gene)");
      System.out.println("  -> a range of numbers, e.g. \"5-29\" (again 1-based index), then all genes in those columns will be looked up.");
      System.out.println("\n$NUMITER is the number of randomizations. For entire genome, recommended 2000000, which supports a genome-wide");
      System.out.println("  Bonferroni-corrected p value of 0.05. (comes from 1/0.05 * 20000 genes * factor 5 for enough precision = 2e6)");
      System.out.println("  For a ~1000 cancer genes, we need only 1e5 (100000) randomizations to get the same cancer_gene_set-wide p 0.05 cutoff.");
      System.out.println("\n$RANDSEED is the random seed; used for reproducibility. Default = 42." );
      System.out.println("\n$MINVARIANTS is min. number of variants (together in cancer+control) a gene must have, or it will not be tested at all." );
      System.out.println("  This is a poor man's version of the 'i-stat' filter, which skips test that can't give a result.");
      System.out.println("  Recommended: 3 (conservative).");
      // System.out.println("\n$MAXT Do you want to calculate the Westfall & Young maxT multiple-testing corrected p-values? Pass 'true' or 'false'.");
      // System.out.println("  Requires more memory; slows calculation slightly; causes the table rows not be output until all genes are done.");
      // System.out.println("  Default: false.");  // note: not thoroughly tested, currently.
      System.out.println("\n$OUTFILESUFFIX is the output file name suffix. It's okay if multiple jobs write to same file; the randomizer handles this well.");
      System.out.println("  This argument is optional; leave empty to just use the default file name \"args[1]__args[2]__geneRandomizerOutput.txt\"");
      System.out.println("\n$RANDOMIZE_VERBOSE is the output file name which will contain statistics (counts) for individual iterations of the randomization.");
      System.out.println("  This does not work well with multiple jobs writing to same file! This argument is optional; leave empty to not output.");
      System.out.println("  Pass \"stdout\" to write to standard output. Warning - this file could get very large for many iterations.");
      System.out.println("");
      System.out.println("\nExamples:");
      System.out.println("java AlfredRandomizer myData.txt isTumorPatient whichCluster ^nc_gene_.+ 10000000 42 3 _allGenes");
      System.out.println("java AlfredRandomizer myData.txt 2 whichCluster (TP53|BRCA1|BRCA2) 10000000 42 3");
      System.out.println("java AlfredRandomizer myData.txt isBRCAvsRest whichCluster 50-1050 2000000 42 4 _canGenes stdout");
      return;
    }
    
  
    AlfredRandomizer germ = new AlfredRandomizer(
            args[0],  // input file path (may be relative to current directory)
            "",       // patient ID column - it is currently ignored
            args[1],  // class labels column - by name (not regex!) or index
            args[2],  // cluster ID column - by name (not regex!) or index
            args[3]   // a regex matching gene columns, or single column index, or range
    );
    
    int nIter = Integer.parseInt( args[4] );
    int randSeed = Integer.parseInt( args[5] );

    // note that this "# of variants" can also be a fractional score, in case we decide to weigh variants somehow
    // these variants count across both the cancer AND the control ca
    double minNumVariants = Double.parseDouble(args[6]);
    
    // boolean calcMaxT = args[7].equalsIgnoreCase("true");  // this parameter is not operational
    boolean calcMaxT = false;
    String randLogFile = args.length>=9 ? args[8] : null;
    
    String outputFilenameShort = null, outputFilename = null;
    //outputFilename = args[1]+"__"+args[2]+"__geneRandomizerOutput.txt";
    outputFilenameShort = germ.classColumnName+"__"+germ.clusterColumnName; //+"__geneRandomizerOutput.txt";
    outputFilenameShort = outputFilenameShort.replace('\\', '_').replace('/', '_').replace('*', '_');
    
    System.out.printf("[" + outputFilenameShort + "] ***** Starting randomization with user parameters: %s\n", Arrays.toString(args));
    outputFilename = "geneRandOut__" + outputFilenameShort + ( args.length>=8 ? args[7] : "" )  +  ".txt";
    
    //System.out.println( "genotype\t" + GeneRandomizationResult.getTabSepHeader() );
    
    // Main loop is gene by gene. This could, in principle, be parallelized, if we see it's too slow.

    List<GeneRandomizationResult> geneRandResults = new ArrayList<>();
    
    int nTestedHyp = 0;
    
    for ( int i=0; i<germ.numGenes; i++ ) {
      
      String geneName = germ.genotypeNames.get(i);
      
      System.out.print("[" + outputFilenameShort + "] ");
      // here check the minNumVariants
      double nPos = germ.countVariantsForGene( germ.genotypes[i] , true);
      double nNeg = germ.countVariantsForGene( germ.genotypes[i] , false);
      System.out.printf( "Gene %s: %.1f variants in positive and %.1f in negative examples; total %.1f. ",
              geneName, nPos, nNeg, nPos+nNeg);
              
      if ( nPos+nNeg < minNumVariants ) {
        System.out.printf(" This is insufficient for statistical testing (required total variants: %.0f). Skipping gene. \n", minNumVariants );
        continue;
      }
      long millis = System.currentTimeMillis();
      GeneRandomizationResult aGeneResult = germ.randomizeOneGene(germ.genotypes[i],
              germ.classLabels, 
              nIter, 
              randSeed,   // all genes are randomized with the same seed (effectively doing the same randomizations of the healthy/diseased people) 
              calcMaxT || randLogFile!=null     // if calculating the Westfall-Young maxT ~OR~ if outputting randomization log, will need to keep the radnomized instances in memory
      );
      aGeneResult.genotypeName = germ.genotypeNames.get(i);
      geneRandResults.add(aGeneResult);   // save this for later. Could be large in memory if we selected calcMaxT; otherwise very compact.
      
      // aGeneResult.observedAllStats // for individual clusters
      
      nTestedHyp++;
      System.out.printf( "Testing for gene %s done in %.1f seconds for %d iterations.\n", geneName,
              (System.currentTimeMillis()-millis)/1000.0, nIter  );
    }

    
    
    
    
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // <editor-fold defaultstate="collapsed" desc=" (OPTIONAL) asked to output the data on individual randomizations? ">
    
    if ( randLogFile!=null  ) {
      
      StringBuilder sb = new StringBuilder();
      Formatter rows;
      
      if ( randLogFile.equalsIgnoreCase("stdout") ) {
        rows = new Formatter(System.out, Locale.US);
        rows.format("\n");
      } else {
        rows = new Formatter(randLogFile, "ISO-8859-1", Locale.US);
      }
      
      // ordering of the per-cluster columns -- output the clusters in a nice/sorted order (ascending/alphabetical) 
      TreeMap<String,Integer> sortedClustNames = new TreeMap<>();
      for ( int i=0; i<germ.clusterMemberIndices.length; i++ ) {
        sortedClustNames.put(  germ.clusterNamesList.get(i), i );
      }
      int[] sortedClIndices = new int[ 1 + sortedClustNames.size() ];
      sortedClIndices[0] = 0;  // clIdx==0 is the data for all clusters, as stored in RandOutStats.observedAllStats[]
      int k=1;
      for ( Integer clustId : sortedClustNames.values() ) {
        sortedClIndices[k] = clustId+1; // adding +1 because the clusters are stored in RandOutStats.observedAllStats[] by 1-based index (0 is reseved for the summary across clusters)
        k++;
      }      
      
      // rows.format("genotype\titeration\tallCl_testStat\tallCl_sumPos\tallCl_numPos\tallCl_sumNeg\tallCl_numNeg");
      /*
      for ( int i=0; i<germ.clusterMemberIndices.length; i++ ) {
        rows.format( "\tclName=%s_testStat\tclName=%s_sumPos\tclName=%s_numPos\tclName=%s_sumNeg\tclName=%s_numNeg", 
                //germ.clusterNamesList.get(i), germ.clusterNamesList.get(i), germ.clusterNamesList.get(i), germ.clusterNamesList.get(i), germ.clusterNamesList.get(i)
                // i+1,i+1,i+1,i+1,i+1
        );
      }*/
      rows.format("genotype\titeration\t");
      for ( int clIdx : sortedClIndices ) {
        if ( clIdx==0 ) { // should come first in this array
          rows.format("allCl_testStat\tallCl_sumPos\tallCl_numPos\tallCl_sumNeg\tallCl_numNeg");          
        } else {
          String aClName = germ.clusterNamesList.get( clIdx-1 );   // ... because they're stored 0-based in the NamesList
          rows.format("\tclName=%s_testStat\tclName=%s_sumPos\tclName=%s_numPos\tclName=%s_sumNeg\tclName=%s_numNeg", 
                  aClName, aClName, aClName, aClName, aClName
                  );
        }
      }
      rows.format("\n");
      rows.flush();
      
      for ( int i=0; i<geneRandResults.size(); i++ ) {  //  this can be smaller than germ.numGenes, if not all genes were tested
        
        GeneRandomizationResult aGeneResult = geneRandResults.get(i);
        String geneName = aGeneResult.genotypeName;
        
        rows.format("%s\t%d", geneName, -1);
        // for ( int clIdx = 0; clIdx < germ.clusterMemberIndices.length + 1; clIdx++ ) {  
        for ( int clIdx : sortedClIndices ) {  
          // clIdx==0 is the data for all clusters;
          RandOutStats ros = aGeneResult.observedAllStats[clIdx];
          rows.format("\t%.5f\t%.2f\t%d\t%.2f\t%d", ros.testStat, ros.sumPos, ros.numPos, ros.sumNeg, ros.numNeg );
        }
        rows.format("\n"); rows.flush();
        
        
        for ( int j=0; j < aGeneResult.randomizedAllStats.length; j++ ) {
          rows.format( "%s\t%d", geneName, j+1 );
          // for ( int clIdx = 0; clIdx < germ.clusterMemberIndices.length + 1; clIdx++ ) {  
          for ( int clIdx : sortedClIndices ) {  
            // clIdx==0 is the data for all clusters;
            RandOutStats ros = aGeneResult.randomizedAllStats[j][clIdx];
            rows.format("\t%.5f\t%.2f\t%d\t%.2f\t%d", ros.testStat, ros.sumPos, ros.numPos, ros.sumNeg, ros.numNeg );
          }
          rows.format("\n"); rows.flush();
        }
        
      } // gene by gene 

      rows.close();  // in case writing to file
    }
    // </editor-fold>  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    
    
    
    
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // <editor-fold defaultstate="collapsed" desc=" (OPTIONAL) following part is ran only if specifying 'calcMaxT' ">
    
    // now, optionally find the maxT-corrected p-values for multiple testing (Westfall and Young 1993).
    if ( calcMaxT ) {
      double[] genePvalsRight = new double[nTestedHyp];  // these arrays contain only the hypotheses we actually tested
      double[] geneTestStats = new double[nTestedHyp];
      int curHyp = 0;
      // for some genes we performed no test (if not enough variants there) 
      // ... so the #_geneRandResults is smaller than #_genes, and the indices may well not match!
      
      for ( int i=0; i<geneRandResults.size(); i++ ) {  
        //genePvalsRight[i] = aGeneResult.pRight;
        //geneTestStats[i] = aGeneResult.observed;
        geneTestStats[curHyp] = geneRandResults.get(i).observed;
        genePvalsRight[curHyp] = geneRandResults.get(i).pRight;
        curHyp++; 
      }
      // (1) find ordering of genes by the test statistic
      int[] sortedIdxAsc = WekaUtilsSortOnly.sort(geneTestStats);  // sort smallest to largest
      int[] sortedIdxDesc = new int[nTestedHyp];
      for ( int i=0; i<nTestedHyp; i++ ) sortedIdxDesc[i] = sortedIdxAsc[ sortedIdxAsc.length - 1 - i ] ;
      
      List<GeneRandomizationResult> geneRandResultsSorted = new ArrayList<>();
      for ( int i=0; i<nTestedHyp; i++ ) {
        geneRandResultsSorted.add( geneRandResults.get( sortedIdxDesc[i] ) );  
      }
      // (2) given the order above, calculate the successive maxima of test statistics in every randomization
      // http://white.stanford.edu/~knk/Psych216A/FinalProjects/Moqian_Wenying/Psych216_multi_test_project_v4.pdf
      double[][] u = new double[nIter][nTestedHyp];

      for ( int b=0; b<nIter; b++ ) {  // ... iter by iter (across all genes)
        for ( int j=geneRandResultsSorted.size()-1; j>=0; j-- ) {  // ... gene by gene - in sorted order - backwards 
                                                                   // (Lowest test statistic first!)
          if ( j==geneRandResultsSorted.size()-1 ) {
            u[b][j] = geneRandResultsSorted.get(j).randomizedStat[b]; 
          } else {
            double here = geneRandResultsSorted.get(j).randomizedStat[b];
            double prev = u[b][j+1];  // geneRandResultsSorted.get(j+1).randomizedStat[b];
            u[b][j] = Math.max(here, prev); 
            //u[b][j] = geneRandResultsSorted.get(j).randomizedStat[b];
          }
        }
      }
      // (3) estimate the adjusted p-values
      double[] geneMaxT_corrPval = new double[geneRandResultsSorted.size()];  // this array is also going to be ordered as the sorted genes/hypotheses
      for ( int j=0; j<nTestedHyp; j++ ) {  // loop over genes(=hypotheses) ... in u[][], the genes are sorted by test statistic
        int counterForGene=0;
        for ( int b=0; b<nIter; b++ ) {
          if ( u[b][j] >= geneRandResultsSorted.get(j).observed )  // ... we check how often is higher than real
            counterForGene++;
        }  
        geneMaxT_corrPval[j] = (double) counterForGene / nIter;   
      }
      // (4) enforce the monotonicity constraint
      double[] geneMaxT_corrPvalMono = new double[ geneMaxT_corrPval.length ];
      for ( int i=0; i<geneRandResultsSorted.size(); i++ ) {  // ... gene by gene - in sorted order
        if (i==0) {
          geneMaxT_corrPvalMono[i] = geneMaxT_corrPval[i];
        } else {
          geneMaxT_corrPvalMono[i] = Math.max( geneMaxT_corrPval[i], geneMaxT_corrPvalMono[i-1] ); 
        }
      }
      // (5) write out the corrected p-values into geneRandResults[] objects
      for ( int i=0; i<nTestedHyp; i++ ) {  // ... gene by gene - here not sorted
        geneRandResultsSorted.get(i).pRightMaxT = geneMaxT_corrPvalMono[i];
      }
    }  // ....calc MaxT or not?
    
    // </editor-fold>  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    
    
    
    
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ finally, output the randomization results summary to a file ~~~~~~~~~~~~~~~~~~~~~~~
    
    List<String> rowsToOutput = new ArrayList<>();    
    for ( int i=0; i<nTestedHyp; i++ ) {
      rowsToOutput.add( geneRandResults.get(i).genotypeName + "\t" + geneRandResults.get(i).toString() + "\n" );
      // System.out.println( geneRandResults.get(i).genotypeName + "\t" + geneRandResults.get(i).toString() );
    }
    // try to write to a file, taking care not to overlap with another process
    // String outputFilename = args[0] + "_randTest.txt";
    
    boolean skipHeader = false;
    if ( new File( outputFilename ).exists()  ) {  // don't write the header multiple times... annoying
      skipHeader = true;
    }
    
    RandomAccessFile stream = null;
    while (stream==null) {
      try {
        stream = new RandomAccessFile(outputFilename, "rw");   // if file not writeable, maybe something else is locking it (in Windows)
      } catch (FileNotFoundException fileNotFoundException) {
        try {
          Thread.sleep(1000);  // wait a bit until trying to lock file again
        } catch (InterruptedException ie) {
          Thread.currentThread().interrupt();
        }  
      }
    }
    stream.seek( stream.length() );
    FileChannel channel = stream.getChannel();
    FileLock lock = null;
    while ( lock == null ) {  // this loop can go forever - if the other process that has the lock doesn't release it... 
      try {
        lock = channel.tryLock();
      } catch (final OverlappingFileLockException e) {
        try {
          Thread.sleep(1000);  // wait a bit until trying to lock file again
        } catch (InterruptedException ie) {
          Thread.currentThread().interrupt();
        }
      }
    }
    if ( ! skipHeader ) {
      String aRow = "genotype\t" + GeneRandomizationResult.getTabSepHeader() + "\n";
      channel.write( ByteBuffer.wrap(aRow.getBytes()) );
      // stream.writeChars( "genotype\t" + GeneRandomizationResult.getTabSepHeader() + "\n" );
    }
    for ( String aRow : rowsToOutput ) {
      // byte[] strBytes = aRow.getBytes();
      //ByteBuffer buffer = ByteBuffer.allocate(strBytes.length);
      //buffer.put(strBytes);
      // buffer.flip();
      channel.write( ByteBuffer.wrap(aRow.getBytes()) );
      // stream.writeChars( aRow ); // doesn't work 
    }
    //channel.force(false);
    lock.release();  // this can never be null because of the loop above
    stream.close();
    channel.close();    
   
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~ end output the randomization results summary to a file ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    
    
  }
  
  
  
  
  
  
  
  
  
  // ===================================================================================================================
  // ======================= supporting code: sorting, random numbers, Fisher's exact test ============================
  
  
  /**
   * A fast random number generator. Not thread safe!
   * Code from http://www.javamex.com/tutorials/random_numbers/java_util_random_subclassing.shtml )
   */
  public static class XorShiftRandom extends Random {
    private long seed;
    public XorShiftRandom() {
      seed = System.nanoTime();
    }
    public XorShiftRandom(long seed) {
      this.seed = seed;
    }
    /** Not thread-safe. */
    @Override
    protected int next(int nbits) {
      long x = this.seed;
      x ^= (x << 21);
      x ^= (x >>> 35);
      x ^= (x << 4);
      this.seed = x;
      x &= ((1L << nbits) -1);
      return (int) x;
    }

  }  // end class XorShiftRandom
  
  
  /**
   * Performs the Fisher Exact Test to determine if there is an association between two categorical variables.
   * 
   * Adapted from JavaScript code by Oyvind Langsrud, see
   * <a href="http://www.langsrud.com/fisher.htm"> http://www.langsrud.com/fisher.htm </a>
   * 
   * @author Fran Supek (fran.supek[AT]irb.hr)
   */
  public static class FisherExactTest {

    /**
     * Reference: "Lanczos, C. 'A precision approximation of the gamma function', J. SIAM Numer. Anal., B, 1, 86-96, 1964."
     * Translation of  Alan Miller's FORTRAN-implementation. See http://lib.stat.cmu.edu/apstat/245
     */
    private double lnfact(int n) {

      if (n<=1) return(0.0);

      int z = n+1;
      //return(lngamm(n+1));

      double x = 0;
      x += 0.1659470187408462e-06/(z+7);
      x += 0.9934937113930748e-05/(z+6);
      x -= 0.1385710331296526    /(z+5);
      x += 12.50734324009056     /(z+4);
      x -= 176.6150291498386     /(z+3);
      x += 771.3234287757674     /(z+2);
      x -= 1259.139216722289     /(z+1);
      x += 676.5203681218835     /(z);
      x += 0.9999999999995183;
      return(Math.log(x)-5.58106146679532777-z+(z-0.5)*Math.log(z+6.5));

    }

    //private double lnbico(double n, double k) {
    private double lnbico(int n, int k) {
      return( lnfact(n)-lnfact(k)-lnfact(n-k) );
    }


    private int sn11,sn1_,sn_1,sn;
    //private double sn11,sn1_,sn_1,sn;
    private double sprob;


    private double hyper0(int n11i, int n1_i, int n_1i, int ni) {
      if ( ! ( n1_i != 0 || n_1i != 0 || ni != 0 ) ) {
        if ( ! (n11i % 10 == 0) ) {
          if ( n11i==sn11+1 ) {
            sprob *= ( (double) (sn1_-sn11)/(n11i) )*( (double) (sn_1-sn11)/(n11i+sn-sn1_-sn_1) );
            sn11 = n11i;
            return sprob;
          }
          if (n11i==sn11-1 ) {
            sprob *= ( (double) (sn11)/(sn1_-n11i) )*( (double) (sn11+sn-sn1_-sn_1)/(sn_1-n11i) );
            sn11 = n11i;
            return sprob;
          }
        }
        sn11 = n11i;
      } else {
        sn11 = n11i;
        sn1_ = n1_i;
        sn_1 = n_1i;
        sn   = ni;
      }

      sprob = Math.exp(lnbico(sn1_,sn11)+lnbico(sn-sn1_,sn_1-sn11)-lnbico(sn,sn_1));

      return sprob;
    }


    public double left, right, twotail;

    final public double exact22(int n11_, int n12_, int n21_, int n22_) {

      if ( n11_<0 || n12_<0 || n21_<0 || n22_<0 )
        throw new IllegalArgumentException("Counts passed to exact22() must all" +
                " be nonnegative.");

      int n1_ = n11_ + n12_;
      int n_1 = n11_ + n21_;
      int n   = n11_ + n12_ + n21_ + n22_;

      //double prob = exact(n11_, n1_, n_1, n);
      double sleft, sright, sless, slarg;
      double prob;

      int i, j;

      int max = n1_;
      if (n_1<max) max=n_1;

      int min = n1_ + n_1 - n;
      if (min < 0) min=0;
      if (min == max) {

        sless = 1; sright= 1;
        sleft = 1; slarg = 1;
        prob = 1;

      } else {

        prob = hyper0(n11_,n1_,n_1,n);
        double probReduced = 0.99999999*prob;

        sleft = 0;
        double p = hyper0(min,0,0,0);
        for (i = min+1; p < probReduced; i++) {
          sleft += p;
          p = hyper0(i,0,0,0);
        }
        i--;
        if (p < 1.00000001*prob) sleft += p;
          else i--;

        sright = 0;
        p = hyper0(max,0,0,0);
        for (j = max-1; p < probReduced; j--) {
          sright += p;
          p = hyper0(j,0,0,0);
        }
        j++;
        if (p<1.00000001*prob) sright += p;
          else j++;

        if ( Math.abs(i-n11_) < Math.abs(j-n11_) ) {
          sless = sleft;
          slarg = 1 - sleft + prob;
        } else {
          sless = 1 - sright + prob;
          slarg = sright;
        }

      }

      this.left    = sless;
      this.right   = slarg;
      this.twotail = sleft + sright;
      if (twotail>1) twotail=1;
      return prob;
    }


  }
  
  
  
  
  /**
   * Quicksort code, copied from weka.core.Utils of Weka 3-7-10. 
   * 
   * @author Fran Supek, CRG.
   */
  public static class WekaUtilsSortOnly {

    /**
     * Tests if the given value codes "missing".
     *
     * @param val the value to be tested
     * @return true if val codes "missing"
     */
    public static boolean isMissingValue(double val) {

      return Double.isNaN(val);
    }

    /**
     * Returns the value used to code a missing value.  Note that
     * equality tests on this value will always return false, so use
     * isMissingValue(double val) for testing..
     *
     * @return the value used as missing value.
     */
    public static double missingValue() {

      return Double.NaN;
    }


    /**
     * Returns index of minimum element in a given
     * array of doubles. First minimum is returned.
     *
     * @param doubles the array of doubles
     * @return the index of the minimum element
     */
    public static /*@pure@*/ int minIndex(double[] doubles) {

      double minimum = 0;
      int minIndex = 0;

      for (int i = 0; i < doubles.length; i++) {
        if ((i == 0) || (doubles[i] < minimum)) {
      minIndex = i;
      minimum = doubles[i];
        }
      }

      return minIndex;
    }

    /**
     * Replaces all "missing values" in the given array of double values with
     * MAX_VALUE.
     *
     * @param array the array to be modified.
     */
    public static void replaceMissingWithMAX_VALUE(double[] array) {

      for (int i = 0; i < array.length; i++) {
        if (isMissingValue(array[i])) {
          array[i] = Double.MAX_VALUE;
        }
      }
    }


    /**
     * Sorts a given array of integers in ascending order and returns an 
     * array of integers with the positions of the elements of the original 
     * array in the sorted array. The sort is stable. (Equal elements remain
     * in their original order.)
     *
     * @param array this array is not changed by the method!
     * @return an array of integers with the positions in the sorted
     * array.
     */
    public static /*@pure@*/ int[] sort(int[] array) {

      int[] index = initialIndex(array.length);
      int[] newIndex = new int[array.length];
      int[] helpIndex;
      int numEqual;

      quickSort(array, index, 0, array.length - 1);

      // Make sort stable
      int i = 0;
      while (i < index.length) {
        numEqual = 1;
        for (int j = i + 1; ((j < index.length)
                 && (array[index[i]] == array[index[j]]));
         j++) {
      numEqual++;
        }
        if (numEqual > 1) {
      helpIndex = new int[numEqual];
      for (int j = 0; j < numEqual; j++) {
        helpIndex[j] = i + j;
      }
      quickSort(index, helpIndex, 0, numEqual - 1);
      for (int j = 0; j < numEqual; j++) {
        newIndex[i + j] = index[helpIndex[j]];
      }
      i += numEqual;
        } else {
      newIndex[i] = index[i];
      i++;
        }
      }
      return newIndex;
    }

    /**
     * Sorts a given array of doubles in ascending order and returns an
     * array of integers with the positions of the elements of the
     * original array in the sorted array. NOTE THESE CHANGES: the sort
     * is no longer stable and it doesn't use safe floating-point
     * comparisons anymore. Occurrences of Double.NaN are treated as 
     * Double.MAX_VALUE.
     *
     * @param array this array is not changed by the method!
     * @return an array of integers with the positions in the sorted
     * array.  
     */
    public static /*@pure@*/ int[] sort(/*@non_null@*/ double[] array) {

      int[] index = initialIndex(array.length);
      if (array.length > 1) {
        array = (double[])array.clone();
        replaceMissingWithMAX_VALUE(array);
        quickSort(array, index, 0, array.length - 1);
      }
      return index;
    }


    /**
     * Initial index, filled with values from 0 to size - 1.
     */
    private static int[] initialIndex(int size) {

      int[] index = new int[size];
      for (int i = 0; i < size; i++) {
        index[i] = i;
      }
      return index;
    }

    /**
     * Sorts left, right, and center elements only, returns resulting center as pivot.
     */
    private static int sortLeftRightAndCenter(double[] array, int[] index, int l, int r) {

      int c = (l + r) / 2;
      conditionalSwap(array, index, l, c);
      conditionalSwap(array, index, l, r);
      conditionalSwap(array, index, c, r);
      return c;
    }

    /**
     * Swaps two elements in the given integer array.
     */
    private static void swap(int[] index, int l, int r) {

      int help = index[l];
      index[l] = index[r];
      index[r] = help;
    }

    /**
     * Conditional swap for quick sort.
     */
    private static void conditionalSwap(double[] array, int[] index, int left, int right) {

      if (array[index[left]] > array[index[right]]) {
        int help = index[left];
        index[left] = index[right];
        index[right] = help;
      }
    }

    /**
     * Partitions the instances around a pivot. Used by quicksort and kthSmallestValue.
     *
     * @param array the array of doubles to be sorted
     * @param index the index into the array of doubles
     * @param l the first index of the subset 
     * @param r the last index of the subset 
     *
     * @return the index of the middle element
     */
    private static int partition(double[] array, int[] index, int l, int r,
                                 double pivot) {

      r--;
      while (true) {
        while ((array[index[++l]] < pivot));
        while ((array[index[--r]] > pivot));
        if (l >= r) {
          return l;
        }
        swap(index, l, r);
      }
    }

    /**
     * Partitions the instances around a pivot. Used by quicksort and kthSmallestValue.
     *
     * @param array the array of integers to be sorted
     * @param index the index into the array of integers
     * @param l the first index of the subset 
     * @param r the last index of the subset 
     *
     * @return the index of the middle element
     */
    private static int partition(int[] array, int[] index, int l, int r) {

      double pivot = array[index[(l + r) / 2]];
      int help;

      while (l < r) {
        while ((array[index[l]] < pivot) && (l < r)) {
          l++;
        }
        while ((array[index[r]] > pivot) && (l < r)) {
          r--;
        }
        if (l < r) {
          help = index[l];
          index[l] = index[r];
          index[r] = help;
          l++;
          r--;
        }
      }
      if ((l == r) && (array[index[r]] > pivot)) {
        r--;
      } 

      return r;
    }

    /**
     * Implements quicksort with median-of-three method and explicit sort for
     * problems of size three or less. 
     *
     * @param array the array of doubles to be sorted
     * @param index the index into the array of doubles
     * @param left the first index of the subset to be sorted
     * @param right the last index of the subset to be sorted
     */
    //@ requires 0 <= first && first <= right && right < array.length;
    //@ requires (\forall int i; 0 <= i && i < index.length; 0 <= index[i] && index[i] < array.length);
    //@ requires array != index;
    //  assignable index;
    private static void quickSort(/*@non_null@*/ double[] array, /*@non_null@*/ int[] index, 
                                 int left, int right) {

      int diff = right - left;

      switch (diff) {
      case 0 :

        // No need to do anything
        return;
      case 1 :

        // Swap two elements if necessary
        conditionalSwap(array, index, left, right);
        return;
      case 2 :

        // Just need to sort three elements
        conditionalSwap(array, index, left, left + 1);
        conditionalSwap(array, index, left, right);
        conditionalSwap(array, index, left + 1, right);
        return;
      default :

        // Establish pivot
        int pivotLocation = sortLeftRightAndCenter(array, index, left, right);

        // Move pivot to the right, partition, and restore pivot
        swap(index, pivotLocation, right - 1); 
        int center = partition(array, index, left, right, array[index[right - 1]]);
        swap(index, center, right - 1);

        // Sort recursively
        quickSort(array, index, left, center - 1);
        quickSort(array, index, center + 1, right);
      }
    }

    /**
     * Implements quicksort according to Manber's "Introduction to
     * Algorithms".
     *
     * @param array the array of integers to be sorted
     * @param index the index into the array of integers
     * @param left the first index of the subset to be sorted
     * @param right the last index of the subset to be sorted
     */
    //@ requires 0 <= first && first <= right && right < array.length;
    //@ requires (\forall int i; 0 <= i && i < index.length; 0 <= index[i] && index[i] < array.length);
    //@ requires array != index;
    //  assignable index;
    private static void quickSort(/*@non_null@*/ int[] array, /*@non_null@*/  int[] index, 
                                  int left, int right) {

      if (left < right) {
        int middle = partition(array, index, left, right);
        quickSort(array, index, left, middle);
        quickSort(array, index, middle + 1, right);
      }
    }



  }
  
  // ===================================================================================================================
  

  
  
}

