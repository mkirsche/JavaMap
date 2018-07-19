import java.util.*;
import java.io.*;
public class alignment {
	// k is the length of kmers used for sketches
	static int k = 15;
	// w is the window size used for sketches
	static int w = 5;
	// map maps base pairs to integers ('A' = 1, 'C' = 2, 'G' = 3, 'T' = 4)
	static int[] map;
	// EPS is the maximum relative distance between shared kmers that can form a chain
	static int EPS = 500;
	// f is the proportion of kmers which are considered repeats
	static double f = 0.001;
	// batchSize is how many reads to put in the database at once
	static int batchSize = 21000;
	// minOverlap is the minimum number of base pairs which an overlap needs to span
	static int minOverlap = 100;
	// c is the minimum number of matching minimizers in an overlap
	static int c = 4;
	// posBits is the number of bits needed to store the position within a read
	static int posBits = 20;
	// pos1Mask will have 1s in positions which are used to encode query position in hit (more below)
	static long pos1Mask;
	// pos2Mask will have 1s in positions which are used to encode target position in hit (more below)
	static long pos2Mask;
	// writeFile is whether or not to write the alignments to a file
	static boolean writeFile = true;
	// overallFilter is whether or not to use a single repeat filter for all batches
	static boolean overallFilter = false;
	// Method for initializing global variables
	static void init()
	{
		map = new int[256];
		map['A'] = 0; map['C'] = 1; map['G'] = 2; map['T'] = 3;
		pos1Mask = ((1L<<posBits) - 1) << posBits;
		pos2Mask = (1<<posBits) - 1;
	}
public static void main(String[] args) throws IOException
{
	// Variables used for timing metrics
	long minimizerTime = 0;
	long chainTime = 0;
	long startTime = System.currentTimeMillis();
	PrintWriter out = new PrintWriter(System.out);
	if(writeFile) out = new PrintWriter(new File("/home/mkirsche/ecoli/oxford.paf"));
	
	// Initialize global variables
	init();
	
	String fn = "/home/mkirsche/ecoli/oxford.fasta";
	Scanner input = new Scanner(new FileInputStream(new File(fn)));
	
	// Get reads/names from fasta file
	ArrayList<String> names = new ArrayList<String>();
	ArrayList<String> reads = new ArrayList<String>();
	while(input.hasNext())
	{
		String name = input.next().substring(1);
		input.nextLine();
		names.add(name);
		String read = input.nextLine();
		reads.add(read);
	}
	int readsProcessed = 0; // Number of reads added to minimizer database so far
	int n = reads.size(); // Total number of reads
	
	// A list of hash values which show up too frequently and are considered repeats
	HashSet<Long> toRemove = new HashSet<Long>();
	
	// If overallFilter is turned on, an initial pass over all reads is made 
	// to find the frequency threshold for repeated minimizers
	if(overallFilter)
	{
		System.out.println("Counting minimizer frequencies");
		
		// A map of hash value to frequency
		HashMap<Long, Integer> overallFreq = new HashMap<Long, Integer>();
		for(int i = 0; i<n; i++)
		{
			String read = reads.get(i);
			long curTime = System.currentTimeMillis();
			long[] cur = getMinimizers(read);
			minimizerTime += System.currentTimeMillis() - curTime;
			
			// Increase the frequency of all minimizers in this read by 1
			for(long x : cur)
			{
				int keepBits = posBits + 1;
				long mVal = (x >> keepBits);
				
				Integer oldFreq = overallFreq.get(mVal);
				overallFreq.put(mVal, 1 + (oldFreq == null ? 0 : oldFreq));
			}
		}
		System.out.println("Distinct minimizer count: " + overallFreq.size());
		
		// Get a sorted list of all frequencies
		ArrayList<Integer> freqs = new ArrayList<Integer>();
		for(long x : overallFreq.keySet())
		{
			freqs.add(overallFreq.get(x));
		}
		Collections.sort(freqs);
		
		// Get the cutoff based on the upper proportion of the sorted list
		int repeatedMinimizers = (int)(f * freqs.size());
		int repeatThreshold = freqs.get(freqs.size() - 1 - repeatedMinimizers);
		
		// Add all repeat minimizers to a set
		for(long x : overallFreq.keySet())
		{
			if(overallFreq.get(x) > repeatThreshold)
				toRemove.add(x);
		}
		
		// No longer need frequencies since repeats are stored in a set
		overallFreq = null; 
		System.out.println("Keep minimizers with freq <= " + repeatThreshold);
	}

	// Process the reads in batches to avoid the database getting too large  
	int numBatches = (n + batchSize - 1) / batchSize;
	for(int iter = 0; iter<n; iter += batchSize)
	{
		int batch = (iter/batchSize + 1);
		System.out.println("Processing batch " + batch + " of " + numBatches);
		
		// Initialize database
		MinimizerDatabase mdb = new MinimizerDatabase();
		System.out.println("Computing minimizer database");
		
		// Iterate over reads in the batch and add their minimizers
		for(int j = iter; j < n && j < iter + batchSize; j++)
		{
			String read = reads.get(j);
			long curTime = System.currentTimeMillis();
			long[] cur = getMinimizers(read);
			minimizerTime += System.currentTimeMillis() - curTime;
			for(long x : cur)
			{
				// The number of bits which represent position/strand
				int keepBits = posBits + 1;
				
				// The bits which represent position/strand
				long psMask = ((1L << keepBits) - 1);
				
				// Representation of sequence/position/strand
				long sps = (x & psMask) + (((long)readsProcessed) << keepBits);
				
				long minimizerValue = x >> keepBits;
				if(overallFilter && toRemove.contains(minimizerValue)) continue;
				mdb.put(minimizerValue, sps);
			}
			readsProcessed++;
		}
		
		// Sort and add map to database for fast queries
		System.out.println("Finalizing minimizer database");
		mdb.finish();
		
		// If the repeat filtering is per-batch, perform that
		int repThreshold = -1;
		if(!overallFilter)
		{
			System.out.println("Filtering repeat minimizers");
			ArrayList<Integer> freqs = new ArrayList<Integer>();
			int tot = 1;
			for(int i = 1; i<mdb.size; i++)
			{
				if(mdb.hashes[i] != mdb.hashes[i-1])
				{
					freqs.add(tot);
					tot = 1;
				}
				else tot++;
			}
			freqs.add(tot);
			Collections.sort(freqs);
			int repeatedMinimizers = (int)(f * freqs.size());
			repThreshold = freqs.get(freqs.size() - 1 - repeatedMinimizers);
			System.out.println("Keep minimizers with freq <= " + repThreshold);
		}
		
		// Query all reads against the datbase
		System.out.println("Querying reads against database");
		for(int i = 0; i<n; i++)
		{
			if(i > 0 && i%5000 == 0)
			{
				System.out.println("Processed " + i + " queries");
			}
			
			String read = reads.get(i);
			
			// Get the minimizers for the query
			long curTime = System.currentTimeMillis();
			long[] cur = getMinimizers(read);
			minimizerTime += System.currentTimeMillis() - curTime;
			
			// Construct a list of hits - minimizer matches to database
			LongList hs = new LongList();
			for(long x : cur)
			{
				int keepBits = posBits + 1;
				int xstrand = (int)(x&1);
				int xpos = (int)((x & ((1L << keepBits) - 1)) >> 1);
				long minimizerValue = x >> keepBits;
				
				// Query the database for the minimizerValue
				Integer mdbPos = mdb.startRange.get(minimizerValue);
				
				// No hits found - ignore this minimizer
				if(mdbPos == null) continue;
				
				// If this minimizer occurs too much, ignore it
				if(!overallFilter && mdbPos + repThreshold < mdb.size 
					&& mdb.hashes[mdbPos + repThreshold] == minimizerValue)
				{
					continue;
				}
				
				// Create a hit for each occurrence of the minimizer
				while(mdbPos < mdb.size && mdb.hashes[mdbPos] == minimizerValue)
				{
					// Database entry is (seqId, targetPos, strand)
					long y = mdb.poss[mdbPos];
					mdbPos++;

					// Encode hit as (seqId, strand, queryPos, targetPos)
					int ystrand = (int)(y&1);
					int ypos = (int)((y & ((1L << keepBits) - 1)) >> 1);
					int seqId = (int)(y >> keepBits);
					if(seqId == i) continue;
					int nStrand = xstrand ^ ystrand;
					long h = getHit(seqId, nStrand, xpos, ypos);
					hs.add(h);
				}
			}
			
			// Sort hits based on (seq, strand, position difference)
			hs.sort(new Comparator<Long>() {
				public int compare(Long a, Long b)
				{
					return compareHits(a, b);
				}
			});
			
			// Cluster hits based on same seq/strand and
			// having position differences within EPS
			int b = 0;
			long[] ha = hs.a;
			for(int e = 0; e<hs.size; e++)
			{
				if(e == hs.size - 1  || getHitSeq(ha[e]) != getHitSeq(ha[e+1]) 
					|| getHitStrand(ha[e]) != getHitStrand(ha[e+1]) 
					|| getHitDiff(ha[e+1]) - getHitDiff(ha[e]) > EPS)
				{
					if(e >= b + c - 1)
					{
						curTime = System.currentTimeMillis();
						int strand = getHitStrand(hs.a[e]);
						int[] lis = alignmentChain(hs, b, e, strand > 0);
						chainTime += System.currentTimeMillis() - curTime;
						if(lis[0] != -1)
						{
							// Get id of target sequence from hit
							int j = getHitSeq(hs.a[e]);
							
							// Print query name/length/start/end
							out.printf("%s\t%s\t%s\t%s\t", 
									names.get(i),
									reads.get(i).length(),
									lis[0],
									lis[1]);
							
							// Print strand and target name/length/start/end
							out.printf("%s\t%s\t%s\t%s\t%s\t", 
									strand == 0 ? '+' : '-',
									names.get(j),
									reads.get(j).length(),
									lis[2],
									lis[3]);
							
							// Print kmer match count, k, and mapping quality
							out.printf("%s\t%s\t%s\n", 
									lis[4],
									k,
									255);
							
						}
					}
					b = e+1;
				}
			}
		}
	}
	long endTime = System.currentTimeMillis();
	System.out.println("Total time (ms): " + (endTime - startTime));
	System.out.println("Time computing minimizers (ms): " + minimizerTime);
	System.out.println("Time chaining hits (ms): " + chainTime);
	out.close();
}
/*
 * Chains together database hits:
 * 
 * Input: List of hits as well as range [b, e] of hits which came from same 
 * sequence and relative strand, and whose difference in target vs. query 
 * positions are within EPS basepairs, and what the relative strand is
 * 
 * Output: Information about the maximal chain of alignments - that is, the
 * largest possible set of alignments such that when sorted by 
 * increasing query position, their target positions are also sorted 
 * (increasing if strand = 0 and decreasing if strand = 1)
 * 
 * Output is {queryStart, queryEnd, targetStart, targetEnd, number of hits}
 */
static int[] alignmentChain(LongList hits, int b, int e, boolean reverse)
{
	// Get hits in range [b, e] and sort by pos1 (query position)
	int n = e - b + 1;
	long[] hitList = new long[n];
	for(int i = 0; i<n; i++)
	{
		hitList[i] = hits.a[i+b] & ((1L << (2*posBits)) - 1);
	}
	Arrays.sort(hitList);
	
	// Split pos1 and pos2 into separate arrays for faster access
	int[] p1s = new int[n], p2s = new int[n];
	for(int i = 0; i<n; i++)
	{
		p1s[i] = getHitP1(hitList[i]);
		p2s[i] = getHitP2(hitList[i]);
	}
	
	// Solve Longest Increasing Subsequence with binary search
	// fastDp[i] will be the index x such that the LIS of all
	// elements up through p2s[x] (which must include p2s[x]) is
	// equal to i+1, with ties being broken based on prioritizing smaller
	// values of p2s[x].  This array will be updated as the algorithm
	// sweeps through the input from left to right.
	int[] fastDp = new int[n+1];
	Arrays.fill(fastDp, -1);
	
	// last[x] is index of the element coming before x in LIS ending at x
	int[] last = new int[n+1];
	Arrays.fill(last, -1);
	fastDp[0] = 0;
	for(int i = 1; i<n; i++)
	{
		// Binary search for how long the LIS ending at i can be
		int lo = -1, hi = i;
		while(lo < hi - 1)
		{
			int mid = (lo+hi)/2;
			boolean good = false;
			if(fastDp[mid] != -1)
			{
				good |= !reverse && (p2s[i] > p2s[fastDp[mid]]);
				good |= reverse && (p2s[i] < p2s[fastDp[mid]]);
			}
			if(good) lo = mid;
			else hi = mid;
		}
		
		// If it's smaller than the existing element at fastDp[lo+1], update it
		if(fastDp[lo+1] == -1 
			|| (!reverse && (p2s[i] <= p2s[fastDp[lo+1]])) 
			|| (reverse && (p2s[i] >= p2s[fastDp[lo+1]])))
		{
			fastDp[lo+1] = i;
			last[i] = (lo < 0) ? -1 : fastDp[lo];
		}
	}
	
	// Count how high the LIS can get before the -1's in fastDp are reached
	int lis = 1, maxi = fastDp[0];
	while(fastDp[lis] != -1)
	{
		lis++;
		maxi = fastDp[lis-1];
	}
	
	// Trace back to find the first element - not necessarily fastDp[0]
	int col = maxi;
	while(last[col] != -1) col = last[col];
	
	// Make sure there are at least c minimizer matches
	if(lis < c) return new int[] {-1};
	int endPos2 = p2s[maxi] + k - 1;
	int endPos = p1s[maxi] + k - 1;
	int startPos = p1s[col];
	int startPos2 = p2s[col];
	
	// Make sure startPos2 < endPos2 regardless of strand
	if(startPos2 > endPos2)
	{
		int tmp = startPos2;
		startPos2 = endPos2;
		endPos2 = tmp;
	}
	
	// Make sure the total length of the chain is at least minOverlap
	if(endPos2 - startPos2 + 1 < minOverlap) return new int[] {-1};
	
	return new int[] {startPos, endPos, startPos2, endPos2, lis};
}
/*
 * The following few methods rely on representing a "hit" as a 64-bit integer
 * A hit is when a minimizer in a query matches an entry in a minimizer database
 * The hit consists of the target sequence, the relative strand,
 * as well as the positions in the query (p1) and target (p2).
 * 
 * getHit takes in these 4 values and outputs a 64-bit integer encoding them
 */
static long getHit(int seq, int strand, int p1, int p2)
{
	return (((long)seq) << (2 * posBits+1)) 
		+ (((long)strand) << (2 * posBits)) 
		+ (((long)p1) << posBits) + p2;
}
/*
 * Extracts the strand from an encoded hit
 */
static int getHitStrand(long hit)
{
	return (int)((hit >> (2*posBits)) & 1);
}
/*
 * Extracts the target sequence id from an encoded hit
 */
static int getHitSeq(long hit)
{
	return (int)(hit >> (2*posBits+1));
}
/*
 * Extracts the query position from an encoded hit
 */
static int getHitP1(long hit)
{
	return (int)((hit&pos1Mask) >> posBits);
}
/*
 * Extracts the target position from an encoded hit
 */
static int getHitP2(long hit)
{
	return (int)(hit&pos2Mask);
}
/*
 * Extracts the difference in target vs. query positions from an encoded hit
 * This is the amount the query must be shifted so the shared kmers line up
 * 
 *  Note that if the sequences are on the same strand this is p1 - p2, 
 *  while if they are on different strands it is p1 + p2.
 */
static int getHitDiff(long hit)
{
	int st = getHitStrand(hit);
	if(st == 0) return getHitP1(hit) - getHitP2(hit);
	else return getHitP1(hit) + getHitP2(hit);
}
/*
 * Compares two encoded hits based on (in order)
 * 1. The target sequence id
 * 2. The relative strand
 * 3. The difference in query vs. target positions
 */
static int compareHits(long a, long b)
{
	int doublePosShift = 2*posBits;
	long shiftA = (a >> doublePosShift);
	long shiftB = (b >> doublePosShift);
	if(shiftA != shiftB) return Long.compare(shiftA, shiftB);
	return getHitDiff(a) - getHitDiff(b);
}
/*
 * Gets a list of (k, w)-minimizers from a string s
 */
static long[] getMinimizers(String s)
{
	// Constants
	long lmax = Long.MAX_VALUE; // "Infinity"
	long mask = (1L << (2*k)) - 1; // Bitmask of 2k 1's
	int hashShift = posBits + 1; // Bits the hash value is shifted in output
	
	// Queue for keeping current window - stored as a circular array of size w
	long[] qPos = new long[w]; // Encoding of position/hash/strand at a position
	long[] qHash = new long[w]; // The minimum hash value among the two strands
	int qIdx = 0; // The index where the next hash value will go in the queue
	
	// Keeping track of minimum value in window
	long minPos = 0; // The encoded position info of the minimum hash
	long minVal = lmax; // The value of the minimum hash in the window
	int mindex = -1; // The index in the queue of the minimum hash
	
	LongList res = new LongList();
	
	// Sliding window values of hashes of string and its reverse complement
	long[] hash = new long[2];
	
	int n = s.length();
	for(int i = 0; i<n; i++)
	{
		// Update both hashes
		hash[1] = (hash[1] >> 2) | ((3^map[s.charAt(i)]) << (2*k-2));
		hash[0] = ((hash[0] << 2) & mask) | map[s.charAt(i)];
		
		// Which strand has a smaller kmer at this position
		int z = hash[0] < hash[1] ? 0 : 1;
		
		// Hash after going through invertible hash function
		long xHash = hash(hash[z], mask);
		
		// Encoding of hash/pos/strand
		long pos = ((i-k+1)<<1) | z | (xHash << hashShift);
		
		if(i >= k-1)
		{
			// Start adding to queue once we have a full kmer
			qPos[qIdx] = pos;
			qHash[qIdx] = xHash;
		}
		if(i == w + k - 2)
		{
			// End of 1st window - add all values which are min value or better
			for(int j = qIdx + 1; j<w; j++)
				if(minVal == qHash[j] && minPos != qPos[j]) res.add(qPos[j]);
			for(int j = 0; j<qIdx; j++)
				if(minVal == qHash[j] && minPos != qPos[j]) res.add(qPos[j]);
		}
		if(xHash <= minVal)
		{
			// New minimum value, so add old minimum to answer list
			if(i >= w + k) res.add(minPos);
			
			// Update minimum to point to this value
			minPos = pos;
			minVal = xHash;
			mindex = qIdx;
		}
		else if(mindex == qIdx)
		{
			// We just overwrote the minimum - add it and figure out new minimum
			if(i >= w + k - 2)
			{
				res.add(minPos);
			}
			// Find new min - prioritize more recently added in case of tie
			minVal = lmax;
			for(int j = qIdx + 1; j<w; j++)
				if(minVal >= qHash[j])
				{
					mindex = j;
					minPos = qPos[j];
					minVal = qHash[j];
				}
			for(int j = 0; j<=qIdx; j++)
				if(minVal >= qHash[j])
				{
					mindex = j;
					minPos = qPos[j];
					minVal = qHash[j];
				}
			// Add any hashes which are equal to the minimum
			if(i >= w + k - 2)
			{
				for(int j = qIdx + 1; j<w; j++)
					if(minVal == qHash[j] && minPos != qPos[j])
					{
						res.add(qPos[j]);
					}
				for(int j = 0; j<=qIdx; j++)
					if(minVal == qHash[j] && minPos != qPos[j])
					{
						res.add(qPos[j]);
					}
			}
		}
		
		// Move queue index
		qIdx++;
		if(qIdx == w) qIdx = 0;
		
	}
	if(minVal != lmax) res.add(minPos);
	long[] ans = new long[res.size];
	for(int i = 0; i<ans.length; i++) ans[i] = res.a[i];
	return ans;
}
/*
 * Invertible hash function used in Minimap to avoid prioritizing poly-A tails
 */
static long hash(long val, long m)
{
	long x = (~val + (val << 21)) & m;
	x = x ^ (x >> 24);
	x = (x + (x<<3) + (x<<8)) & m;
	x = x ^ (x >> 14);
	x = (x + (x<<2) + (x<<4)) & m;
	x = x ^ (x >> 28);
	x = (x + (x << 31)) & m;
	return x;
}
/*
 * Encodes a hash value, position, and strand as a single 64-bit integer
 */
static long getMinimizer(long hash, int pos, int strand)
{
	return strand + (pos << 1) + (hash << (posBits + 1));
}
/*
 * A database of Minimizers represented as integer encoding seq/pos/strand
 */
static class MinimizerDatabase
{
	HashMap<Long, Integer> startRange;
	long[] hashes;
	long[] poss;
	int size;
	MinimizerDatabase()
	{
		size = 0;
		hashes = new long[16];
		poss = new long[16];
	}
	void put(long hash, long pos)
	{
		if(size == hashes.length)
		{
			long[] nh = new long[size<<1];
			long[] np = new long[size<<1];
			for(int i = 0; i<size; i++)
			{
				nh[i] = hashes[i];
				np[i] = poss[i];
			}
			hashes = nh;
			poss = np;
		}
		hashes[size] = hash;
		poss[size] = pos;
		size++;
	}
	/*
	 * Sort database and map minimizers to the start of range where they occur
	 */
	void finish()
	{
		long time = System.currentTimeMillis();
		long[][] a = new long[][] {hashes, poss};
		long[][] sorted = sort(a);
		hashes = sorted[0];
		poss = sorted[1];
		long ntime = System.currentTimeMillis();
		System.out.println("Time sorting database (ms): " + (ntime - time));
		startRange = new HashMap<Long, Integer>();
		for(int i = 0; i<size; i++)
		{
			if(i == 0 || hashes[i] != hashes[i-1]) startRange.put(hashes[i], i);
		}
		long finalTime = System.currentTimeMillis();
		System.out.println("Time hashing database (ms): " + (finalTime - ntime));
	}
	/*
	 * Merge sort which sorts two lists according to the order of the first list
	 */
	long[][] sort(long[][] a)
	{
		int n = a[0].length;
		if(n == 1) return a;
		int n1 = n/2, n2 = n - n1;
		long[][] a1 = new long[2][n1];
		long[][] a2 = new long[2][n2];
		for(int i = 0; i<n1; i++)
		{
			a1[0][i] = a[0][i];
			a1[1][i] = a[1][i];
		}
		for(int i = 0; i<n2; i++)
		{
			a2[0][i] = a[0][i+n1];
			a2[1][i] = a[1][i+n1];
		}
		a1 = sort(a1);
		a2 = sort(a2);
		return merge(a1, a2);
	}
	/*
	 * Merge function used in merge sort
	 */
	long[][] merge(long[][] a1, long[][] a2)
	{
		int n1 = a1[0].length, n2 = a2[0].length;
		int n = n1 + n2;
		long[][] res = new long[2][n];
		int i = 0, j = 0, k = 0;
		while(i < n1 || j < n2)
		{
			if(i == n1)
			{
				res[0][k] = a2[0][j];
				res[1][k++] = a2[1][j++];
			}
			else if(j == n2)
			{
				res[0][k] = a1[0][i];
				res[1][k++] = a1[1][i++];
			}
			else if(a1[0][i] <= a2[0][j])
			{
				res[0][k] = a1[0][i];
				res[1][k++] = a1[1][i++];
			}
			else
			{
				res[0][k] = a2[0][j];
				res[1][k++] = a2[1][j++];
			}
		}
		return res;
	}
}
/*
 * A faster ArrayList<Long> with limited functionality
 */
static class LongList
{
	int size;
	long[] a;
	LongList(int n)
	{
		a = new long[n];
		size = 0;
	}
	LongList()
	{
		a = new long[16];
		size = 0;
	}
	void add(long x)
	{
		if(size == a.length)
		{
			long[] b = new long[size<<1];
			for(int i = 0; i<size; i++) b[i] = a[i];
			a = b;
		}
		a[size++] = x;
	}
	void sort(Comparator<Long> cl)
	{
		Long[] res = new Long[size];
		for(int i = 0; i<size; i++) res[i] = a[i];
		Arrays.sort(res, cl);
		for(int i = 0; i<size; i++) a[i] = res[i];
	}
}
}
