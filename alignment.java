import java.util.*;
import java.io.*;
public class alignment {
	static int k = 15, w = 5;
	static int[] map;
	static int EPS = 500;
	static double f = 0.001;
	static int batchSize = 6000;
	static int minOverlap = 100;
	static int c = 4;
	static int posBits = 20;
	static boolean writeFile = true;
	static boolean overallFilter = true;
public static void main(String[] args) throws IOException
{
	
	long minimizerTime = 0;
	long startTime = System.currentTimeMillis();
	PrintWriter out = new PrintWriter(System.out);
	if(writeFile) out = new PrintWriter(new File("/home/mkirsche/ecoli/oxford.paf"));
	map = new int[256];
	map['A'] = 0; map['C'] = 1; map['G'] = 2; map['T'] = 3;
	//String ss = "AGTGCGTGAGCGTGCGAACGT";
	//String tt = "ACGTTCGCACGCTCACGCACT";
	//HashSet<Long> set1 = getMinimizers(ss);
	//HashSet<Long> set2 = getMinimizers(tt);
	//for(Long a : set1) System.out.println(Long.toBinaryString(a));
	//System.out.println();
	//for(Long b : set2) System.out.println(Long.toBinaryString(b));
	String fn = "/home/mkirsche/ecoli/oxford.fasta";
	Scanner input = new Scanner(new FileInputStream(new File(fn)));
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
	HashSet<Integer> toRemove = new HashSet<Integer>();
	if(overallFilter)
	{
		System.out.println("Counting minimizer frequencies");
		HashMap<Integer, Integer> overallFreq = new HashMap<Integer, Integer>();
		for(int i = 0; i<n; i++)
		{
			String read = reads.get(i);
			long curTime = System.currentTimeMillis();
			HashSet<Long> cur = getMinimizers(read);
			minimizerTime += System.currentTimeMillis() - curTime;
			for(long x : cur)
			{
				int keepBits = posBits + 1;
				int minimizerValue = (int)(x >> keepBits);
				overallFreq.put(minimizerValue, 
						overallFreq.containsKey(minimizerValue) ? (1 + overallFreq.get(minimizerValue)) : 1);
			}
		}
		System.out.println("Number of distinct minimizers: " + overallFreq.size());
		ArrayList<Integer> freqMinimizers = new ArrayList<Integer>();
		for(int x : overallFreq.keySet())
		{
			freqMinimizers.add(overallFreq.get(x));
		}
		Collections.sort(freqMinimizers);
		int repeatedMinimizers = (int)(f * freqMinimizers.size());
		int repeatThreshold = freqMinimizers.get(freqMinimizers.size() - 1 - repeatedMinimizers);
		for(int x : overallFreq.keySet())
		{
			if(overallFreq.get(x) > repeatThreshold)
				toRemove.add(x);
		}
		overallFreq.clear();
		System.out.println("Only keeping minimizers occurring at most " + repeatThreshold + " times");
	}
	
	for(int iter = 0; iter<n; iter += batchSize)
	{
		System.out.println("Processing batch " + (iter/batchSize + 1) + " of " + (n + batchSize - 1) / batchSize);
		HashMap<Long, ArrayList<Long>> map = new HashMap<>();
		System.out.println("Computing minimizer database");
		for(int j = iter; j < n && j < iter + batchSize; j++)
		{
			String read = reads.get(j);
			long curTime = System.currentTimeMillis();
			HashSet<Long> cur = getMinimizers(read);
			minimizerTime += System.currentTimeMillis() - curTime;
			for(long x : cur)
			{
				int keepBits = posBits + 1;
				long key = (x & ((1L << keepBits) - 1)) + (((long)readsProcessed) << keepBits);
				long minimizerValue = x >> keepBits;
				if(overallFilter && toRemove.contains((int)minimizerValue)) continue;
				
				if(!map.containsKey(minimizerValue)) map.put(minimizerValue, new ArrayList<Long>());
				map.get(minimizerValue).add(key);
			}
			readsProcessed++;
		}
		int repeatThreshold = -1;
		if(!overallFilter)
		{
			System.out.println("Filtering repeat minimizers");
			ArrayList<Integer> freqMinimizers = new ArrayList<Integer>();
			for(long x : map.keySet()) freqMinimizers.add(map.get(x).size());
			Collections.sort(freqMinimizers);
			int repeatedMinimizers = (int)(f * freqMinimizers.size());
			repeatThreshold = freqMinimizers.get(freqMinimizers.size() - 1 - repeatedMinimizers);
			System.out.println("Only keeping minimizers occurring at most " + repeatThreshold + " times");
		}
		System.out.println("Querying reads against database");
		for(int i = 0; i<n; i++)
		{
			if(i > 0 && i%5000 == 0) System.out.println("Processed " + i + " queries");
			String read = reads.get(i);
			long curTime = System.currentTimeMillis();
			HashSet<Long> cur = getMinimizers(read);
			minimizerTime += System.currentTimeMillis() - curTime;
			ArrayList<Hit> hs = new ArrayList<Hit>();
			for(long x : cur)
			{
				int keepBits = posBits + 1;
				int xstrand = (int)(x&1);
				int xpos = (int)((x & ((1L << keepBits) - 1)) >> 1);
				long minimizerValue = x >> keepBits;
				ArrayList<Long> matches = map.get(minimizerValue);
				if(matches == null) continue;
				if(!overallFilter && matches.size() > repeatThreshold) continue;
				for(long y : matches)
				{
					int ystrand = (int)(y&1);
					int ypos = (int)((y & ((1L << keepBits) - 1)) >> 1);
					int seqId = (int)(y >> keepBits);
					if(seqId == i) continue;
					int nStrand = xstrand ^ ystrand;
					Hit h = new Hit(seqId, nStrand, xpos, ypos);
					hs.add(h);
				}
			}
			Collections.sort(hs);
			int b = 0;
			for(int e = 0; e<hs.size(); e++)
			{
				if(e == hs.size() - 1 || hs.get(e).seq != hs.get(e+1).seq || hs.get(e).strand != hs.get(e+1).strand || hs.get(e+1).diff - hs.get(e).diff > EPS)
				{
					if(e >= b + c - 1)
					{
						int[] lis = alignmentChain(hs, b, e, hs.get(e).strand > 0);
						if(lis[0] != -1)
						{
							String queryName = names.get(i);
							int queryLength = reads.get(i).length();
							int queryStart = lis[0];
							int queryEnd = lis[1];
							char relativeStrand = hs.get(e).strand == 1 ? '+' : '-';
							int j = hs.get(e).seq;
							String targetName = names.get(j);
							int targetLength = reads.get(j).length();
							int targetStart = lis[2];
							int targetEnd = lis[3];
							int residueMatches = lis[4];
							int alignmentBlockLength = k;
							int mappingQuality = 255;
							
							out.printf("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n", 
									queryName, queryLength, queryStart, queryEnd,
									relativeStrand, targetName, targetLength, targetStart, targetEnd,
									residueMatches, alignmentBlockLength, mappingQuality);
							
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
	out.close();
}
static int[] alignmentChain(ArrayList<Hit> hits, int b, int e, boolean reverse)
{
	int n = e - b + 1;
	Hit[] hitList = new Hit[n];
	for(int i = 0; i<n; i++) hitList[i] = hits.get(i+b);
	Arrays.sort(hitList, new HitComparator());
	int[][] dp = new int[n][n];
	dp[0][0] = 1;
	for(int i = 1; i<n; i++)
	{
		for(int j = 0; j<i; j++) dp[i][j] = dp[i-1][j];
		dp[i][i] = 1;
		for(int j = 0; j<i; j++)
		{
			boolean good = !reverse && hitList[i].pos2 > hitList[j].pos2;
			good |= reverse && hitList[i].pos2 < hitList[j].pos2;
			if(good)
				dp[i][i] = Math.max(dp[i][i], 1 + dp[i-1][j]);
		}
	}
	int max = 0, maxi = -1;
	for(int i = 0; i<n; i++)
	{
		if(dp[n-1][i] > max)
		{
			max = dp[n-1][i];
			maxi = i;
		}
	}
	int lis = max;
	if(lis < c) return new int[] {-1};
	int endPos2 = hitList[maxi].pos2 + k - 1;
	int endPos = hitList[maxi].pos + k - 1;
	int startPos = hitList[maxi].pos;
	int startPos2 = hitList[maxi].pos2;
	int col = maxi;
	int row = n-1;
	rowloop: while(row > 0)
	{
		if(dp[row][col] == dp[row-1][col])
		{
			row--;
			continue;
		}
		for(int j = col-1; j >= 0; j--)
		{
			boolean good = !reverse && hitList[col].pos2 > hitList[j].pos2 && dp[row][col] == 1 + dp[row-1][j];
			good |= reverse && hitList[col].pos2 < hitList[j].pos2 && dp[row][col] == 1 + dp[row-1][j];
			if(good)
			{
				row--;
				col = j;
				endPos2 = Math.max(endPos2, hitList[j].pos2 + k - 1);
				endPos = Math.max(endPos, hitList[j].pos + k - 1);
				startPos = Math.min(startPos, hitList[j].pos);
				startPos2 = Math.min(startPos2, hitList[j].pos2);
				continue rowloop;
			}
		}
		break;
	}
	if(endPos2 - startPos2 + 1 < minOverlap) return new int[] {-1};
	return new int[] {startPos, endPos, startPos2, endPos2, lis};
}
static class HitComparator implements Comparator<Hit>
{
    public int compare(Hit a, Hit b)
    {
        return a.pos - b.pos;
    }
}
static class Hit implements Comparable<Hit>
{
	int seq;
	int strand;
	int pos;
	int pos2;
	int diff;
	Hit(int se, int ss, int pp, int pp2)
	{
		seq = se;
		strand = ss;
		pos = pp;
		pos2 = pp2;
		diff = (strand == 0) ? (pos - pos2) : (pos + pos2);
	}
	public int compareTo(Hit o) {
		if(seq != o.seq) return seq - o.seq;
		if(strand != o.strand) return strand - o.strand;
		if(diff != o.diff) return diff - o.diff;
		if(pos2 != o.pos2) return pos2 - o.pos2;
		return 0;
	}
	
}
static HashSet<Long> getMinimizers(String s)
{
	long mask = (1L << (2*k)) - 1;
	HashSet<Long> res = new HashSet<Long>();
	long hash = 0, revHash = 0;
	int n = s.length();
	long[] hashes = new long[n-k+1];
	long[] revHashes = new long[n-k+1];
	for(int i = 0; i<k; i++)
	{
		hash = (hash << 2) | map[s.charAt(i)];
		revHash = (revHash << 2) | (3^map[s.charAt(k-1-i)]);
	}
	hashes[0] = hash(hash, 2*k);
	revHashes[0] = hash(revHash, 2*k);
	for(int i = k; i<n; i++)
	{
		revHash = (revHash >> 2) | ((3^map[s.charAt(i)]) << (2*k-2));
		hash = ((hash << 2) & mask) | map[s.charAt(i)];
		hashes[i-k+1] = hash(hash, 2*k);
		revHashes[i-k+1] = hash(revHash, 2*k);
	}
	for(int i = 0; i<=n-k-w; i++)
	{
		long min = Long.MAX_VALUE;
		for(int j = 0; j<w; j++)
		{
			if(hashes[i+j] != revHashes[i+j])
			{
				min = Math.min(min, Math.min(hashes[i+j], revHashes[i+j]));
			}
		}
		for(int j = 0; j<w; j++)
		{
			if(hashes[i+j] != revHashes[i+j])
			{
				if(hashes[i+j] == min) res.add(getMinimizer(min, i+j, 0));
				if(revHashes[i+j] == min) res.add(getMinimizer(min, i+j, 1));
			}
		}
	}
	return res;
}
static long hash(long val, int p)
{
	long m = (1L<<p) - 1;
	long x = (~val + (val << 21)) & m;
	x = x ^ (x >> 24);
	x = (x + (x<<3) + (x<<8)) & m;
	x = x ^ (x >> 14);
	x = (x + (x<<2) + (x<<4)) & m;
	x = x ^ (x >> 28);
	x = (x + (x << 31)) & m;
	return x;
}
static long getMinimizer(long hash, int pos, int strand)
{
	return strand + (pos << 1) + (hash << (posBits + 1));
}
}
