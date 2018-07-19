import java.util.*;
import java.io.*;
public class containment {
	static boolean local = false;
public static void main(String[] args) throws IOException
{
	if(!local && args.length != 3)
	{
		System.out.println("Usage: java FileIntersect <pafFn> <readsFn> <outFn>");
		return;
	}
	String pafFn = "/home/mkirsche/ecoli/oxford.paf";
	String readFn = "/home/mkirsche/ecoli/oxford.fasta";
	String outFn = "/home/mkirsche/ecoli/oxford_filtered_miniasmclone.fasta";

	if(args.length > 0)
	{
		pafFn = args[0];
		readFn = args[1];
		outFn = args[2];
	}
	
	Scanner pafInput = new Scanner(new FileInputStream(new File(pafFn)));
	Scanner readInput = new Scanner(new FileInputStream(new File(readFn)));
	HashSet<String> contained = new HashSet<String>();
	while(pafInput.hasNext())
	{
		PafEntry cur = new PafEntry(pafInput.nextLine());
		// Check if contained
		int[] ls = new int[] {cur.queryLength, cur.targetLength};
		int[] bs = new int[] {cur.queryStart, cur.targetStart};
		int[] es = new int[] {cur.queryEnd, cur.targetEnd};
		int whichContained = contained(ls, bs, es);
		if(whichContained == 1) contained.add(cur.queryName);
		else if(whichContained == 2) contained.add(cur.targetName);
	}
	System.out.println(contained.size());
	PrintWriter out = new PrintWriter(new File(outFn));
	int tot = 0;
	while(readInput.hasNext())
	{
		tot++;
		String name = readInput.next().substring(1);
		readInput.nextLine();
		String read = readInput.nextLine();
		if(contained.contains(name)) continue;
		out.println(name);
	}
	out.close();
	System.out.println("Filtered " + contained.size() + " of " + tot + " reads");
}
static int o = 1000; // Max overhang length
static double r = 0.8; // Max overhang to mapping length ratio
/*
 * Returns 0 if neither is contained, 1 if first read is contained, and 2 if second read is contained
 */
static int contained(int[] ls, int[] bs, int[] es)
{
	int overhang = Math.min(bs[0], bs[1]) + Math.min(ls[0] - es[0], ls[1] - es[1]); // 
	int mapLen = Math.max(es[0] - bs[0], es[1] - bs[1]);
	if(overhang > Math.min(o, mapLen * r)) return 0; // Internal match but not contained
	if(bs[0] <= bs[1] && ls[0] - es[0] <= ls[1] - es[1]) return 1;
	if(bs[1] <= bs[0] && ls[1] - es[1] <= ls[0] - es[0]) return 2;
	return 0;
}
static class PafEntry
{
	String queryName, targetName;
	int queryStart, queryEnd, targetStart, targetEnd;
	char relativeStrand;
	int targetLength, queryLength;
	int residueMatches, alignmentBlockLength, mappingQuality;
	PafEntry(String s)
	{
		StringTokenizer st = new StringTokenizer(s, "\t");
		queryName = st.nextToken();
		queryLength = Integer.parseInt(st.nextToken());
		queryStart = Integer.parseInt(st.nextToken());
		queryEnd = Integer.parseInt(st.nextToken());
		relativeStrand = st.nextToken().charAt(0);
		targetName = st.nextToken();
		targetLength =  Integer.parseInt(st.nextToken());
		targetStart = Integer.parseInt(st.nextToken());
		targetEnd = Integer.parseInt(st.nextToken());
		residueMatches = Integer.parseInt(st.nextToken());
		alignmentBlockLength = Integer.parseInt(st.nextToken());
		mappingQuality = Integer.parseInt(st.nextToken());
		if(queryStart > queryEnd)
		{
			int tmp = queryStart;
			queryStart = queryEnd;
			queryEnd = tmp;
		}
		if(targetStart > targetEnd)
		{
			int tmp = targetStart;
			targetStart = targetEnd;
			targetEnd = tmp;
		}
	}
}
}
