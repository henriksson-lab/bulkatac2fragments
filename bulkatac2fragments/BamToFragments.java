package bulkatac2fragments;


import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileWriter;
import java.io.IOException;
import java.io.InputStreamReader;
import java.io.PrintWriter;
import java.util.LinkedList;
import java.util.PriorityQueue;
import java.util.TreeMap;


public class BamToFragments {
	
	//Note. this is broken. need to postpone writing until pair found
	
	//One priority queue? to keep R1. https://docs.oracle.com/javase/8/docs/api/java/util/PriorityQueue.html
	//Constant time to check if paired read in there
	//Take first item to check safe position
	//If first item too far away then drop it. add a warning read to set. drop R2 when it appears
	
	//Put outgoing in a queue. output while still safe to do so
	
	
	public static void read(File fBAM, String bc, PrintWriter pw) throws IOException {

		TreeMap<String, String[]> previousRead=new TreeMap<String, String[]>();

		String cmd[]= {
				"sh", "-c",
				"samtools view -bf 0x2 "+fBAM+" | bedtools bamtobed -i stdin"
		};
		
		Process p = Runtime.getRuntime().exec(cmd);
		BufferedReader inp = new BufferedReader( new InputStreamReader(p.getInputStream()) );

		
		String lastChr="";
		int lastStart=0, lastEnd=0;
		
		//Example line:
		//1       11840   11991   SRR10984490.12149842/2  1       +
		
		String line;
		int readRecords=0;
		while((line=inp.readLine())!=null) {
			//Update user about progress
			readRecords++;
			if(readRecords%1000000 == 0){
				//Calculate progress
				System.out.println("records so far: "+readRecords+"    num reads to be matched up: "+previousRead.size());
			}
			
			String[] parts=line.split("\t", 0);
			
			if(parts.length<2) {
				System.out.println(line);
				throw new IOException("Unexpected input. Example line above");
			}
			
			String readname=parts[3];
			
			//Remove trailing /1 or /2 for PE
			int thes=readname.indexOf('/');
			if(thes==-1)
				throw new RuntimeException("Missing / in read name. Is it really paired end?");
			readname=readname.substring(0,thes);
			
			
			//Try to match up with the other read
			if(previousRead.containsKey(readname)) {
				
				String[] otherRead=previousRead.remove(readname);

				//Only keep reads that are on the same chromosome
				String chr1=parts[0];
				String chr2=otherRead[0];
				if(chr1.equals(chr2)) {
					int start=Math.min(Integer.parseInt(parts[1]), Integer.parseInt(otherRead[1]));
					int end=Math.max(Integer.parseInt(parts[2]), Integer.parseInt(otherRead[2]));
					
					//Deduplicate reads. Here ignoring the count, could be improved in the future
					if(!chr1.equals(lastChr) || lastStart!=start || lastEnd!=end) {
						pw.println(chr1+"\t"+start+"\t"+end+"\t"+bc+"\t"+1);
						
						lastChr=chr1;
						lastStart=start;
						lastEnd=end;
					}					
				}
				
			} else {
				previousRead.put(readname, parts);
			}

		}
		inp.close();
		
		System.out.println("total records read: "+readRecords+"    num reads to be matched up: "+previousRead.size());
	}
	
	
	
	public static void printHelp() {
		System.out.println("This software takes paired-end ATAC-seq files and generes 10x-like fragment files.");
		System.out.println("Thus they can be analyzed using single-cell workflows (such as signac or archr)");
		System.out.println();
		System.out.println("java -jar bam2fragments.jar bam2bed input.bam [-o output.tsv] [-b barcode_to_write]]");
		System.out.println("java -jar bam2fragments.jar mergebed atac_fragments.tsv *.tsv");
		System.out.println();
		//System.out.println("Later merge all files using: bedtools merge -i *.atac_fragments.tsv > atac_fragments.tsv");
		System.out.println("Compress it: bgzip -@ 8 -i atac_fragments.tsv");
		//System.out.println("Index it: tabix -p vcf atac_fragments.tsv.gz");
		System.exit(0);
	}
	
	
	public static void main(String[] args) {

		if(args.length==0) {
			printHelp();
		}
		
		if(args[0].equals("bam2bed") && args.length>1) {
			
			File fin=new File(args[1]);
			File fout=new File(fin.getParentFile(), fin.getName()+".atac_fragments.tsv");
			String bc=fin.getName();
			
			for(int i=2;i<args.length;i++) {
				if(args[i].equals("-o")) {
					fout=new File(args[i+1]);
					i++;
				} else if(args[i].equals("-b")) {
					bc=args[i+1];		
					i++;
				} else {
					System.out.println("Parse error on "+args[i]);
					System.exit(0);
				}
			}
			
			try {
				PrintWriter pw=new PrintWriter(new BufferedWriter(new FileWriter(fout)));
				
				read(fin, bc, pw);
				
				pw.close();
			} catch (IOException e) {
				e.printStackTrace();
			}
			
		} else if(args[0].equals("mergebed") && args.length>1) {
			
			File fout=new File(args[1]);
			
			LinkedList<File> finList=new LinkedList<File>();
			
			for(int i=2;i<args.length;i++) {
				finList.add(new File(args[i]));
			}

			
			try {
				PrintWriter pw=new PrintWriter(new BufferedWriter(new FileWriter(fout)));
				
				merge(pw, finList);
				
				pw.close();
			} catch (IOException e) {
				e.printStackTrace();
			}
			
		} else if(args[0].equals("newmergebed") && args.length>1) {
			
			File fout=new File(args[1]);
			
			LinkedList<String> finList=new LinkedList<String>();
			
			for(int i=2;i<args.length;i++) {
				finList.add(args[i]);
			}

			
			try {
				PrintWriter pw=new PrintWriter(new BufferedWriter(new FileWriter(fout)));
				
				NewMergeBed.merge(finList, new PrintWriter(new BufferedWriter(pw)));
				
				pw.close();
			} catch (IOException e) {
				e.printStackTrace();
			}
		} else {
			System.out.println("Parse error on "+args[0]);
			System.exit(0);
		}
		
		
		
	}


	private static void merge(PrintWriter pw, LinkedList<File> finList) throws IOException {
		
		PriorityQueue<PrioritizedReader> q=new PriorityQueue<PrioritizedReader>();
		
		for(File f:finList) {
			PrioritizedReader r=new PrioritizedReader(f);
			if(r.readLine()) {
				q.add(r);
				System.out.println("Including "+f);
			}
		}
		
		System.out.println("Merging all");
		int readRecords=0;
		while(!q.isEmpty()) {
			PrioritizedReader r=q.poll();
			pw.println(r.line);
			if(r.readLine())
				q.add(r);
			readRecords++;
			if(readRecords%1000000 == 0){
				//Calculate progress
				System.out.println("records so far: "+readRecords+"    num files still going: "+q.size());
			}

		}
		
		System.out.println("done");

		
	}

}
