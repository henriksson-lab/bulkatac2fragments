package bulkatac2fragments;


import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.io.InputStreamReader;
import java.io.PrintWriter;
import java.util.TreeMap;


public class BamToFragments {
	
	
	
	public static void read(File fBAM, String bc, PrintWriter pw) throws IOException {

		TreeMap<String, String[]> previousRead=new TreeMap<String, String[]>();

		Process p = Runtime.getRuntime().exec("samtools view -bf 0x2 "+fBAM+" | bedtools bamtobed -i stdin");
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
					if(!chr1.equals(lastChr) && lastStart!=start && lastEnd!=end) {
						pw.println(chr1+"\t"+start+"\t"+end+"\t"+bc+"\t"+1);
					}					
				}
				
			} else {
				previousRead.put(readname, parts);
			}

		}
		inp.close();
		
		System.out.println("total records read: "+readRecords+"    num reads to be matched up: "+previousRead.size());
	}
	
	
	
	
	
	
	public static void main(String[] args) {

		if(args.length==0) {
			System.out.println("This software takes paired-end ATAC-seq files and generes 10x-like fragment files.");
			System.out.println("Thus they can be analyzed using single-cell workflows (such as signac or archr)");
			System.out.println();
			System.out.println("Arguments: input.bam [output.tsv [barcode_to_write]]");
			System.out.println();
			System.out.println("Later merge all files using: bedtools merge -i *.atac_fragments.tsv > atac_fragments.tsv");
			System.out.println("Compress it: bgzip -@ 8 -i atac_fragments.tsv");
			//System.out.println("Index it: tabix -p vcf atac_fragments.tsv.gz");
			System.exit(0);
		}
		
		
		File fin=new File(args[0]);
		
		File fout=new File(args[0]+".atac_fragments.tsv");
		if(args.length>1)
			fout=new File(args[1]);
		
		String bc=fin.getName();
		if(args.length>2)
			bc=args[2];		
		
		try {
			PrintWriter pw=new PrintWriter(new BufferedWriter(new FileWriter(fout)));
			
			read(fin, bc, pw);
			
			pw.close();
		} catch (IOException e) {
			e.printStackTrace();
		}
		
		
	}

}
