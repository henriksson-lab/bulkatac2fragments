package bulkatac2fragments;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import java.io.PrintStream;
import java.io.PrintWriter;
import java.util.List;
import java.util.PriorityQueue;
import java.util.StringTokenizer;

public class NewMergeBed {
	
	
	public static class OneBed implements Comparable<OneBed> {
		
		BufferedReader br;
		String line;
		
		String chr;
		int from;
		int to;
		
		
		String name;
		
		
		public OneBed(File f) throws FileNotFoundException {
			br = new BufferedReader(new FileReader(f));			
			name=f.getName();
		}
		
		public boolean nextLine() throws IOException {
			
			line=br.readLine();
			if(line==null) {
				return false;
			} else {
				
				StringTokenizer stok=new StringTokenizer(line,"\t");
				
				chr=stok.nextToken();
				from=Integer.parseInt(stok.nextToken());
				to=Integer.parseInt(stok.nextToken());
				
				return true;
			}
		}
		
		
		
		@Override
		public int compareTo(OneBed o) {
			int c=chr.compareTo(o.chr);
			if(c<0) {
				return -1;
			} else if(c>0) {
				return 1;
			} else {
				if(from<o.from) {
					return -1;
				} else if(from>o.from) {
					return 1;
				} else {
					
					if(to<o.to) {
						return -1;
					} else if(to>o.to) {
						return 1;
					} else {
						return 0;
					}
					
				}
			}
		}
		
		
		public StringBuffer getLine() {
			StringBuffer sb=new StringBuffer();
			sb.append(chr);
			sb.append("\t");
			sb.append(from);
			sb.append("\t");
			sb.append(to);
			sb.append("\t");
			sb.append(name);
			sb.append("\t");
			sb.append(1);  //number of fragments
			return sb;
		}
		
	}
	
	
	
	
	public static void merge(List<String> args, PrintWriter ps) throws IOException {
		

		int numline=0;
		
		//Place all files in queue
		PriorityQueue<OneBed> q = new PriorityQueue<NewMergeBed.OneBed>();
		for(String s:args) {
			OneBed b=new OneBed(new File(s));
			if(b.nextLine()) {
				q.add(b);
			}
		}
		
		//Get all lines in order
		while(!q.isEmpty()) {
			OneBed b=q.poll();
			ps.println(b.getLine());
			if(b.nextLine()) {
				q.add(b);
			}
			numline++;
			
			if(numline%1000000==0) {
				System.out.println("numline "+numline+ "  current files "+q.size());
			}
		}
		
	}
	
	
	
	/*
	public static void main(String[] args) throws IOException {

		
		
		//Place all files in queue
		PriorityQueue<OneBed> q = new PriorityQueue<NewMergeBed.OneBed>();
		for(String s:args) {
			OneBed b=new OneBed(new File(s));
			if(b.nextLine()) {
				q.add(b);
			}
		}
		
		//Get all lines in order
		while(!q.isEmpty()) {
			OneBed b=q.poll();
			System.out.println(b.getLine());
			if(b.nextLine()) {
				q.add(b);
			}
		}
		
		//System.out.println("done");
		
		
		
	}*/

}
