import std.stdio;
import std.algorithm.searching:until;
import std.algorithm.sorting:sort;
import std.algorithm.iteration:uniq,filter;
import std.algorithm:map,count;
import std.range:array;
import std.path:buildPath;
import std.conv:to;
import std.parallelism:parallel;
import dhtslib.sam:SAMFile;

import csv;
void main(string[] args)
{
	if(args.length!=5){
		writeln("usage: ./depthgauge [input tsv] [number of first sample column] [bam folder] [output tsv]");
	}else{
		auto t =Table(args[1],args[3],args[2].to!int-1);
		getDepths(t,args[3]);
		File f = File(args[4],"w");
		t.write(f);
		f.close;
	}

}

auto depth_at_pos(ref SAMFile bam,int chr,uint pos){
	//return bam[chr][pos..pos+1].makePileup(true,pos,pos,false).front.coverage;
	return bam.query(chr,pos,pos+1).count;
}

void getDepths(ref Table t,string prefix){
	foreach(j,sample;parallel(t.samples)){
		auto bam = SAMFile(buildPath(prefix,sample~".bam"),0);
		foreach(i,rec;t.records.sort.array){
			t.matrix[i][j]=depth_at_pos(bam,rec.chr,rec.pos);
		}
		writeln(sample);
	}
}




