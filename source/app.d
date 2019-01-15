import std.stdio;
import std.algorithm.searching:until;
import std.path:buildPath;
import std.conv:to;

import bio.std.hts.bam.reader:BamReader;
import bio.std.hts.bam.pileup;

import csv;
void main(string[] args)
{
	if(args.length!=5){
		writeln("usage: ./depthgauge [input tsv] [number of first sample column] [bam folder] [output tsv]");
	}else{
		auto t =Table(args[1],args[2].to!int-1);
		getDepths(t,args[3]);
		File f = File(args[4],"w");
		t.write(f);
		f.close;
	}

}

auto depth_at_pos(ref BamReader bam,string chr,uint pos){
	return bam[chr][pos..pos+1].makePileup(true,pos,pos,false).front.coverage;
}

void getDepths(ref Table t,string prefix){
	foreach(j,sample;t.samples){
		auto bam = new BamReader(buildPath(prefix,sample~".bam"));
		foreach(i,rec;t.records){
			t.matrix[i][j]=depth_at_pos(bam,rec.chr.idup,rec.pos);
		}
	}
}




