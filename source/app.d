import std.stdio;
import std.algorithm.searching:until;
import std.algorithm.iteration:filter;
import std.algorithm:map;
import std.range:array;
import std.path:buildPath;
import std.conv:to;

import bio.std.hts.bam.reader:BamReader;
import bio.std.hts.bam.multireader:MultiBamReader;
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

auto depth_at_pos(ref MultiBamReader bam,string sample,string chr,uint pos){
	return bam[chr][pos..pos+1].filter!(x=>x["RG"]==sample).makePileup(true,pos,pos,false).front.coverage;
}

void getDepths(ref Table t,string prefix){
	auto bam = new MultiBamReader(t.samples.map!(x=>buildPath(prefix,x~".bam")).array);
	foreach(j,sample;t.samples){
		foreach(i,rec;t.records){
			t.matrix[i][j]=depth_at_pos(bam,sample,rec.chr.idup,rec.pos);
		}
	}
}




