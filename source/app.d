import std.stdio;
import std.algorithm.searching:until;
import std.algorithm.sorting:sort;
import std.algorithm:map;
import std.algorithm.iteration:uniq;
import std.path:buildPath;
import std.conv:to;
import std.range;
import std.parallelism:parallel,defaultPoolThreads;

import bio.std.hts.bam.reader:BamReader;
import bio.std.hts.bam.pileup;
import bio.std.hts.bam.region;

import csv;
void main(string[] args)
{
    //defaultPoolThreads(10);
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

auto depth_at_pos(ref BamReader bam,string chr,uint pos){
	return bam[chr][pos..pos+1].makePileup(true,pos,pos,false).front.coverage;
}

void getDepths(ref Table t,string prefix){
    BamRegion[] regions;
    auto bam = new BamReader(buildPath(prefix,t.samples.front~".bam"));
    foreach(rec;t.records){
        BamRegion reg;
        reg.ref_id=rec.chr;
        reg.start=rec.pos;
        reg.end=rec.pos+1;
        regions~=reg;
    }
	foreach(j,sample;t.samples){
		bam = new BamReader(buildPath(prefix,sample~".bam"));
        auto reads=bam.getReadsOverlapping(regions);
        auto nr=t.records.sort.array;
		foreach(i,rec;enumerate(regions.sort)){
            //reads.map!(x=>x.ref_id).uniq.writeln;
            while(rec.ref_id!=reads.front.ref_id){
                reads.popFront;
            }
            //writeln(nr[i],rec,takeUntil!"a.ref_id!=b"(reads,rec.ref_id).makePileup(false,rec.start,rec.start,false).front.coverage);
			t.matrix[i][j]=takeUntil!"a.ref_id!=b"(reads,rec.ref_id).makePileup(false,rec.start,rec.start,false).front.coverage;
		} 
        writeln(sample);
	}
}




